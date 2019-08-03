{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.MixedModel where

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH
import qualified Numeric.SparseDenseConversions as SD

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P

import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Numeric.LinearAlgebra.Class   as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import           Numeric.LinearAlgebra.Sparse   ( (##)
                                                , (#^#)
                                                , (<#)
                                                , (#>)
                                                , (-=-)
                                                , (-||-)
                                                )


import qualified Numeric.LinearAlgebra         as LA
--import qualified Numeric.NLOPT                 as NL



import qualified Data.List                     as L
import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

type FixedPredictors = LA.Matrix Double
type Observations = LA.Vector Double

data RegressionModel = RegressionModel FixedPredictors Observations

-- for group k, what group effects are we modeling?
-- first Bool is for intercept, then (optional) vector,
-- of same length as X has columns, to indicate which
-- predictors in X get random slopes
-- NB:  If X has a constant column, there is redundancy
-- here between the Bool in the tuple and the first bool
-- in the vector.
data LevelSpec = LevelSpec { nCategories :: Int
                           , groupIntercept :: Bool
                           , groupSlopes :: Maybe (VB.Vector Bool)
                           } deriving (Show, Eq)

type Levels = VB.Vector LevelSpec

data MixedModel = MixedModel RegressionModel Levels

type RandomEffectModelMatrix = SLA.SpMatrix Double
type CovarianceVector = LA.Vector Double
type Lambda = SLA.SpMatrix Double

data RandomEffectCalculated = RandomEffectCalculated RandomEffectModelMatrix (CovarianceVector -> Lambda)

type PermutationMatrix = SLA.SpMatrix Double -- should this be (SLA.SpMatrix Int) ??

data DevianceType = ML | REML deriving (Show, Eq)

effectsForLevel :: LevelSpec -> Int
effectsForLevel (LevelSpec _ b vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForLevel #-}

colsForLevel :: LevelSpec -> Int
colsForLevel l = nCategories l * effectsForLevel l
{-# INLINABLE colsForLevel #-}

-- classify row into its levels
-- the vector has a level number for each level
type RowClassifier = Int -> VB.Vector Int

--type LevelEffectsSpec = [(Bool, Maybe (LA.Vector Bool))]

type SemC r = (MonadIO (P.Sem r), P.Member (P.Error T.Text) r)
runPIRLS_M
  :: P.Sem '[{-P.Error SomeException,-}
             P.Error T.Text, P.Lift IO] a
  -> IO (Either T.Text a)
runPIRLS_M = P.runM . P.runError {-. P.runErrorAsAnother (T.pack . show)-}

setCovarianceVector :: Levels -> Double -> Double -> LA.Vector Double
setCovarianceVector levels diag offDiag =  FL.fold fld levels
 where
  fld = FL.Fold
    (\bs l ->
      let e = effectsForLevel l
      in  bs ++ L.replicate e diag ++ replicate (e * (e - 1) `div` 2) offDiag
    )
    []
    LA.fromList


makeZ
  :: FixedPredictors
  -> Levels
  -> RowClassifier
  -> RandomEffectModelMatrix
makeZ mX levels rc =
  let
    (nO, nP) = LA.size mX
    k        = VB.length levels -- number of levels
    levelSize n = nCategories (levels VB.! n)
    q          = FL.fold FL.sum $ fmap colsForLevel levels -- total number of columns in Z
    -- build list of items to put in Z, level by level
    obsIndices = Seq.fromFunction nO id
    intercept startingCol level =
      let spIndex obsIndex = (obsIndex, rc obsIndex VB.! level)
          spEntry obsIndex =
            let spi = spIndex obsIndex in (fst spi, startingCol + snd spi, 1)
      in  fmap spEntry obsIndices
    slope startingCol level predIndex =
      let spIndex obsIndex = (obsIndex, rc obsIndex VB.! level)
          spValue obsIndex = mX `LA.atIndex` (obsIndex, predIndex)
          spEntry obsIndex =
            let spi = spIndex obsIndex
            in  (fst spi, startingCol + snd spi, spValue obsIndex)
      in  fmap spEntry obsIndices
    zEntriesForSlope startingCol level vb =
      let predIndices = fmap fst $ VB.filter snd $ VB.zip
            (VB.generate (VB.length vb) id)
            vb
          startingCols = VB.generate
            (VB.length predIndices)
            (\n -> startingCol + n * (levelSize level))
          newStart = startingCol + (VB.length predIndices * (levelSize level))
          entriesF (pi, sc) = slope sc level pi
      in  ( FL.fold FL.mconcat $ fmap entriesF $ VB.zip predIndices startingCols
          , newStart
          )
    zEntriesForLevel startingCol level (LevelSpec qL b vbM) =
      let (interceptZs, newStartI) = if b
            then (intercept startingCol level, startingCol + qL)
            else (Seq.empty, startingCol)
          (slopeZs, newStartS) = case vbM of
            Nothing -> (Seq.empty, newStartI)
            Just vb -> zEntriesForSlope newStartI level vb
      in  (interceptZs <> slopeZs, newStartS)
    zFold = FL.Fold
      (\(zs, sc, l) x ->
        let (newZs, newStart) = zEntriesForLevel sc l x
        in  (zs <> newZs, newStart, l + 1)
      )
      (Seq.empty, 0, 0)
      id
    (zEntries, numCols, _) = FL.fold zFold levels
  in
    SLA.fromListSM (nO, q) zEntries

-- TODO: Add levels checks
-- TODO: Add Lambda checks
checkProblem :: SemC r => MixedModel
             -> RandomEffectCalculated
             -> P.Sem r ()
checkProblem (MixedModel (RegressionModel mX vY) _) (RandomEffectCalculated smZ _)  = do
  let (n, p)     = LA.size mX
      yRows      = LA.size vY
      (zRows, q) = SLA.dim smZ
  when (zRows /= n)
    $  P.throw @T.Text
    $  (T.pack $ show n)
    <> " cols in X but "
    <> (T.pack $ show zRows)
    <> " cols in Z"
  when (yRows /= n)
    $  P.throw @T.Text
    $  (T.pack $ show n)
    <> " cols in X but "
    <> (T.pack $ show yRows)
    <> " entries in Y"
  return ()


-- S is the diagonal matrix of covariances in theta
-- T is unit-lower-triangular of off-diagonal covariances
-- this function produces a function mapping the vector theta to the matrix Lambda = ST
makeSTF
  :: Levels
  -> (CovarianceVector -> (SLA.SpMatrix Double, SLA.SpMatrix Double))
makeSTF levels
  = let
      f
        :: Num a
        => LevelSpec
        -> ([a] -> [(Int, Int, a)], [a] -> [(Int, Int, a)], Int, Int)
      f l@(LevelSpec n _ _) =
        let
          e = effectsForLevel l
          s' th =
            fmap (\(x, v) -> (x, x, v)) $ zip (take e $ L.iterate (+ 1) 0) th
          tlt' th = fmap (\((r, c), v) -> (r, c, v)) $ zip
            ([ (r, c) | c <- [0 .. (e - 1)], r <- [(c + 1) .. (e - 1)] ])
            th
          t' th =
            tlt' th ++ fmap (\x -> (x, x, 1)) (L.take e $ L.iterate (+ 1) 0)
        in
          (s', t', e, n)
      makers = fmap f levels
      diagOffset o = fmap (\(r, c, v) -> (r + o, c + o, v))
      diagCopiesFrom i e n ivs =
        mconcat $ fmap (flip diagOffset ivs) $ L.take n $ L.iterate (+ e) i
    in
      \thV ->
        let thL = LA.toList thV
            fld = FL.Fold
              (\(sL, tL, thL', offset) (s', t', e, n) ->
                ( sL ++ diagCopiesFrom offset e n (s' thL')
                , tL ++ diagCopiesFrom offset e n (t' $ L.drop e thL')
                , L.drop (e + (e * (e - 1) `div` 2)) thL'
                , offset + (e * n)
                )
              )
              ([], [], thL, 0)
              id
            (s, t, _, qT) = FL.fold fld makers
        in  (SLA.fromListSM (qT, qT) s, SLA.fromListSM (qT, qT) t)

-- For each level, i, with p_i random effects,
-- Lambda is built from lower-triangular blocks
-- which are p_i x p_i.  Theta contains the values
-- for each block, diagonals first, then off diagonals,
-- in col major order.
makeLambda
  :: Levels
  -> (CovarianceVector -> Lambda)
makeLambda levels =
  let blockF vTH =
        FL.Fold (\(bs, vx) l' -> ((l', vx) : bs,  VS.drop (effectsForLevel l') vx))
        ([], vTH)
        (reverse . fst)
      blockData vTh = FL.fold (blockF vTh) levels
      templateBlock l vTh = 
        let pl = effectsForLevel l
            lTh = VS.toList vTh
            diags vx = fmap (\(x, v) -> (x, x, v)) $ zip (take pl $ L.iterate (+ 1) 0) vx
            lts  vx =  fmap (\((r, c), v) -> (r, c, v)) $ zip
              ([ (r, c) | c <- [0 .. (pl - 1)], r <- [(c + 1) .. (pl - 1)] ])
              vx
        in SLA.fromListSM (pl, pl) $ diags lTh ++ lts (drop pl lTh)
      perLevel (l, vx) = replicate (nCategories l) $ templateBlock l vx
      allDiags vTh = concat $ fmap perLevel $ blockData vTh
  in (\vTh -> SLA.fromBlocksDiag (allDiags vTh))

xTxPlusI :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTxPlusI smX =
  (SLA.transposeSM smX SLA.## smX) SLA.^+^ (SLA.eye $ SLA.ncols smX)

logDetTriangularSM :: RealFloat a => SLA.SpMatrix a -> a
logDetTriangularSM smX =
  let f (_, _, x) = log x
      diag (r, c, _) = (r == c)
  in  FL.fold (FL.premap f FL.sum) $ L.filter diag $ SLA.toListSM smX

type CholmodFactor = (CH.ForeignPtr CH.Common -- ^ pre-allocated CHOMOD common space
                     , CH.ForeignPtr CH.Factor -- ^ precomputed pattern work on Z
                     , PermutationMatrix -- ^ permutation matrix from above
                     )

profiledDeviance2
  :: CholmodFactor
  -> DevianceType
  -> MixedModel
  -> RandomEffectCalculated
  -> CovarianceVector -- ^ theta
  -> IO (Double, LA.Vector Double, LA.Vector Double , LA.Vector Double) -- ^ (pd, beta, u, b) 
profiledDeviance2 (fpc, fpf, smP)
  dt
  (MixedModel (RegressionModel mX vY) _)
  (RandomEffectCalculated smZ mkLambda)
  vTh =
  do
    let -- (smS, smT) = mkST vTh
      lambda = mkLambda vTh
      smZS       = smZ SLA.## lambda --(smT SLA.## smS)
      smZSt      = SLA.transpose smZS
      n          = LA.size vY
      (_, p)     = LA.size mX
      (_, q)     = SLA.dim smZ
      upperTriangular r c _ = (r <= c)
      lowerTriangular r c _ = (r >= c)
    -- Cholesky factorize to get L_theta *and* update factor for solving with it
    CH.spMatrixFactorize fpc fpf CH.SquareSymmetricLower
      $ SLA.filterSM lowerTriangular
      $ xTxPlusI
      $ smZS
    -- compute Rzx
    let smZStX = smZSt SLA.## (SD.toSparseMatrix mX)
    -- TODO: these can probably be combined as a single function in CholmodExtras, saving some copying of data
    smPZStX <- CH.solveSparse fpc fpf CH.CHOLMOD_P smZStX -- P(Z*)'X
    smRzx   <- CH.solveSparse fpc fpf CH.CHOLMOD_LD smPZStX -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
    smLth   <- CH.choleskyFactorSM fpf fpc -- this has to happen after the solves because it unfactors the factor...
  -- compute Rx
    let xTxMinusRzxTRzx =
          (LA.tr mX) LA.<> mX - (SD.toDenseMatrix $ (SLA.transposeSM smRzx) SLA.## smRzx)
    let smRx = SD.toSparseMatrix $ LA.chol $ LA.trustSym $ xTxMinusRzxTRzx
        smUT =
          SLA.filterSM upperTriangular
          $   (SLA.transpose smLth -||- smRzx)
          -=- (SLA.zeroSM p q -||- smRx)
    let smLT = SLA.transpose smUT
    let svBu = (smP SLA.## smZSt) SLA.#> (SD.toSparseVector vY)
        svBl = SD.toSparseVector $ (LA.tr mX) LA.#> vY
    let svB =
          SLA.fromListSV (q + p)
          $  (SLA.toListSV svBu)
          ++ (fmap (\(i, x) -> (i + q, x)) (SLA.toListSV svBl))
    svX :: SLA.SpVector Double <- SLA.luSolve smLT smUT svB
    let svPu    = SLA.takeSV q svX
        svu     = (SLA.transpose smP) SLA.#> svPu -- I could also do this via a cholmod solve
        svb     = lambda SLA.#> svu -- (smT SLA.## smS) SLA.#> svu
        vBeta   = SD.toDenseVector $ SLA.takeSV p $ SLA.dropSV q svX
        vDev    = vY - (mX LA.#> vBeta) - (SD.toDenseVector $ smZ SLA.#> svb)
        rTheta2 = (vDev LA.<.> vDev) + (svu SLA.<.> svu)
    let logLth        = logDetTriangularSM smLth
        (dof, logDet) = case dt of
          ML   -> (realToFrac n, logLth)
          REML -> (realToFrac (n - p), logLth + (logDetTriangularSM smRx))
        pd = (2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof)))
    return (pd, vBeta, SD.toDenseVector $ svu, SD.toDenseVector $ svb)


