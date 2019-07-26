{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.MixedModel where

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

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH

import qualified Numeric.LinearAlgebra         as LA
--import qualified Numeric.NLOPT                 as NL



import qualified Data.List                     as L
import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

type LevelSpec = (Int, Bool, Maybe (VB.Vector Bool)) -- level sizes and effects

data DevianceType = ML | REML deriving (Show, Eq)

effectsForLevel :: LevelSpec -> Int
effectsForLevel (_, b, vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForLevel #-}

colsForLevel :: LevelSpec -> Int
colsForLevel l@(qL, _, _) = qL * effectsForLevel l
{-# INLINABLE colsForLevel #-}

type Levels = VB.Vector LevelSpec

-- classify row into its levels
-- the vector has a level number for each level
type RowClassifier = Int -> VB.Vector Int

-- for group k, what group effects are we modeling?
-- first Bool is for intercept, then (optional) vector,
-- of same length as X has columns, to indicate which
-- predictors in X get random slopes
-- NB:  If X has a constant column, there is redundancy
-- here between the Bool in the tuple and the first bool
-- in the vector.
--type LevelEffectsSpec = [(Bool, Maybe (LA.Vector Bool))]

type SemC r = (MonadIO (P.Sem r), {- P.Member (P.Error SomeException) r,-} P.Member (P.Error T.Text) r)
runPIRLS_M
  :: P.Sem '[{-P.Error SomeException,-}
             P.Error T.Text, P.Lift IO] a
  -> IO (Either T.Text a)
runPIRLS_M = P.runM . P.runError {-. P.runErrorAsAnother (T.pack . show)-}

makeZ
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double)
  => LA.Matrix a
  -> Levels
  -> RowClassifier
  -> SLA.SpMatrix a
makeZ mX levels rc =
  let
    (nO, nP) = LA.size mX
    k        = VB.length levels -- number of levels
    levelSize n = let (s, _, _) = levels VB.! n in s

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
    zEntriesForLevel startingCol level (qL, b, vbM) =
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

asDense sp =
  LA.matrix (snd $ SLA.dimSM sp) $ fmap (\(_, _, x) -> x) $ SLA.toDenseListSM sp

asDenseV :: SLA.SpVector Double -> LA.Vector Double
asDenseV = VS.convert . SLA.toVectorDense


toSparseMatrix
  :: (LA.Container LA.Vector a, RealFrac a) => LA.Matrix a -> SLA.SpMatrix a
toSparseMatrix x
  = let
      (r, c)     = LA.size x
      colIndices = take c $ L.iterate (+ 1) 0
      rowIndices = take r $ L.iterate (+ 1) 0
      indexRow rI =
        fmap (\(cI, val) -> (rI, cI, val)) . zip colIndices . LA.toList
      indexedRows =
        concat . fmap (\(rI, rV) -> indexRow rI rV) . zip rowIndices . LA.toRows
    in
      SLA.fromListSM (r, c) $ indexedRows x

toSparseVector
  :: (LA.Container LA.Vector a, RealFrac a) => LA.Vector a -> SLA.SpVector a
toSparseVector v = SLA.fromListDenseSV (LA.size v) $ LA.toList v

-- TODO: Add levels checks
checkProblem
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double, SemC r)
  => LA.Matrix a -- ^ X
  -> LA.Vector a -- ^ y
  -> SLA.SpMatrix a -- ^ Z
  -> P.Sem r ()
checkProblem mX vY smZ = do
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
makeSTF
  :: (Num a, VS.Storable a)
  => Levels
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a))
makeSTF levels
  = let
      f
        :: Num a
        => LevelSpec
        -> ([a] -> [(Int, Int, a)], [a] -> [(Int, Int, a)], Int, Int)
      f l@(n, _, _) =
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

xTxPlusI :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTxPlusI smX =
  (SLA.transposeSM smX) SLA.## smX SLA.^+^ (SLA.eye $ SLA.ncols smX)

logDetTriangularSM :: RealFloat a => SLA.SpMatrix a -> a
logDetTriangularSM smX =
  let f (_, _, x) = log x
      diag (r, c, _) = (r == c)
  in  FL.fold (FL.premap f FL.sum) $ L.filter diag $ SLA.toListSM smX

profiledDeviance2
  :: CH.ForeignPtr CH.Common
  -> CH.ForeignPtr CH.Factor -- ^ precomputed pattern work on Z
  -> SLA.SpMatrix Double -- ^ permutation matrix from above
  -> DevianceType
  -> (  LA.Vector Double
     -> (SLA.SpMatrix Double, SLA.SpMatrix Double)
     ) -- ^ make S and T from theta
  -> LA.Matrix Double -- ^ X
  -> LA.Vector Double -- ^ y
  -> SLA.SpMatrix Double -- ^ Z
  -> LA.Vector Double -- ^ theta
  -> IO (Double, LA.Vector Double, LA.Vector Double) -- ^ (pd, beta, b) 
profiledDeviance2 fpc fpf smP dt mkST mX vY smZ vTh = do
  let (smS, smT) = mkST vTh
      smZS       = smZ SLA.## (smT SLA.## smS)
      smZSt      = SLA.transpose smZS
      n          = LA.size vY
      (_, p)     = LA.size mX
      (_, q)     = SLA.dim smZ
  -- Cholesky factorize to get L_theta *and* update factor for solving with it
  CH.spMatrixFactorize fpc fpf CH.SquareSymmetricLower $ xTxPlusI $ smZS
  -- compute Rzx
  let smZStX = smZSt SLA.## (toSparseMatrix mX)
  -- TODO: these can probably be combined as a single function in CholmodExtras, saving some copying of data
  smPZStX <- CH.solveSparse fpc fpf CH.CHOLMOD_P smZStX -- P(Z*)'X
  smRzx   <- CH.solveSparse fpc fpf CH.CHOLMOD_LD smPZStX -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
  smLth   <- CH.choleskyFactorSM fpf fpc -- this has to happen after the solves because it unfactors the factor...
  -- compute Rx
  let xTxMinusRzxTRzx =
        (LA.tr mX) LA.<> mX - (asDense $ (SLA.transposeSM smRzx) SLA.## smRzx)
  let smRx = toSparseMatrix $ LA.chol $ LA.trustSym $ xTxMinusRzxTRzx
      smUT = (SLA.transpose smLth -||- smRzx) -=- (SLA.zeroSM p q -||- smRx)
      smLT = SLA.transpose smUT
      svBu = (smP SLA.## smZSt) SLA.#> (toSparseVector vY)
      svBl = toSparseVector $ (LA.tr mX) LA.#> vY
  let svB =
        SLA.fromListSV (q + p)
          $  (SLA.toListSV svBu)
          ++ (fmap (\(i, x) -> (i + q, x)) (SLA.toListSV svBl))
  svX :: SLA.SpVector Double <- SLA.luSolve smLT smUT svB
  let svPu          = SLA.takeSV q svX
      svu           = (SLA.transpose smP) SLA.#> svPu -- I could also do this via a cholmod solve
      svb           = (smT SLA.## smS) SLA.#> svu
      vBeta         = asDenseV $ SLA.takeSV p $ SLA.dropSV q svX
      vDev          = vY - (mX LA.#> vBeta) - (asDenseV $ smZS SLA.#> svu)
      rTheta2       = (vDev LA.<.> vDev) + (svu SLA.<.> svu)
      logLth        = logDetTriangularSM smLth
      (logDet, dof) = case dt of
        ML   -> (realToFrac n, logLth)
        REML -> (realToFrac (n - p), logLth + (logDetTriangularSM smRx))
      pd = 2 * logDet + dof * (1 + (2 * pi * rTheta2 / dof))
  return (pd, vBeta, asDenseV $ svb)


