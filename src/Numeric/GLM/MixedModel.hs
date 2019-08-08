{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.GLM.MixedModel where

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH
import qualified Numeric.SparseDenseConversions
                                               as SD
import qualified Numeric.GLM.Types             as GLM

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P

import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA
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
data GroupFitSpec = GroupFitSpec { nCategories :: Int
                                 , groupIntercept :: Bool
                                 , groupSlopes :: Maybe (VB.Vector Bool)
                                 } deriving (Show, Eq)

type GroupFitSpecs = VB.Vector GroupFitSpec

data MixedModel = MixedModel RegressionModel GroupFitSpecs

type RandomEffectModelMatrix = SLA.SpMatrix Double
type CovarianceVector = LA.Vector Double
type Lambda = SLA.SpMatrix Double

data RandomEffectCalculated = RandomEffectCalculated RandomEffectModelMatrix (CovarianceVector -> Lambda)

type PermutationMatrix = SLA.SpMatrix Double -- should this be (SLA.SpMatrix Int) ??

data DevianceType = ML | REML deriving (Show, Eq)

effectsForGroup :: GroupFitSpec -> Int
effectsForGroup (GroupFitSpec _ b vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForGroup #-}

colsForGroup :: GroupFitSpec -> Int
colsForGroup l = nCategories l * effectsForGroup l
{-# INLINABLE colsForGroup #-}

-- classify row into its levels
-- the vector has a category number for each level
-- e.g., rc :: RowClassifier
-- rc (rowNumber :: Int) == [cat1 :: Int, cat2 :: Int]
-- catN :: RowClass
type RowClassifier = Int -> VB.Vector Int

categoryNumber :: RowClassifier -> Int -> Int -> Int
categoryNumber rowClassifier rowIndex levelIndex =
  rowClassifier rowIndex VB.! levelIndex

data ParameterEstimates =
  ParameterEstimates
  { pEstimate :: LA.Vector Double
  , pCovariance :: LA.Matrix Double
  }

type FixedParameterEstimates = ParameterEstimates
type GroupParameterEstimates = VB.Vector FixedParameterEstimates

type SemC r = (MonadIO (P.Sem r), P.Member (P.Error T.Text) r)
runPIRLS_M :: P.Sem '[P.Error T.Text, P.Lift IO] a -> IO (Either T.Text a)
runPIRLS_M = P.runM . P.runError {-. P.runErrorAsAnother (T.pack . show)-}

setCovarianceVector :: GroupFitSpecs -> Double -> Double -> LA.Vector Double
setCovarianceVector groupFS diag offDiag = FL.fold fld groupFS
 where
  fld = FL.Fold
    (\bs l ->
      let e       = effectsForGroup l
          entries = e * (e + 1) `div` 2
          col n = n `div` e
          row n = col n + n `mod` e
          set n = if (row n == col n) then diag else offDiag
      in  bs ++ fmap set (take entries $ iterate (+ 1) 0)
    )
    []
    LA.fromList


{-
Z is the random effects model matrix
like X, the fixed effect model matrix,
it has a row for each observation.
Z has a column for each random effect:
A column of ones for a random intercept and
a column with the predictor from X for a random
slope.
-}
makeZ
  :: FixedPredictors
  -> GroupFitSpecs
  -> RowClassifier
  -> RandomEffectModelMatrix
makeZ mX groupFSs rc =
  let
    (nO, nP) = LA.size mX
    k        = VB.length groupFSs -- number of levels
    groupSize n = nCategories (groupFSs VB.! n)
    q = FL.fold FL.sum $ fmap colsForGroup groupFSs -- total number of columns in Z
    predictor rowIndex fixedEffectIndex =
      mX `LA.atIndex` (rowIndex, fixedEffectIndex)
    -- construct Z for level as a fold over rows
    entries :: Int -> Int -> GroupFitSpec -> Int -> [(Int, Int, Double)]
    entries colOffset groupNumber groupFS rowIndex
      = let
          entryOffset =
            colOffset
              + ( effectsForGroup groupFS
                * categoryNumber rc rowIndex groupNumber
                )
          (intercept, slopeOffset) = if groupIntercept groupFS
            then ([(rowIndex, entryOffset, 1)], entryOffset + 1)
            else ([], entryOffset)
          slopeF :: FL.Fold (Bool, Int) [(Int, Int, Double)]
          slopeF = FL.Fold
            (\(mes, ziCol) (b, fixedEffectIndex) -> if b
              then
                ( (rowIndex, ziCol, predictor rowIndex fixedEffectIndex) : mes
                , ziCol + 1
                )
              else (mes, ziCol)
            )
            ([], slopeOffset)
            fst
          slopes = case groupSlopes groupFS of
            Just slopeV -> FL.fold slopeF $ VB.zip slopeV (VB.generate nP id)
            Nothing     -> []
        in
          intercept ++ slopes
    ziFold groupNumber groupOffset groupFS =
      let q = colsForGroup groupFS
      in  FL.Fold (\mes n -> mes ++ entries groupOffset groupNumber groupFS n)
                  []
                  id
    groupNumbers = VB.generate (VB.length groupFSs) id
    groupOffsets = VB.prescanl' (+) 0 $ fmap colsForGroup groupFSs
    zis          = concat $ FL.fold
      ( sequenceA
      $ fmap (\(gNumber, gOffset, gSpec) -> ziFold gNumber gOffset gSpec)
      $ VB.zip3 groupNumbers groupOffsets groupFSs
      )
      [0 .. (nO - 1)]
  in
    SLA.fromListSM (nO, q) zis

-- TODO: Add levels checks
-- TODO: Add Lambda checks
checkProblem :: SemC r => MixedModel -> RandomEffectCalculated -> P.Sem r ()
checkProblem (MixedModel (RegressionModel mX vY) _) (RandomEffectCalculated smZ _)
  = do
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


-- For each level, i, with p_i random effects,
-- Lambda is built from lower-triangular blocks
-- which are p_i x p_i.  Theta contains the values
-- for each block, in col major order
makeLambda :: GroupFitSpecs -> (CovarianceVector -> Lambda)
makeLambda groupFSs =
  let blockF vTH = FL.Fold
        (\(bs, vx) l' -> ((l', vx) : bs, VS.drop (effectsForGroup l') vx))
        ([], vTH)
        (reverse . fst)
      blockData vTh = FL.fold (blockF vTh) groupFSs
      templateBlock l vTh =
        let pl  = effectsForGroup l
            lTh = VS.toList vTh
            lts vx = fmap (\((r, c), v) -> (r, c, v)) $ zip
              ([ (r, c) | c <- [0 .. (pl - 1)], r <- [c .. (pl - 1)] ])
              vx
        in  SLA.fromListSM (pl, pl) $ lts lTh
      perGroup (l, vx) = replicate (nCategories l) $ templateBlock l vx
      allDiags vTh = concat $ fmap perGroup $ blockData vTh
  in  (\vTh -> SLA.fromBlocksDiag (allDiags vTh))

xTxPlusI :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTxPlusI smX = (smX SLA.#~^# smX) SLA.^+^ (SLA.eye $ SLA.ncols smX)

logDetTriangularSM :: RealFloat a => SLA.SpMatrix a -> a
logDetTriangularSM smX =
  let f (_, _, x) = log x
      diag (r, c, _) = (r == c)
  in  FL.fold (FL.premap f FL.sum) $ L.filter diag $ SLA.toListSM smX

logDetTriangularM :: (RealFloat a, LA.Container VS.Vector a) => LA.Matrix a -> a
logDetTriangularM mX =
  let f n = log $ mX `LA.atIndex` (n, n)
      (rows, _) = LA.size mX
  in  FL.fold (FL.premap f FL.sum) [0 .. (rows - 1)]

type CholmodFactor = (CH.ForeignPtr CH.Common -- ^ pre-allocated CHOMOD common space
                     , CH.ForeignPtr CH.Factor -- ^ precomputed pattern work on Z
                     , PermutationMatrix -- ^ permutation matrix from above
                     )

cholmodAnalyzeProblem
  :: SemC r => RandomEffectCalculated -> P.Sem r CholmodFactor
cholmodAnalyzeProblem (RandomEffectCalculated smZ _) = do
  cholmodC <- liftIO CH.allocCommon
  liftIO $ CH.startC cholmodC
  liftIO $ CH.setFinalLL 1 cholmodC -- I don't quite get this.  We should be able to solve LDx = b either way??
  (cholmodF, smP) <-
    liftIO $ CH.spMatrixAnalyzeWP cholmodC CH.SquareSymmetricLower $ xTxPlusI
      smZ
--  liftIO $ CH.printFactor cholmodF "After analyze, before factorize" cholmodC
--  liftIO $ putStrLn "smP=" >> (LA.disp 1 $ SD.toDenseMatrix smP)
  return (cholmodC, cholmodF, smP)

data ProfiledDevianceVerbosity = PDVNone | PDVSimple | PDVAll deriving (Show, Eq)

profiledDeviance
  :: ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> MixedModel
  -> RandomEffectCalculated
  -> CovarianceVector -- ^ theta
  -> IO
       ( Double
       , LA.Vector Double
       , SLA.SpVector Double
       , SLA.SpVector Double
       ) -- ^ (pd, beta, u, b) 
profiledDeviance verbosity cf dt mm@(MixedModel (RegressionModel mX vY) _) reCalc@(RandomEffectCalculated smZ mkLambda) vTh
  = do
    let n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        lambda = mkLambda vTh
        smZS   = smZ SLA.## lambda --(smT SLA.## smS)
        smZSt  = SLA.transpose smZS
    CholeskySolutions smLth smRzx mRx svBeta svu <- cholmodCholeskySolutions
      cf
      mm
      reCalc
      vTh
    when (verbosity == PDVAll) $ do
      putStrLn "Z"
      LA.disp 2 $ SD.toDenseMatrix smZ
      putStrLn "Lambda"
      LA.disp 2 $ SD.toDenseMatrix lambda
      putStrLn "Lth"
      LA.disp 2 $ SD.toDenseMatrix smLth
      putStrLn "Rzx"
      LA.disp 2 $ SD.toDenseMatrix smRzx
      putStrLn "Rx"
      LA.disp 2 mRx
    let svb     = lambda SLA.#> svu -- (smT SLA.## smS) SLA.#> svu
        vBeta   = SD.toDenseVector svBeta
        vDev    = vY - (mX LA.#> vBeta) - (SD.toDenseVector $ smZS SLA.#> svu)
        rTheta2 = (vDev LA.<.> vDev) + (svu SLA.<.> svu)
    let logLth        = logDetTriangularSM smLth
        (dof, logDet) = case dt of
          ML   -> (realToFrac n, logLth)
          REML -> (realToFrac (n - p), logLth + (logDetTriangularM mRx))
        pd = (2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof)))
    when (verbosity == PDVAll) $ do
      let (_, _, smP) = cf
      putStrLn $ "smP="
      LA.disp 1 $ SD.toDenseMatrix smP
      putStrLn $ "rTheta^2=" ++ show rTheta2
      putStrLn $ "2 * logLth=" ++ show (2 * logLth)
      putStrLn $ "2 * logDet=" ++ show (2 * logDet)
    when (verbosity == PDVAll || verbosity == PDVSimple) $ do
      putStrLn $ "pd(th=" ++ (show vTh) ++ ") = " ++ show pd
    return (pd, vBeta, svu, svb)

upperTriangular :: Int -> Int -> a -> Bool
upperTriangular r c _ = (r <= c)

lowerTriangular :: Int -> Int -> a -> Bool
lowerTriangular r c _ = (r >= c)

data CholeskySolutions = CholeskySolutions { lTheta :: SLA.SpMatrix Double
                                           , rZX :: SLA.SpMatrix Double
                                           , rX :: LA.Matrix Double
                                           , svBeta :: SLA.SpVector Double
                                           , svu :: SLA.SpVector Double
                                           } deriving (Show)

cholmodCholeskySolutions
  :: CholmodFactor
  -> MixedModel
  -> RandomEffectCalculated
  -> CovarianceVector
  -> IO CholeskySolutions
cholmodCholeskySolutions cholmodFactor mixedModel randomEffCalcs vTh = do
  let (cholmodC, cholmodF, _) = cholmodFactor
      (MixedModel             (RegressionModel mX vY) levels  ) = mixedModel
      (RandomEffectCalculated smZ                     mkLambda) = randomEffCalcs
      lambda                  = mkLambda vTh
      smZS                    = smZ SLA.## lambda
      smZSt                   = SLA.transpose smZS
      n                       = LA.size vY
      (_, p)                  = LA.size mX
      (_, q)                  = SLA.dim smZ
  -- Cholesky factorize to get L_theta *and* update factor for solving with it
  CH.spMatrixFactorize cholmodC cholmodF CH.SquareSymmetricLower
    $ SLA.filterSM lowerTriangular
    $ xTxPlusI
    $ smZS
    -- compute Rzx

  let svZSty = smZSt SLA.#> SD.toSparseVector vY
  smPZSty <- CH.solveSparse cholmodC
                            cholmodF
                            CH.CHOLMOD_P
                            (SD.svColumnToSM svZSty)
  svCu <- SLA.toSV <$> CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPZSty
  let smZStX = smZSt SLA.#~# (SD.toSparseMatrix mX)
  smPZStX <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smZStX -- P(Z*)'X
  smRzx   <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPZStX -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
  --  let smRxTRx = (SD.toSparseMatrix $ LA.tr mX LA.<> mX) SLA.^-^ (SLA.transposeSM smRzx SLA.## smRzx)
  --      beta = SLA.luSolve (SD.toSparseVector (mX LA.#> vY) SLA.^-^ (smRzx SLA.#> svCu)
  -- compute Rx
  let xTxMinusRzxTRzx =
        (LA.tr mX) LA.<> mX - (SD.toDenseMatrix $ smRzx SLA.#^# smRzx)
      mRx              = LA.chol $ LA.trustSym $ xTxMinusRzxTRzx
      vXty             = (LA.tr mX) LA.#> vY
      vRzxtCu          = SD.toDenseVector $ (SLA.transposeSM smRzx) SLA.#> svCu
      vXtyMinusRzxtCu  = vXty - vRzxtCu
      betaSols         = LA.cholSolve mRx (LA.asColumn vXtyMinusRzxtCu)
      vBeta            = head $ LA.toColumns $ betaSols
      svBeta           = SD.toSparseVector vBeta
      svRzxBeta        = smRzx SLA.#> svBeta
--      smRzxBeta        = SD.svColumnToSM svRzxBeta
      smCuMinusRzxBeta = SD.svColumnToSM (svCu SLA.^-^ svRzxBeta)
  smPu  <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Lt smCuMinusRzxBeta
  svu   <- SLA.toSV <$> (CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Pt $ smPu)
  -- NB: This has to happen after the solves because it unfactors the factor...
  -- NB: It does *not* undo the analysis
  smLth <- CH.choleskyFactorSM cholmodF cholmodC
  return $ CholeskySolutions smLth smRzx mRx svBeta svu


