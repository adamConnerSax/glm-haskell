{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE DeriveDataTypeable  #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.MixedModel where

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH
import qualified Numeric.SparseDenseConversions
                                               as SD
import qualified Numeric.GLM.Types             as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Data.IndexedSet               as IS

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Knit.Effect.Logger            as P

--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import qualified Control.Exception             as X
import           Control.Monad                  ( when, join )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Data.Array                    as A
import           Data.Function                  ((&))
import qualified Data.List                     as L
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA
import           Data.Typeable                  (Typeable)

import           Numeric.LinearAlgebra.Sparse   ( (#>) )

import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL
import           System.IO.Unsafe               ( unsafePerformIO )



import qualified Data.Map                      as M
import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

import System.IO (hSetBuffering, stdout, BufferMode(..))

data GLMError = OtherGLMError T.Text | NonGLMError T.Text deriving (Show, Typeable)
instance X.Exception GLMError

type Effects r = (P.Member (P.Error GLMError) r, P.LogWithPrefixesLE r)
type EffectsIO r = ( Effects r
                   , P.Member (P.Embed IO) r)
                   

type GLMEffects a = P.Sem '[P.Logger P.LogEntry, P.PrefixLog, P.Error GLMError, P.Embed IO, P.Final IO] a

runEffectsIO :: GLMEffects a -> IO (Either GLMError a)
runEffectsIO action = action 
  & P.filteredLogEntriesToIO P.logAll
  & P.errorToIOFinal
  & P.embedToFinal
  & P.runFinal


unsafePerformEffects :: GLMEffects a -> a
unsafePerformEffects action = action
  & runEffectsIO
  & fmap (either X.throwIO return)
  & join
  & unsafePerformIO


type FixedPredictors = LA.Matrix Double
type Observations = LA.Vector Double

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

makeGroupFitSpec
  :: (Ord b, Enum b, Bounded b, Effects r)
  => Int -- ^ number of items in group
  -> GLM.FixedEffects b
  -> GLM.IndexedEffectSet b
  -> P.Sem r GroupFitSpec
makeGroupFitSpec n fixedEffects indexedGroupEffects = do
  let indexedFixedEffects = GLM.indexedFixedEffectSet fixedEffects
  when (not $ IS.subset indexedGroupEffects indexedFixedEffects) $ P.throw
    $ OtherGLMError "group contains effects not in fixed effects in \"makeGroupFitSpec\""
  let modeled        = isJust . IS.index indexedGroupEffects
      groupIntercept = modeled GLM.Intercept
      modeledNI x = if x == GLM.Intercept then False else modeled x
      groupSlopesL = fmap modeledNI $ IS.members indexedFixedEffects
      groupSlopes  = if True `L.elem` groupSlopesL
        then Just (VB.fromList groupSlopesL)
        else Nothing
  return $ GroupFitSpec n groupIntercept groupSlopes

type FitSpecByGroup g = M.Map g GroupFitSpec --VB.Vector GroupFitSpec

fitSpecByGroup
  :: (Ord g, Show g, Ord b, Enum b, Bounded b, Show b, Effects r)
  => GLM.FixedEffects b
  -> GLM.EffectsByGroup g b
  -> GLM.RowClassifier g
  -> P.Sem r (FitSpecByGroup g)
fitSpecByGroup fixedEffects ebg rowClassifier = do
  let lookupGroupSize (grp, ie) =
        M.lookup grp (GLM.groupSizes rowClassifier) >>= return . (grp, , ie)
      makeSpec (grp, n, ige) =
        makeGroupFitSpec n fixedEffects ige >>= return . (grp, )
      lookupError =
        "Failed to find a group in the group effect set ("
          <> (T.pack $ show ebg)
          <> ") in  the rowClassifier ("
          <> (T.pack $ show rowClassifier)
  groupInfoList <-
    maybe (P.throw $ OtherGLMError lookupError) return $ traverse lookupGroupSize $ M.toList ebg
  fmap M.fromList $ traverse makeSpec groupInfoList

data RegressionModel b = RegressionModel (GLM.FixedEffects b) FixedPredictors Observations deriving (Show, Eq)
data MixedModel b g = MixedModel (RegressionModel b) (FitSpecByGroup g) deriving (Show, Eq)

data GeneralizedLinearMixedModel b g =
  LMM (MixedModel b g)
  | GLMM (MixedModel b g) WMatrix GLM.ObservationDistribution deriving (Show, Eq)

mixedModel :: GeneralizedLinearMixedModel b g -> MixedModel b g
mixedModel (LMM x     ) = x
mixedModel (GLMM x _ _) = x

generalized :: GeneralizedLinearMixedModel b g -> Bool
generalized (LMM _) = False
generalized (GLMM _ _ _) = True

type RandomEffectModelMatrix = SLA.SpMatrix Double
type CovarianceVec = LA.Vector Double
type LambdaMatrix = SLA.SpMatrix Double

type WMatrix = LA.Vector Double -- constant weights
type UMatrix = SLA.SpMatrix Double
type VMatrix = LA.Matrix Double
type LMatrix = SLA.SpMatrix Double -- Cholesky Z block
type RzxMatrix = SLA.SpMatrix Double
type RxxMatrix = LA.Matrix Double

type BetaVec = LA.Vector Double
type UVec = SLA.SpVector Double
type MuVec = LA.Vector Double
type EtaVec = LA.Vector Double

data BetaU = BetaU { vBeta :: LA.Vector Double, svU :: SLA.SpVector Double } deriving (Show)
{-
--type DenseBetaU = BetaU LA.Vector
--type SparseBetaU = BetaU SLA.SpVector

toDenseBetaU :: SparseBetaU -> DenseBetaU
toDenseBetaU (BetaU svBeta svU) =
  BetaU (SD.toDenseVector svBeta) (SD.toDenseVector svU)

toSparseBetaU :: DenseBetaU -> SparseBetaU
toSparseBetaU (BetaU vBeta vU) =
  BetaU (SD.toSparseVector vBeta) (SD.toSparseVector vU)
-}
scaleBetaU :: Double -> BetaU -> BetaU
scaleBetaU x (BetaU b u) = BetaU (LA.scale x b) (SLA.scale x u)

addBetaU :: BetaU -> BetaU -> BetaU
addBetaU (BetaU bA uA) (BetaU bB uB) = BetaU (LA.add bA bB) (uA SLA.^+^ uB)


data RandomEffectCalculated = RandomEffectCalculated RandomEffectModelMatrix (CovarianceVec -> LambdaMatrix)

type PMatrix = SLA.SpMatrix Double -- should this be (SLA.SpMatrix Int) ??

data DevianceType = ML | REML deriving (Show, Eq)

effectsForGroup :: GroupFitSpec -> Int
effectsForGroup (GroupFitSpec _ b vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForGroup #-}

colsForGroup :: GroupFitSpec -> Int
colsForGroup l = nCategories l * effectsForGroup l
{-# INLINABLE colsForGroup #-}

data ParameterEstimates =
  ParameterEstimates
  { pEstimate :: LA.Vector Double
  , pCovariance :: LA.Matrix Double
  }

type FixedParameterEstimates = ParameterEstimates
type GroupParameterEstimates = VB.Vector FixedParameterEstimates

data MinimizeDevianceVerbosity = MDVNone | MDVSimple

minimizeDeviance
  :: EffectsIO r
  => MinimizeDevianceVerbosity
  -> DevianceType
  -> GeneralizedLinearMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec -- ^ initial guess for theta
  -> P.Sem r
       ( CovarianceVec
       , Double
       , Double
       , BetaU
       , LA.Vector Double
       , CholeskySolutions
       ) -- ^ (theta, profiled_deviance, sigma2, beta, u, b, cholesky blocks)
minimizeDeviance verbosity dt glmm reCalc@(RandomEffectCalculated smZ mkLambda) th0 =
  P.wrapPrefix "minimizeDeviance" $ do
    P.logLE P.Info "Starting..."
    when (generalized glmm && (dt == REML)) $ P.throw $ OtherGLMError "Can't use REML for generalized models."
    P.logLE P.Diagnostic "Beginning Cholmod analysis of Z matrix structure."
    cholmodAnalysis <- cholmodAnalyzeProblem reCalc
    P.logLE P.Diagnostic "Finished Cholmod analysis of Z matrix structure."
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      pd :: EffectsIO r => CovarianceVec -> P.Sem r ( Double
                                                    , Double
                                                    , BetaU
                                                    , SLA.SpVector Double
                                                    , CholeskySolutions
                                                    )
      pd x = profiledDeviance pdv cholmodAnalysis dt glmm reCalc x
      obj x = unsafePerformEffects $ fmap (\(d, _, _, _, _) -> d) $ pd x
      stop                = NL.ObjectiveAbsoluteTolerance 1e-6 NL.:| []
      MixedModel _ levels = mixedModel glmm
      thetaLB             = thetaLowerBounds levels
      algorithm           = NL.BOBYQA obj [thetaLB] Nothing
      problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
      expThetaLength      = FL.fold
        (FL.premap
         (\l -> let e = effectsForGroup l in e + (e * (e - 1) `div` 2))
         FL.sum
        )
        levels
    when (LA.size th0 /= expThetaLength)
      $  P.throw
      $  OtherGLMError
      $  "guess for theta has "
      <> (T.pack $ show $ LA.size th0)
      <> " entries but should have "
      <> (T.pack $ show expThetaLength)
      <> "."
    let eSol = NL.minimizeLocal problem th0
    case eSol of
      Left  result                       -> P.throw (OtherGLMError $ T.pack $ show result)
      Right (NL.Solution pdS thS result) -> do
        liftIO ( putStrLn
                 $  "Solution ("
                 ++ show result
                 ++ ") reached! At th="
                 ++ show thS
               )
        (pdVal, sigma2, betaU, svb, cs) <- pd thS
        return (thS, pdVal, sigma2, betaU, SD.toDenseVector svb, cs)

-- ugh.  But I dont know a way in NLOPT to have bounds on some not others.
thetaLowerBounds :: FitSpecByGroup g -> NL.Bounds
thetaLowerBounds groupFSM =
  let negInfinity :: Double = negate $ 1 / 0
  in  NL.LowerBounds $ setCovarianceVector groupFSM 0 negInfinity --FL.fold fld levels

setCovarianceVector :: FitSpecByGroup g -> Double -> Double -> CovarianceVec
setCovarianceVector groupFSM diag offDiag = FL.fold fld groupFSM
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
  :: forall g r
   . (Enum g, Bounded g, A.Ix g, Effects r)
  => FixedPredictors
  -> FitSpecByGroup g
  -> GLM.RowClassifier g
  -> P.Sem r RandomEffectModelMatrix
makeZ mX groupFSM rc = do
  let
    maybeError :: T.Text -> Maybe a -> P.Sem r a
    maybeError err = maybe (P.throw $ OtherGLMError $ err) return 
    (nO, nP) = LA.size mX
    k        = FL.fold FL.length groupFSM -- number of levels
--    groupSize g = fmap nCategories $ M.lookup g groupFSM 
    q        = FL.fold FL.sum $ fmap colsForGroup groupFSM -- total number of columns in Z
    predictor rowIndex fixedEffectIndex =
      mX `LA.atIndex` (rowIndex, fixedEffectIndex)
    -- construct Z for level as a fold over rows
    entries :: Int -> g -> GroupFitSpec -> Int -> P.Sem r [(Int, Int, Double)]
    entries colOffset group groupFS rowIndex = do
      categoryN <- maybeError "Error in categorNumberFromRowIndex"
                   $ GLM.categoryNumberFromRowIndex rc rowIndex group
      let
        entryOffset = colOffset + (effectsForGroup groupFS * categoryN)
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
      return $ intercept ++ slopes
    ziFoldM
      :: g -> Int -> GroupFitSpec -> FL.FoldM (P.Sem r) Int [(Int, Int, Double)]
    ziFoldM group groupOffset groupFS =
      let q = colsForGroup groupFS
      in  FL.FoldM
            (\mes n -> fmap (mes ++) $ entries groupOffset group groupFS n)
            (return [])
            return
    groups       = IS.members $ GLM.groupIndices rc -- [minBound .. maxBound] --Numbers = VB.generate (VB.length groupFSs) id
    groupOffsets = FL.prescan FL.sum $ fmap colsForGroup $ M.elems groupFSM
  zisM <- FL.foldM
    ( sequenceA
    $ fmap (\(group, gOffset, gSpec) -> ziFoldM group gOffset gSpec)
    $ zip3 groups groupOffsets (M.elems groupFSM)
    )
    [0 .. (nO - 1)]
  return $ SLA.fromListSM (nO, q) $ concat zisM

-- TODO: Add levels checks
-- TODO: Add Lambda checks
checkProblem :: Effects r => MixedModel b g -> RandomEffectCalculated -> P.Sem r ()
checkProblem (MixedModel (RegressionModel _ mX vY) _) (RandomEffectCalculated smZ _)
  = do
    let (n, _)     = LA.size mX
        yRows      = LA.size vY
        (zRows, _) = SLA.dim smZ
    when (zRows /= n)
      $  P.throw
      $  OtherGLMError
      $  (T.pack $ show n)
      <> " cols in X but "
      <> (T.pack $ show zRows)
      <> " cols in Z"
    when (yRows /= n)
      $  P.throw
      $  OtherGLMError
      $  (T.pack $ show n)
      <> " cols in X but "
      <> (T.pack $ show yRows)
      <> " entries in Y"
    return ()
-- For each level, i, with p_i random effects,
-- Lambda is built from lower-triangular blocks
-- which are p_i x p_i.  Theta contains the values
-- for each block, in col major order
makeLambda :: FitSpecByGroup g -> (CovarianceVec -> LambdaMatrix)
makeLambda groupFSM =
  let blockF vTH = FL.Fold
        (\(bs, vx) l' -> ((l', vx) : bs, VS.drop (effectsForGroup l') vx))
        ([], vTH)
        (reverse . fst)
      blockData vTh = FL.fold (blockF vTh) groupFSM
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

type CholmodFactor = ( CH.ForeignPtr CH.Common -- ^ pre-allocated CHOMOD common space
                     , CH.ForeignPtr CH.Factor -- ^ precomputed pattern work on Z
                     , PMatrix -- ^ permutation matrix from above
                     )

cholmodAnalyzeProblem
  :: EffectsIO r => RandomEffectCalculated -> P.Sem r CholmodFactor
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

-- this one is IO since we'll need to unsafePerformIO it
profiledDeviance
  :: EffectsIO r
  => ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> GeneralizedLinearMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec -- ^ theta
  -> P.Sem r 
       ( Double
       , Double
       , BetaU
       , SLA.SpVector Double
       , CholeskySolutions
       ) -- ^ (pd, sigma^2, beta, u, b) 
profiledDeviance verbosity cf dt glmm reCalc@(RandomEffectCalculated smZ mkLambda) vTh
  = P.wrapPrefix "profiledDeviance" $ do
    liftIO $ hSetBuffering stdout NoBuffering
    let MixedModel (RegressionModel _ mX vY) _ = mixedModel glmm
        n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        lambda = mkLambda vTh
        smZS   = smZ SLA.## lambda --(smT SLA.## smS)
--        smZSt  = SLA.transpose smZS
    P.logLE P.Diagnostic "Starting. Calling getCholeskySolutions."
    (cs@(CholeskySolutions smLth smRzx mRxx), betaU) <- getCholeskySolutions
      cf
      glmm
      GLM.UseCanonical -- FIX
      reCalc
      vTh
    when (verbosity == PDVAll) $ liftIO $ do
      putStrLn "Z"
      LA.disp 2 $ SD.toDenseMatrix smZ
      putStrLn "Lambda"
      LA.disp 2 $ SD.toDenseMatrix lambda
      putStrLn "Lth"
      LA.disp 2 $ SD.toDenseMatrix smLth
      putStrLn "Rzx"
      LA.disp 2 $ SD.toDenseMatrix smRzx
      putStrLn "Rxx"
      LA.disp 2 mRxx
    let
      svb           = lambda SLA.#> (svU betaU)
--      vBetaU        = toDenseBetaU svBetaU
      vEta          = denseLinearPredictor mX smZS betaU
      (pd, rTheta2) = profiledDeviance' verbosity
                                        dt
                                        glmm
                                        reCalc
                                        vTh
                                        cs
                                        vEta
                                        (svU betaU)
      dof = case dt of
        ML   -> realToFrac n
        REML -> realToFrac (n - p)
      sigma2 = rTheta2 / dof
    when (verbosity == PDVAll) $ liftIO $ do
      let (_, _, smP) = cf
      putStrLn $ "smP="
      LA.disp 1 $ SD.toDenseMatrix smP
      putStrLn $ "rTheta^2=" ++ show rTheta2
      let logLth = logDetTriangularSM smLth
          logDet = case dt of
            ML   -> logLth
            REML -> logLth + logDetTriangularM mRxx
      putStrLn $ "2 * logLth=" ++ show (2 * logLth)
      putStrLn $ "2 * logDet=" ++ show (2 * logDet)
    when (verbosity == PDVAll || verbosity == PDVSimple) $ liftIO $ do
      putStrLn $ "pd(th=" ++ (show vTh) ++ ") = " ++ show pd
    return (pd, sigma2, betaU, svb, cs)

profiledDeviance'
  :: ProfiledDevianceVerbosity
  -> DevianceType
  -> GeneralizedLinearMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> CholeskySolutions
  -> EtaVec
  -> UVec
  -> (Double, Double) -- profiled deviance, deviance + u'u 
profiledDeviance' pdv dt glmm re vTh chol vEta svU =
  let MixedModel (RegressionModel _ mX vY) _ = mixedModel glmm
      (          RandomEffectCalculated smZ                    mkLambda) = re
      (observationDistribution, vW) = case glmm of
        LMM _       -> (GLM.Normal, (LA.fromList $ L.replicate (LA.size vY) 1.0)) 
        GLMM _ vW' od -> (od, vW')
      logLth        = logDetTriangularSM (lTheta chol) 
      dev2 = GLM.deviance observationDistribution GLM.UseCanonical vW vY vEta
      rTheta2       = dev2 + (SLA.norm2Sq svU)
      n             = LA.size vY
      (_  , p     ) = LA.size mX
      (dof, logDet) = case dt of
        ML   -> (realToFrac n, logLth)
        REML -> (realToFrac (n - p), logLth + (logDetTriangularM $ rXX chol))
      pd = (2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof)))
  in  (pd, rTheta2)


type LinearPredictor = LA.Vector Double
type SparseLinearPredictor = SLA.SpVector Double

denseLinearPredictor
  :: FixedPredictors -> SLA.SpMatrix Double -> BetaU -> LinearPredictor
denseLinearPredictor mX smZS betaU =
  mX LA.#> (vBeta betaU)
  + (SD.toDenseVector $ smZS SLA.#> svU betaU)

sparseLinearPredictor
  :: FixedPredictors -> SLA.SpMatrix Double -> BetaU -> SparseLinearPredictor
sparseLinearPredictor mX smZS betaU =
  SD.toSparseVector (mX LA.#> vBeta betaU)
  SLA.^+^ (smZS SLA.#> svU betaU)

upperTriangular :: Int -> Int -> a -> Bool
upperTriangular r c _ = (r <= c)

lowerTriangular :: Int -> Int -> a -> Bool
lowerTriangular r c _ = (r >= c)

type LTheta = SLA.SpMatrix Double
type Rzx = SLA.SpMatrix Double
type Rxx = LA.Matrix Double
type Beta = LA.Vector Double

data CholeskySolutions = CholeskySolutions { lTheta :: LTheta
                                           , rZX :: Rzx
                                           , rXX :: Rxx
                                           } deriving (Show)

data NormalEquations = NormalEquationsLMM | NormalEquationsGLMM WMatrix MuVec (SLA.SpVector Double)

normalEquationsRHS
  :: EffectsIO r
  => NormalEquations
  -> CholmodFactor
  -> UMatrix
  -> VMatrix
  -> Observations
  -> P.Sem r (SLA.SpMatrix Double, LA.Vector Double)
normalEquationsRHS NormalEquationsLMM (cholmodC, cholmodF, _) smU mV vY
  = P.wrapPrefix "normalEquationsRHS" $ do
  let svRhsZ = SLA.transpose smU SLA.#> SD.toSparseVector vY
      vRhsX  = LA.tr mV LA.#> vY
  P.logLE P.Diagnostic "Calling CHOLMOD.solveSparse to multiply by P"
  smRhsZ <- liftIO $ CH.solveSparse cholmodC
            cholmodF
            CH.CHOLMOD_P
            (SD.svColumnToSM svRhsZ)
  return (smRhsZ, vRhsX) --(SLA.transpose smU SLA.#> SD.toSparseVector vY, LA.tr mV LA.#> vY)

normalEquationsRHS (NormalEquationsGLMM vW vMu svU) (cholmodC, cholmodF, _) smU mV vY
  = P.wrapPrefix "normalEquationsRHS" $ do
    let vX   = VS.zipWith (*) (VS.map sqrt vW) (vY - vMu)
        svUX = SLA.transpose smU SLA.#> (SD.toSparseVector vX)
    P.logLE P.Diagnostic "Calling CHOLMOD.solveSparse to multiply by P"    
    smPUX <- liftIO $ CH.solveSparse cholmodC
                            cholmodF
                            CH.CHOLMOD_P
                            (SD.svColumnToSM svUX)
    P.logLE P.Diagnostic $ "Finished. dim(smPUX)=" <> (T.pack $ show $ SLA.dim smPUX)
    let smRhsZ = smPUX SLA.^-^ (SD.svColumnToSM svU)
        vRhsX    = LA.tr mV LA.#> vX
    P.logLE P.Diagnostic $ "dim(smRhsZ)=" <> (T.pack $ show $ SLA.dim smRhsZ)        
    return (smRhsZ, vRhsX)

cholmodCholeskySolutionsLMM
  :: EffectsIO r
  => CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> P.Sem r (CholeskySolutions, BetaU)
cholmodCholeskySolutionsLMM cholmodFactor mixedModel randomEffCalcs vTh = do
  let (MixedModel             (RegressionModel _ mX vY) levels  ) = mixedModel
      (RandomEffectCalculated smZ mkLambda) = randomEffCalcs
      lambda = mkLambda vTh
      smZS   = smZ SLA.## lambda
  cholmodCholeskySolutions' cholmodFactor smZS mX NormalEquationsLMM mixedModel

observationDistribution :: GeneralizedLinearMixedModel b g -> GLM.ObservationDistribution
observationDistribution (LMM _) = GLM.Normal
observationDistribution (GLMM _ _ od) = od

linkFunctionType
  :: GeneralizedLinearMixedModel b g -> GLM.UseLink -> GLM.LinkFunctionType
linkFunctionType glmm ul = case ul of
  GLM.UseOther x   -> x
  GLM.UseCanonical -> GLM.canonicalLink $ observationDistribution glmm
  
spUV
  :: GeneralizedLinearMixedModel b g
  -> GLM.UseLink
  -> RandomEffectCalculated
  -> CovarianceVec
  -> LinearPredictor
  -> (UMatrix, VMatrix)

spUV (LMM (MixedModel (RegressionModel _ mX _) _)) _ (RandomEffectCalculated smZ mkLambda) vTh _
  = (smZ SLA.#~# mkLambda vTh, mX)

spUV glmm@(GLMM (MixedModel (RegressionModel _ mX _) _) vW _) ul (RandomEffectCalculated smZ mkLambda) vTh vEta
  = let smZS = smZ SLA.#~# mkLambda vTh
        mWdMudEta = LA.diag $ weights
          mX
          smZS
          vW
          (GLM.linkFunction $ linkFunctionType glmm ul)
          vEta
        smU = SD.toSparseMatrix mWdMudEta SLA.#~# smZS
        mV  = mWdMudEta LA.<> mX
    in  (smU, mV)

type ZStarMatrix = SLA.SpMatrix Double

weights
  :: FixedPredictors
  -> ZStarMatrix
  -> WMatrix
  -> GLM.LinkFunction
  -> LinearPredictor
  -> WMatrix
weights mX smZS vW lf vEta =
  let vdMudEta = VS.map (GLM.derivInv lf) vEta
  in  VS.zipWith (*) (VS.map sqrt vW) vdMudEta

data ShrinkPD = Shrunk EtaVec UVec | NotShrunk Double EtaVec UVec Double deriving (Show)

getCholeskySolutions
  :: EffectsIO r
  => CholmodFactor
  -> GeneralizedLinearMixedModel b g
  -> GLM.UseLink
  -> RandomEffectCalculated
  -> CovarianceVec
  -> P.Sem r (CholeskySolutions, BetaU)
getCholeskySolutions cf (LMM mm) _ reCalc vTh =
  P.wrapPrefix "getCholeskySolutions" $ do
  P.logLE P.Diagnostic "Calling cholmodCholeskySolutionsLMM"
  cholmodCholeskySolutionsLMM cf mm reCalc vTh

-- TODO: There is much low hanging fruit here in terms of not re-calcing things
--
getCholeskySolutions cf glmm@(GLMM mm vW od) ul reCalc vTh = P.wrapPrefix "getCholeskySolutions" $ do
  P.logLE P.Diagnostic "Starting (GLMM)..."
  let MixedModel (RegressionModel _ mX vY) levels = mm
      (RandomEffectCalculated smZ mkLambda) = reCalc
      smZS        = smZ SLA.#~# mkLambda vTh
      n           = LA.size vY
      (_, q)      = SLA.dim smZS
      (_, _, smP) = cf
      vEta0      =  computeEta0 (linkFunctionType glmm ul) vY
      svU0         = SD.toSparseVector $ LA.fromList $ L.replicate q 0 -- FIXME (maybe 0 is right here??)
      lf = GLM.linkFunction $ linkFunctionType glmm ul
      deltaCCS vEta svU = P.wrapPrefix "deltaCCS" $ do
        P.logLE P.Diagnostic "Starting..."
        let (smU, mV) = spUV glmm ul reCalc vTh vEta
            vMu       = VS.map (GLM.invLink lf) vEta
            neqs      = NormalEquationsGLMM vW vMu svU
        (chol, dBetaU) <- cholmodCholeskySolutions' cf smU mV neqs mm
        P.logLE P.Diagnostic "Finished."
        return (dBetaU, chol, smU)
      incS betaU dBetaU x =
        addBetaU betaU (scaleBetaU x dBetaU)
      refineDBetaU vEta svU' dBetaU chol = P.wrapPrefix "refineDBetaU" $ do
        P.logLE P.Diagnostic $ "Starting \nvEta="
          <> (T.pack $ show vEta)
          <> "\nsvU'=" <> (T.pack $ show svU')
          <> "\ndBetaU=" <> (T.pack $ show dBetaU)
        let pdF x y = fst $ profiledDeviance' PDVNone ML glmm reCalc vTh chol x y 
            pd0 = pdF vEta svU'            
            checkOne x =              
              let svU'' = svU' SLA.^+^ SLA.scale x (svU dBetaU) --svBetaU' = incS svBetaU svdBetaU x
                  vdEta = denseLinearPredictor mX smZS dBetaU
                  vEta' = vEta + LA.scale x vdEta
                  pdNew = pdF vEta' svU''
              in  if pdNew < pd0 then Shrunk vEta' svU'' else NotShrunk pd0 vEta' svU'' pdNew
            check triesLeft x = do
              P.logLE P.Diagnostic $ "check (triesLeft=" <> (T.pack $ show triesLeft) <> "; x=" <> (T.pack $ show x) <> ")..."
              case checkOne x of
                ns@(NotShrunk _ _ _ _) ->
                  if triesLeft > 1
                  then (P.logLE P.Diagnostic (T.pack $ show ns) >> check (triesLeft - 1) (x / 2))
                  else P.throw $ OtherGLMError "Too many step-halvings in getCholeskySolutions."
                Shrunk x y -> return (x,y)
        check 10 1.0
      iterateOne vEta svU' = P.wrapPrefix "iterateOne" $ do
        P.logLE P.Diagnostic "Starting..."
        (dBetaU, chol, smU) <- deltaCCS vEta svU'
        refined <- refineDBetaU vEta svU' dBetaU chol
        P.logLE P.Diagnostic "Finished."
        return (refined, chol, dBetaU, smU)
      iterateUntilConverged n convergedOCC vEta svU' = P.wrapPrefix "iterateUntilConverged" $ do
        P.logLE P.Diagnostic $ "Starting (n=" <> (T.pack $ show n) <> ")..."
        ((vEta', svU''), chol, dBetaU, smU) <- iterateOne vEta svU'
        let occ = orthogonalConvergence
                  glmm
                  ul
                  smZS
                  smP
                  (lTheta chol)
                  smU
                  vY
                  vEta'
                  svU''
                  (svU dBetaU)  
        if occ < convergedOCC
          then (P.logLE P.Diagnostic "Finished." >> return (vEta', svU', chol))
          else
          (if n == 1
            then P.throw (OtherGLMError "Too many iterations in getCholeskySolutions.")
            else iterateUntilConverged (n - 1) convergedOCC vEta' svU' 
          )
      finish (vEta, svU', chol) =
        (chol, BetaU (betaFromXZStarEtaU mX smZS svU' (SD.toSparseVector vEta)) svU')
  P.logLE P.Diagnostic $ "vY=" <> (T.pack $ show vY)
  P.logLE P.Diagnostic $ "vEta0=" <> (T.pack $ show vEta0)
  res <- fmap finish $ iterateUntilConverged 10 0.001 vEta0 svU0
  P.logLE P.Diagnostic "Finished."
  return res

orthogonalConvergence
  :: GeneralizedLinearMixedModel b g
  -> GLM.UseLink
  -> ZStarMatrix -- Z*Lambda(theta)
  -> PMatrix
  -> LMatrix
  -> UMatrix
  -> Observations
  -> EtaVec
  -> UVec
  -> UVec
  -> Double
orthogonalConvergence glmm ul smZS smP smL smU vY vEta svU svdU =
  let MixedModel (RegressionModel _ mX vY) _ = mixedModel glmm
      lft = case ul of
        GLM.UseOther x   -> x
        GLM.UseCanonical -> GLM.canonicalLink $ case glmm of
          LMM  _      -> GLM.Normal
          GLMM _ _ od -> od
      vMu        = VS.map (GLM.invLink $ GLM.linkFunction lft) vEta
      svYMinusMu = SD.toSparseVector $ vY - vMu
      x1         = SLA.transpose smL SLA.#> svdU 
      x2 =
        smP
          SLA.#>  (SLA.transpose smU SLA.#> svYMinusMu)
          SLA.^-^ svU
      q' = realToFrac $ snd $ SLA.dim smZS
      n' = realToFrac $ LA.size vY
  in  (SLA.norm2 x1 / sqrt q') / (SLA.norm2 x2 / sqrt (n' - q'))

-- compute a starting point for the linear predictor.  Just assume the given observations but fix any out of bounds
-- That is, assume mu = y* and then eta = g (mu) = g (y*) where
-- y* is y, but with some care at the boundaries.  That is, y but with all elements in the domain of g.
computeEta0 :: GLM.LinkFunctionType -> Observations -> LA.Vector Double
computeEta0 lft vY =
  let lF = GLM.link $ GLM.linkFunction lft
      lowerBound lb x = if x < lb then lb else x
      upperBound ub x = if x > ub then ub else x
      bounded lb ub = lowerBound lb . upperBound ub
      f x = case lft of
        GLM.IdentityLink -> x
        GLM.LogisticLink -> bounded 0.00001 0.99999 x --minimum (maximum (0.0001, x), 0.99999)
        GLM.ExponentialLink -> lowerBound 0.0001 x
  in VS.map (lF . f) vY

-- NB: eta = X <*> Beta + Z* <*> u
-- X <*> Beta = eta - Z* <*> u
-- X'X <*> Beta = X' <*> (eta - Z* <*> u)
-- Beta = inverse(X'X) <*> X' <*> (eta - Z* <*> u)
-- and X'X is symmetric positive definite so we solve the above via Cholesky
betaFromXZStarEtaU
  :: FixedPredictors
  -> SLA.SpMatrix Double
  -> UVec
  -> SLA.SpVector Double
  -> BetaVec
betaFromXZStarEtaU mX smZS svU svEta =
  let vRhs = LA.tr mX LA.#> (SD.toDenseVector $ (smZS SLA.#> svU) SLA.^-^ svEta)
      cXtX = LA.chol $ LA.trustSym $ LA.tr mX LA.<> mX
  in head $ LA.toColumns $ LA.cholSolve cXtX (LA.asColumn vRhs)

{-
glmmCholeskyStep
  :: CholmodFactor
  -> NormalEquations
  -> MixedModel
  -> WMatrix
  -> LinkFunctionType
  -> SLA.SpVector Double
  -> SLA.SpVector Double
  -> IO CholeskySolutions
glmmCholeskyStep mm vW lft svBeta svU =
-}


cholmodCholeskySolutions'
  :: EffectsIO r
  => CholmodFactor
  -> UMatrix -- U (ZS in linear case)
  -> VMatrix -- V (X in linear case)
  -> NormalEquations
  -> MixedModel b g
  -> P.Sem r (CholeskySolutions, BetaU)
cholmodCholeskySolutions' cholmodFactor smU mV nes mixedModel = P.wrapPrefix "cholmodCholeskySolutions'" $ do
  let (cholmodC, cholmodF, smP) = cholmodFactor
      (MixedModel (RegressionModel _ _ vY) _) = mixedModel
      n                         = LA.size vY
      (_, p)                    = LA.size mV
      (_, q)                    = SLA.dim smU
  -- Cholesky factorize to get L_theta *and* update factor for solving with it
  P.logLE P.Diagnostic "Starting..."
  liftIO $ CH.spMatrixFactorize cholmodC cholmodF CH.SquareSymmetricLower
    $ SLA.filterSM lowerTriangular
    $ xTxPlusI
    $ smU
  -- compute Rzx
  P.logLE P.Diagnostic "Calling normal equations."
  (smPRhsZ, vRhsX) <- normalEquationsRHS nes cholmodFactor smU mV vY
  P.logLE P.Diagnostic "Calling solveSparse for svC."
  svC <- SLA.toSV <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPRhsZ)  
  let smUtV = (SLA.transpose smU) SLA.#~# (SD.toSparseMatrix mV)
  P.logLE P.Diagnostic "Calling solveSparse for smPUtV."
  smPUtV <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smUtV -- PU'V
  P.logLE P.Diagnostic "Calling solveSparse for smRzx."
  smRzx  <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPUtV -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
  -- compute Rxx
  -- NB: c is defined as Lu + Rzx(Beta) which allows solving the normal equations in pieces
  let vTvMinusRzxTRzx =
        (LA.tr mV) LA.<> mV - (SD.toDenseMatrix $ smRzx SLA.#^# smRzx)
      mRxx            = LA.chol $ LA.trustSym $ vTvMinusRzxTRzx
  -- now we have the entire Cholesky decomposition.  Solve the normal equations
      vRzxtC          = SD.toDenseVector $ (SLA.transposeSM smRzx) SLA.#> svC
      betaSols        = LA.cholSolve mRxx (LA.asColumn $ vRhsX - vRzxtC)
      vBeta           = head $ LA.toColumns $ betaSols
      svBeta          = SD.toSparseVector vBeta
      svRzxBeta       = smRzx SLA.#> svBeta
      smCMinusRzxBeta = SD.svColumnToSM (svC SLA.^-^ svRzxBeta)
  P.logLE P.Diagnostic "Calling solveSparse for smPu."
  smPu  <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Lt smCMinusRzxBeta
  P.logLE P.Diagnostic "Calling solveSparse for svu."
  svu   <- SLA.toSV <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Pt $ smPu)
  -- NB: This has to happen after the solves because it unfactors the factor...
  -- NB: It does *not* undo the analysis
  P.logLE P.Diagnostic "Calling choleskyFactorSM for smLth."
  smLth <- liftIO $ CH.choleskyFactorSM cholmodF cholmodC
  P.logLE P.Diagnostic "Finished."
  return $ (CholeskySolutions smLth smRzx mRxx, BetaU vBeta svu)



---- Unused
{-
cholmodCholeskySolutions_Old
  :: CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> IO (CholeskySolutions, BetaU)
cholmodCholeskySolutions_Old cholmodFactor mixedModel randomEffCalcs vTh = do
  let (cholmodC, cholmodF, _) = cholmodFactor
      (MixedModel             (RegressionModel _ mX vY) levels  ) = mixedModel
      (RandomEffectCalculated smZ mkLambda) = randomEffCalcs
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
-- NB: cu is defined as Lu + Rzx(Beta) which allows solving the equations in pieces
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
  return (CholeskySolutions smLth smRzx mRx, BetaU svBeta svu)
-}
