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
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Data.IndexedSet               as IS

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Knit.Effect.Logger            as P

--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import qualified Control.Exception             as X
import           Control.Monad                  ( when
                                                , join
                                                )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Data.Array                    as A
import           Data.Function                  ( (&) )
import qualified Data.List                     as L
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA
import           Data.Typeable                  ( Typeable )

import           Numeric.LinearAlgebra.Sparse   ( (#>) )

import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL
import           System.IO.Unsafe               ( unsafePerformIO )



import qualified Data.Map                      as M
import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

import           System.IO                      ( hSetBuffering
                                                , stdout
                                                , BufferMode(..)
                                                )

data GLMError = OtherGLMError T.Text | NonGLMError T.Text deriving (Show, Typeable)
instance X.Exception GLMError

type Effects r = (P.Member (P.Error GLMError) r, P.LogWithPrefixesLE r)
type EffectsIO r = ( Effects r
                   , P.Member (P.Embed IO) r)


type GLMEffects a = P.Sem '[P.Logger P.LogEntry, P.PrefixLog, P.Error GLMError, P.Embed IO, P.Final IO] a

runEffectsVerboseIO :: GLMEffects a -> IO (Either GLMError a)
runEffectsVerboseIO action =
  action
    & P.filteredLogEntriesToIO P.logAll
    & P.errorToIOFinal
    & P.embedToFinal
    & P.runFinal


runEffectsIO :: GLMEffects a -> IO (Either GLMError a)
runEffectsIO action =
  action
    & P.filteredLogEntriesToIO P.nonDiagnostic
    & P.errorToIOFinal
    & P.embedToFinal
    & P.runFinal


unsafePerformEffects :: GLMEffects a -> a -> a
unsafePerformEffects action def =
  action
    & runEffectsIO
    & fmap (either X.throwIO return)
    & join
    & (\x -> X.catch
        x
        (\(e :: X.SomeException) -> putStrLn ("Error: " ++ show e) >> return def
        )
      )
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
  when (not $ IS.subset indexedGroupEffects indexedFixedEffects)
    $ P.throw
    $ OtherGLMError
        "group contains effects not in fixed effects in \"makeGroupFitSpec\""
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
    maybe (P.throw $ OtherGLMError lookupError) return
    $ traverse lookupGroupSize
    $ M.toList ebg
  fmap M.fromList $ traverse makeSpec groupInfoList

data RegressionModelSpec b = RegressionModelSpec { rmsFixedEffects :: (GLM.FixedEffects b)
                                             , rmsFixedPredictors :: FixedPredictors
                                             , rmsObservations :: Observations
                                             } deriving (Show, Eq)

data MixedModelSpec b g = MixedModelSpec { mmsRegressionModelSpec :: (RegressionModelSpec b),
                                           mmsFitSpecByGroup :: (FitSpecByGroup g)
                                         } deriving (Show, Eq)

data LMMOptimizer = LMM_BOBYQA | LMM_NELDERMEAD | LMM_SBPLX | LMM_NEWUOA_BOUND deriving (Show, Eq)

lmmOptimizerToNLOPT
  :: LMMOptimizer
  -> (  NL.Objective
     -> [NL.Bounds]
     -> Maybe NL.InitialStep
     -> NL.LocalAlgorithm
     )
lmmOptimizerToNLOPT LMM_BOBYQA       = NL.BOBYQA
lmmOptimizerToNLOPT LMM_NELDERMEAD   = NL.NELDERMEAD
lmmOptimizerToNLOPT LMM_SBPLX        = NL.SBPLX
lmmOptimizerToNLOPT LMM_NEWUOA_BOUND = NL.NEWUOA_BOUND

data LMMControls = LMMControls { lmmOptimizer :: LMMOptimizer
                               , lmmOptimizerTolerance :: Double
                               } deriving (Show, Eq)

defaultLMMControls :: LMMControls
defaultLMMControls = LMMControls LMM_BOBYQA 1e-6

data PIRLSConvergenceType = PCT_Eta | PCT_Deviance | PCT_Orthogonal deriving (Show, Eq)

data PIRLSConvergenceCriterion = PIRLSConvergenceCriterion { pirlsConvergenceType :: PIRLSConvergenceType
                                                           , pirlsTolerance :: Double
                                                           , pirlsMaxSteps :: Int } deriving (Show, Eq)

data GLMMControls = GLMMControls { glmLink :: GLM.UseLink
                                 , pirlsMaxStepHalvings :: Int
                                 , pirlsConvergenceCriterion :: PIRLSConvergenceCriterion
                                 } deriving (Show, Eq)

defaultGLMMControls :: GLMMControls
defaultGLMMControls = GLMMControls
  GLM.UseCanonical
  10
  (PIRLSConvergenceCriterion PCT_Deviance 1e-7 20)

data LinearMixedModelSpec b g =
  LinearMixedModelSpec { lmmsMixedModelSpec :: MixedModelSpec b g
                       , lmmsControls :: LMMControls
                       } deriving (Show)

data GeneralizedLinearMixedModelSpec b g =
  GeneralizedLinearMixedModelSpec { glmmsLinearMixedModelSpec :: LinearMixedModelSpec b g
                                  , glmmsWeights :: WMatrix -- this should prolly be in LinearMixedModelSpec
                                  , glmmsObservationsDistribution :: GLM.ObservationsDistribution
                                  , glmmsControls :: GLMMControls
                                  } deriving (Show)

data MixedModel b g = LinearMixedModel (LinearMixedModelSpec b g)
                    | GeneralizedLinearMixedModel (GeneralizedLinearMixedModelSpec b g)

mixedModelSpec :: MixedModel b g -> MixedModelSpec b g
mixedModelSpec (LinearMixedModel lmmSpec) = lmmsMixedModelSpec lmmSpec
mixedModelSpec (GeneralizedLinearMixedModel glmmSpec) =
  lmmsMixedModelSpec . glmmsLinearMixedModelSpec $ glmmSpec

regressionModelSpec :: MixedModel b g -> RegressionModelSpec b
regressionModelSpec = mmsRegressionModelSpec . mixedModelSpec

lmmControls :: MixedModel b g -> LMMControls
lmmControls (LinearMixedModel lmmSpec) = lmmsControls lmmSpec
lmmControls (GeneralizedLinearMixedModel glmmSpec) =
  lmmsControls . glmmsLinearMixedModelSpec $ glmmSpec

generalized :: MixedModel b g -> Bool
generalized (LinearMixedModel            _) = False
generalized (GeneralizedLinearMixedModel _) = True

fixedPredictors :: MixedModel b g -> FixedPredictors
fixedPredictors = rmsFixedPredictors . regressionModelSpec

observations :: MixedModel b g -> Observations
observations = rmsObservations . regressionModelSpec


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

data BetaU = BetaU { bu_vBeta :: LA.Vector Double, bu_svU :: SLA.SpVector Double } deriving (Show)

scaleBetaU :: Double -> BetaU -> BetaU
scaleBetaU x (BetaU b u) = BetaU (LA.scale x b) (SLA.scale x u)

addBetaU :: BetaU -> BetaU -> BetaU
addBetaU (BetaU bA uA) (BetaU bB uB) = BetaU (LA.add bA bB) (uA SLA.^+^ uB)

--diffBetaU :: BetaU -> BetaU -> BetaU
--diffBetaU (BetaU bA uA) (BetaU bB uB) = BetaU (bA - bB) (uA SLA.^-^ uB)

data RandomEffectCalculated = RandomEffectCalculated { recModelMatrix:: RandomEffectModelMatrix
                                                     , recMkLambda :: (CovarianceVec -> LambdaMatrix)
                                                     }

makeZS :: RandomEffectCalculated -> CovarianceVec -> SLA.SpMatrix Double
makeZS reCalc vTh =
  let smZ      = recModelMatrix reCalc
      mkLambda = recMkLambda reCalc
  in  smZ SLA.## (mkLambda vTh)

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
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec -- ^ initial guess for theta
  -> P.Sem
       r
       ( CovarianceVec
       , Double
       , Double
       , BetaU
       , LA.Vector Double
       , CholeskySolutions
       ) -- ^ (theta, profiled_deviance, sigma2, beta, u, b, cholesky blocks)
minimizeDeviance verbosity dt mm reCalc th0 =
  P.wrapPrefix "minimizeDeviance" $ do
    P.logLE P.Diagnostic "Starting..."
    when (generalized mm && (dt == REML)) $ P.throw $ OtherGLMError
      "Can't use REML for generalized models."
    P.logLE P.Diagnostic "Beginning Cholmod analysis of Z matrix structure."
    cholmodAnalysis <- cholmodAnalyzeProblem reCalc
    P.logLE P.Diagnostic "Finished Cholmod analysis of Z matrix structure."
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      (objectiveOs, finalOs) = if generalized mm
        then (Optim_GLMM_PIRLS, Optim_GLMM_Final)
        else (Optim_LMM, Optim_LMM)
      pd
        :: EffectsIO r
        => OptimizationStage
        -> CovarianceVec
        -> P.Sem
             r
             ( Double
             , Double
             , BetaU
             , SLA.SpVector Double
             , CholeskySolutions
             )
      pd os x = profiledDeviance pdv cholmodAnalysis dt mm reCalc os x
      obj x = unsafePerformEffects
        (fmap (\(d, _, _, _, _) -> d) $ pd objectiveOs x)
        0
      levels    = mmsFitSpecByGroup $ mixedModelSpec mm
      thetaLB   = thetaLowerBounds levels
      algorithm = (lmmOptimizerToNLOPT $ lmmOptimizer $ lmmControls mm)
        obj
        [thetaLB]
        Nothing
      stop =
        NL.ObjectiveAbsoluteTolerance (lmmOptimizerTolerance $ lmmControls mm)
          NL.:| []
      problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
      expThetaLength = FL.fold
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
      Left result -> P.throw (OtherGLMError $ T.pack $ show result)
      Right (NL.Solution pdS thS result) -> do
        P.logLE P.Info
          $  "Solution ("
          <> (T.pack $ show result)
          <> ") reached! At th="
          <> (T.pack $ show thS)
        (pdVal, sigma2, betaU, svb, cs) <- pd finalOs thS
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
      categoryN <-
        maybeError "Error in categorNumberFromRowIndex"
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
-- TODO: parse, don't validate!  Which I think means typed sizes here.
checkProblem
  :: Effects r => MixedModel b g -> RandomEffectCalculated -> P.Sem r ()
checkProblem mm reCalc = do
  let mX         = rmsFixedPredictors $ regressionModelSpec mm
      vY         = rmsObservations $ regressionModelSpec mm
      smZ        = recModelMatrix reCalc
      (n, _)     = LA.size mX
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
  return () --(MixedModel (RegressionModel _ mX vY) _) (RandomEffectCalculated smZ _)
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

xTx :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTx smX = (smX SLA.#^# smX)

xTxPlusI :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTxPlusI smX = (smX SLA.#^# smX) SLA.^+^ (SLA.eye $ SLA.ncols smX)

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
cholmodAnalyzeProblem reCalc = do
  cholmodC <- liftIO CH.allocCommon
  liftIO $ CH.startC cholmodC
  liftIO $ CH.setFinalLL 1 cholmodC -- I don't quite get this.  We should be able to solve LDx = b either way??
  (cholmodF, smP) <-
    liftIO
    $ CH.spMatrixAnalyzeWP cholmodC CH.SquareSymmetricLower
    $ xTx
    $ recModelMatrix reCalc
--  liftIO $ CH.printFactor cholmodF "After analyze, before factorize" cholmodC
--  liftIO $ putStrLn "smP=" >> (LA.disp 1 $ SD.toDenseMatrix smP)
  return (cholmodC, cholmodF, smP)

data ProfiledDevianceVerbosity = PDVNone | PDVSimple | PDVAll deriving (Show, Eq)
data OptimizationStage = Optim_LMM | Optim_GLMM_PIRLS | Optim_GLMM_Final deriving (Show, Eq)

-- this one is IO since we'll need to unsafePerformIO it
profiledDeviance
  :: EffectsIO r
  => ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> MixedModel b g
  -> RandomEffectCalculated
  -> OptimizationStage
  -> CovarianceVec -- ^ theta
  -> P.Sem
       r
       ( Double
       , Double
       , BetaU
       , SLA.SpVector Double
       , CholeskySolutions
       ) -- ^ (pd, sigma^2, beta, u, b) 
profiledDeviance verbosity cf dt mm reCalc os vTh =
  P.wrapPrefix "profiledDeviance" $ do
    liftIO $ hSetBuffering stdout NoBuffering
    let mX     = rmsFixedPredictors $ regressionModelSpec mm
        smZ    = recModelMatrix reCalc
        lambda = recMkLambda reCalc vTh
        vY     = rmsObservations $ regressionModelSpec mm
        n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        smZS   = makeZS reCalc vTh

    P.logLE P.Diagnostic "Starting. Calling getCholeskySolutions."
    (cs@(CholeskySolutions smLth smRzx mRxx), betaU) <- getCholeskySolutions
      cf
      mm
      reCalc
      os
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
    let svb     = lambda SLA.#> (bu_svU betaU)
        vEta    = denseLinearPredictor mX smZS betaU
        osFinal = if (generalized mm) then Optim_GLMM_Final else Optim_LMM
    (pd, rTheta2) <- profiledDeviance' verbosity
                                       dt
                                       mm
                                       reCalc
                                       vTh
                                       osFinal
                                       cs
                                       (bu_vBeta betaU)
                                       (bu_svU betaU)
    let dof = case dt of
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
  :: EffectsIO r
  => ProfiledDevianceVerbosity
  -> DevianceType
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> OptimizationStage
  -> CholeskySolutions
  -> BetaVec
  -> UVec
  -> P.Sem r (Double, Double) -- profiled deviance, deviance + u'u 
profiledDeviance' pdv dt mm reCalc vTh os chol vBeta svU =
  P.wrapPrefix "profiledDeviance'" $ do
    P.logLE P.Diagnostic
      $  "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vTh)
      <> ")..."
--    P.logLE P.Diagnostic $ "vEta=" <> (T.pack $ show vEta)
--    P.logLE P.Diagnostic $ "svU=" <> (T.pack $ show svU)
    let
      mX       = rmsFixedPredictors $ regressionModelSpec mm
      vY       = rmsObservations $ regressionModelSpec mm
      (od, vW) = case mm of
        LinearMixedModel _ ->
          (GLM.Normal, (LA.fromList $ L.replicate (LA.size vY) 1.0))
        GeneralizedLinearMixedModel glmmSpec ->
          (glmmsObservationsDistribution glmmSpec, glmmsWeights glmmSpec)
      logLth = logDetTriangularSM (lTheta chol)
      smZS   = makeZS reCalc vTh
      vEta   = denseLinearPredictor mX smZS (BetaU vBeta svU)
      vMu = LA.cmap (GLM.invLink $ GLM.linkFunction (linkFunctionType mm)) vEta
      uDotu  = svU SLA.<.> svU
    P.logLE P.Diagnostic
      $  "logLth="
      <> (T.pack $ show logLth)
      <> "; u.u="
      <> (T.pack $ show uDotu)
    (pd, rTheta2) <- case mm of
      LinearMixedModel _ -> do
        when (os /= Optim_LMM)
          $ P.throw
          $ OtherGLMError
              "In profiledDeviance': OptimzationStage must be Optim_LMM for all LMM calls"
        let
          n             = LA.size vY
          vW            = VS.replicate n 1.0
          devResidual   = GLM.deviance (GLM.Normal) vW vY vMu
          rTheta2       = devResidual + uDotu
          (_  , p     ) = LA.size mX
          (dof, logDet) = case dt of
            ML -> (realToFrac n, logLth)
            REML ->
              (realToFrac (n - p), logLth + (logDetTriangularM $ rXX chol))
        P.logLE P.Diagnostic $ "LMM devResid=: " <> (T.pack $ show devResidual)
        return
          $ ((2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof))), rTheta2)
      GeneralizedLinearMixedModel _ -> do
        let devResidual = GLM.deviance od vW vY vMu --GLM.devianceCondAbs od vW vY vMu
            aic         = GLM.aicR od vW vY vMu devResidual
            rTheta2     = devResidual + uDotu
        P.logLE P.Diagnostic
          $  "GLMM: devResid= "
          <> (T.pack $ show devResidual)
          <> "; aicR="
          <> (T.pack $ show aic)
        case os of
          Optim_LMM ->
            P.throw
              $ OtherGLMError
              $ "In profiledDeviance': OptimizationStage was Optim_LMM in a GLMM call"
          Optim_GLMM_PIRLS -> return $ (devResidual + uDotu, rTheta2)
          Optim_GLMM_Final -> return $ (aic + uDotu + 2 * logLth, rTheta2)
    P.logLE P.Diagnostic $ "; pd=" <> (T.pack $ show pd)
    return (pd, rTheta2)


type LinearPredictor = LA.Vector Double
type SparseLinearPredictor = SLA.SpVector Double

denseLinearPredictor
  :: FixedPredictors -> SLA.SpMatrix Double -> BetaU -> LinearPredictor
denseLinearPredictor mX smZS betaU =
  mX LA.#> (bu_vBeta betaU) + (SD.toDenseVector $ smZS SLA.#> bu_svU betaU)

sparseLinearPredictor
  :: FixedPredictors -> SLA.SpMatrix Double -> BetaU -> SparseLinearPredictor
sparseLinearPredictor mX smZS betaU =
  SD.toSparseVector (mX LA.#> bu_vBeta betaU) SLA.^+^ (smZS SLA.#> bu_svU betaU)

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

data NormalEquations = NormalEquationsLMM | NormalEquationsGLMM WMatrix GLM.LinkFunction EtaVec UVec

normalEquationsRHS
  :: EffectsIO r
  => NormalEquations
  -> UMatrix
  -> VMatrix
  -> Observations
  -> P.Sem r (SLA.SpVector Double, LA.Vector Double)
normalEquationsRHS NormalEquationsLMM smU mV vY =
  P.wrapPrefix "normalEquationsRHS" $ do
    let svRhsZ = SLA.transpose smU SLA.#> SD.toSparseVector vY
        vRhsX  = LA.tr mV LA.#> vY
--    P.logLE P.Diagnostic $ "RHS=" <> (T.pack $ show (svRhsZ, vRhsX))
    return (svRhsZ, vRhsX) --(SLA.transpose smU SLA.#> SD.toSparseVector vY, LA.tr mV LA.#> vY)

normalEquationsRHS (NormalEquationsGLMM vVarW lf vEta svU) smU mV vY =
  P.wrapPrefix "normalEquationsRHS" $ do
{-  
    P.logLE P.Diagnostic "Starting..."

    P.logLE P.Diagnostic $ "vVarW=" <> (T.pack $ show vVarW)

    P.logLE P.Diagnostic $ "svU=" <> (T.pack $ show svU)
    P.logLE P.Diagnostic $ "smU=" <> (T.pack $ show smU)
    P.logLE P.Diagnostic $ "mV=" <> (T.pack $ show mV)
    P.logLE P.Diagnostic $ "vMu=" <> (T.pack $ show vMu)
    P.logLE P.Diagnostic $ "vY=" <> (T.pack $ show vY)
-}
    let vMu         = VS.map (GLM.invLink lf) vEta
--        vMuEta      = VS.map (GLM.derivInv lf) vEta
--        vYMinusMu   = vY - vMu
        vWtYMinusMu = VS.zipWith3 (\wt y mu -> sqrt wt * (y - mu)) vVarW vY vMu
  {-
      vWrkResids  = VS.zipWith3 (\y mu muEta -> (y - mu) / muEta) vY vMu vMuEta
      vWrkResp    = VS.zipWith (+) vEta vWrkResids
      vWtWrkResp  = VS.zipWith3
        (\wt muEta wrkResp -> muEta * (sqrt wt) * wrkResp)
        vVarW
        vMuEta
        vWrkResp
-}
--        vZ          = VS.zipWith (*) vMuEta vWtYMinusMu
        vX          = vWtYMinusMu
{-      
    P.logLE P.Diagnostic
      $  "vWtYminusMu="
      <> (T.pack $ show vWtYMinusMu)
      <> "\nvWtWrkResp="
      <> (T.pack $ show vWtWrkResp)
-}
    let svUX   = SLA.transpose smU SLA.#> (SD.toSparseVector vX)
        svRhsZ = svUX SLA.^-^ svU
        vRhsX  = LA.tr mV LA.#> vX
--    P.logLE P.Diagnostic $ "RHS=" <> (T.pack $ show (svRhsZ, vRhsX))
    return (svRhsZ, vRhsX)

cholmodCholeskySolutionsLMM
  :: EffectsIO r
  => CholmodFactor
  -> MixedModelSpec b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> P.Sem r (CholeskySolutions, BetaU)
cholmodCholeskySolutionsLMM cholmodFactor mixedModelSpec reCalc vTh = do
  let mX     = rmsFixedPredictors $ mmsRegressionModelSpec mixedModelSpec
      vY     = rmsObservations $ mmsRegressionModelSpec mixedModelSpec
      levels = mmsFitSpecByGroup mixedModelSpec
      smZS   = makeZS reCalc vTh
  cholmodCholeskySolutions' cholmodFactor
                            smZS
                            mX
                            NormalEquationsLMM
                            mixedModelSpec

observationDistribution :: MixedModel b g -> GLM.ObservationsDistribution
observationDistribution (LinearMixedModel _) = GLM.Normal
observationDistribution (GeneralizedLinearMixedModel glmmSpec) =
  glmmsObservationsDistribution glmmSpec

useLink :: MixedModel b g -> GLM.UseLink
useLink (LinearMixedModel _) = GLM.UseCanonical
useLink (GeneralizedLinearMixedModel glmmSpec) =
  glmLink $ glmmsControls glmmSpec

linkFunctionType :: MixedModel b g -> GLM.LinkFunctionType
linkFunctionType mm = case useLink mm of
  GLM.UseOther x   -> x
  GLM.UseCanonical -> GLM.canonicalLink $ observationDistribution mm

getCholeskySolutions
  :: EffectsIO r
  => CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> OptimizationStage
  -> CovarianceVec
  -> P.Sem r (CholeskySolutions, BetaU)
getCholeskySolutions cf mm@(LinearMixedModel _) reCalc _ vTh =
  P.wrapPrefix "getCholeskySolutions (LMM)" $ do
    P.logLE P.Diagnostic "Starting..."
    P.logLE P.Diagnostic $ "theta=" <> (T.pack $ show vTh)
    P.logLE P.Diagnostic "Calling cholmodCholeskySolutionsLMM"
    res <- cholmodCholeskySolutionsLMM cf (mixedModelSpec mm) reCalc vTh
    P.logLE P.Diagnostic "Finished"
    return res
-- TODO: There is much low hanging fruit here in terms of not re-calcing things
--
getCholeskySolutions cf mm@(GeneralizedLinearMixedModel glmmSpec) reCalc os vTh
  = P.wrapPrefix "getCholeskySolutions (GLMM)" $ do
    P.logLE P.Diagnostic
      $  "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vTh)
      <> ")..."
    let mX          = rmsFixedPredictors $ regressionModelSpec mm
        vY          = rmsObservations $ regressionModelSpec mm
        vW          = glmmsWeights glmmSpec
        od          = glmmsObservationsDistribution glmmSpec
        GLMMControls ul maxHalvings pirlsCC = glmmsControls glmmSpec
        smZS        = makeZS reCalc vTh
        n           = LA.size vY
        (_, q)      = SLA.dim smZS
        (_, _, smP) = cf
        vEtaExact   = computeEta0 (linkFunctionType mm) vY
        svU0        = SD.toSparseVector $ LA.fromList $ L.replicate q 0 -- FIXME (maybe 0 is right here??)
        vBeta0      = betaFromXZStarEtaU mX smZS svU0 vEtaExact -- approximate fit with no random effects
--        vEta0       = denseLinearPredictor mX smZS (BetaU vBeta0 svU0)
{-        
    P.logLE P.Diagnostic
      $  "vEtaExact ="
      <> (T.pack $ show vEtaExact)
      <> "\nsvU0="
      <> (T.pack $ show svU0)
      <> "\nvBeta0="
      <> (T.pack $ show vBeta0)
      <> "\nvEta0="
      <> (T.pack $ show vEta0)
-}
    let
      lf = GLM.linkFunction $ linkFunctionType mm
      pdFunction os chol x y =
        fst <$> profiledDeviance' PDVNone ML mm reCalc vTh os chol x y
      --incS betaU dBetaU x = addBetaU betaU (scaleBetaU x dBetaU)
      iterateUntilConverged n pdCurrent vBeta svU = do
        P.logLE P.Diagnostic
          $  "iterateUntilConverged (n="
          <> (T.pack $ show n)
          <> "; pd="
          <> (T.pack $ show pdCurrent)
          <> ")"

        (pd', vBeta', svU', chol, dBetaU, smU) <- updateEtaBetaU
          cf
          mm
          maxHalvings
          reCalc
          (pdFunction os)
          vBeta
          svU
          vTh
        cc <- case pirlsConvergenceType pirlsCC of
          PCT_Eta -> do
            let vEta  = denseLinearPredictor mX smZS (BetaU vBeta svU)
                vEta' = denseLinearPredictor mX smZS (BetaU vBeta' svU')
                dEta  = vEta' - vEta
            return $ (LA.norm_2 dEta) / (LA.norm_2 vEta')
          PCT_Deviance -> return $ abs $ (pdCurrent - pd') / pd'
          PCT_Orthogonal ->
            P.throw $ OtherGLMError "ConvergeOrthognal not implemented yet."
{-          
        occ <- pirlsConvergence glmm
                                smZS
                                smP
                                (lTheta chol)
                                smU
                                vEta
                                vEta'
                                svU'
                                (bu_svU dBetaU)
-}
        case (cc < (pirlsTolerance pirlsCC), n) of
          (True, _) ->
            P.logLE P.Diagnostic "Finished." >> return (vBeta', svU', chol)
          (False, 1) -> P.throw
            (OtherGLMError "Too many iterations in getCholeskySolutions.")
          (False, m) -> do
            P.logLE P.Diagnostic $ "Not converged.  cc=" <> (T.pack $ show cc)
            iterateUntilConverged (m - 1) pd' vBeta' svU'
    (vBeta', svU', chol) <- iterateUntilConverged (pirlsMaxSteps pirlsCC)
                                                  (1.0 / 0.0) -- Infinity :: Double
                                                  vBeta0
                                                  svU0
--    let vBeta' = betaFromXZStarEtaU mX smZS svU' vEta'
    P.logLE P.Diagnostic "Finished."
    return (chol, BetaU vBeta' svU')

updateEtaBetaU
  :: EffectsIO r
  => CholmodFactor
  -> MixedModel b g
  -> Int
  -> RandomEffectCalculated
  -> (CholeskySolutions -> BetaVec -> UVec -> P.Sem r Double)
  -> BetaVec
  -> UVec
  -> CovarianceVec
  -> P.Sem
       r
       ( Double
       , BetaVec
       , UVec
       , CholeskySolutions
       , BetaU
       , SLA.SpMatrix Double
       )
updateEtaBetaU cf mm maxHalvings reCalc pdFunction vBeta svU vTh =
  P.wrapPrefix "updateEtaBetaU" $ do
    P.logLE P.Diagnostic $ "Starting..." ---(Eta=" <> (T.pack $ show vEta) <> ")"
    (dBetaU, chol, smU) <- compute_dBetaU cf mm reCalc vBeta svU vTh
--    P.logLE P.Diagnostic $ "d" <> (T.pack $ show dBetaU)
    (pdFinal, vBeta', svU', dBetaU') <- refine_dBetaU mm
                                                      maxHalvings
                                                      reCalc
                                                      (pdFunction chol)
                                                      vBeta
                                                      svU
                                                      dBetaU
                                                      vTh
    P.logLE P.Diagnostic
      $  "Finished (||dU||^2 = "
      <> (T.pack $ show $ (bu_svU dBetaU' SLA.<.> bu_svU dBetaU'))
      <> "; ||dBeta||^2="
      <> (T.pack $ show $ (bu_vBeta dBetaU' LA.<.> bu_vBeta dBetaU'))
      <> ")." --- (Eta=" <> (T.pack $ show vEta') <> ")"
    return (pdFinal, vBeta', svU', chol, dBetaU', smU)

data ShrinkPD = Shrunk Double BetaVec UVec BetaU | NotShrunk Double BetaVec UVec Double deriving (Show)

refine_dBetaU
  :: EffectsIO r
  => MixedModel b g
  -> Int -- max step halvings
  -> RandomEffectCalculated
  -> (BetaVec -> UVec -> P.Sem r Double) -- profiled deviance
  -> BetaVec
  -> UVec
  -> BetaU
  -> CovarianceVec
  -> P.Sem r (Double, BetaVec, UVec, BetaU)
refine_dBetaU mm maxHalvings reCalc pdF vBeta svU dBetaU vTh =
  P.wrapPrefix "refineDBetaU" $ do
    P.logLE P.Diagnostic $ "Starting..."
    {-        
          <>  "\nvEta="
          <> (T.pack $ show vEta)
          <> "\nsvU'="
          <> (T.pack $ show svU')
          <> "\ndBetaU="
          <> (T.pack $ show dBetaU)
-}
    let --pdF x y = fst <$> profiledDeviance' PDVNone ML glmm reCalc vTh chol x y
        mX   = rmsFixedPredictors $ regressionModelSpec mm
        vY   = rmsObservations $ regressionModelSpec mm
        smZS = makeZS reCalc vTh
    pd0 <- pdF vBeta svU
    let
      checkOne x = do
        let deltaBetaU = scaleBetaU x dBetaU
            svU'       = svU SLA.^+^ (bu_svU deltaBetaU) --svBetaU' = incS svBetaU svdBetaU x
            vBeta'     = vBeta + (bu_vBeta deltaBetaU)
--            vEta = denseLinearPredictor mX smZS (BetaU vBeta svU)
--            vdEta = denseLinearPredictor mX smZS dBetaU
--            vEta' = vEta + LA.scale x vdEta
        pdNew <- pdF vBeta' svU'
        return $ if pdNew <= pd0
          then Shrunk pdNew vBeta' svU' deltaBetaU
          else NotShrunk pd0 vBeta' svU' pdNew
      check triesLeft x = do
        P.logLE P.Diagnostic
          $  "check (triesLeft="
          <> (T.pack $ show triesLeft)
          <> "; x="
          <> (T.pack $ show x)
          <> ")..."
        co <- checkOne x
        case co of
          ns@(NotShrunk _ _ _ _) -> if triesLeft > 1
            then
              (do
--P.logLE P.Diagnostic (T.pack $ show ns)
                check (triesLeft - 1) (x / 2)
              )
            else
              P.throw
              $  OtherGLMError
              $  "Too many step-halvings in getCholeskySolutions. (dBetaU = "
              <> (T.pack $ show dBetaU)
          sh@(Shrunk pdFinal vBetaFinal svUFinal deltaBetaU) -> do
--            P.logLE P.Diagnostic $ "refined: " <> (T.pack $ show sh)
            return (pdFinal, vBetaFinal, svUFinal, deltaBetaU)
    check maxHalvings 1.0


-- glmer: updateXwts but also with lambda
spUV
  :: MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> EtaVec
  -> (UMatrix, VMatrix)
spUV mm@(LinearMixedModel _) reCalc vTh _ =
  let mX       = rmsFixedPredictors $ regressionModelSpec mm
      smZ      = recModelMatrix reCalc
      mkLambda = recMkLambda reCalc
  in  (smZ SLA.#~# mkLambda vTh, mX)

spUV mm@(GeneralizedLinearMixedModel glmmSpec) reCalc vTh vEta =
  let mX    = rmsFixedPredictors $ regressionModelSpec mm
      vW    = glmmsWeights glmmSpec
      smZS  = makeZS reCalc vTh
      vWUV  = weightsForUV mm vW vEta
      smWUV = SLA.mkDiagonal (VS.length vWUV) $ VS.toList vWUV
      smU   = smWUV SLA.## smZS
      mV    = (LA.diag vWUV) LA.<> mX
  in  (smU, mV)

type ZStarMatrix = SLA.SpMatrix Double

-- glmer: sqrtWrkWt
weightsForUV :: MixedModel b g -> WMatrix -> LinearPredictor -> WMatrix
weightsForUV mm vW vEta =
  let lf        = GLM.linkFunction $ linkFunctionType mm
      vMu       = VS.map (GLM.invLink lf) vEta
      vdMudEta  = VS.map (GLM.derivInv lf) vEta
      vVariance = GLM.scaledVariance (observationDistribution mm) vMu
      vVSW      = GLM.varianceScaledWeights (observationDistribution mm) vW vMu --    
  in  VS.zipWith (\muEta vsw -> muEta * sqrt vsw) vdMudEta vVSW
--    VS.zipWith3 (\w muEta var -> muEta * sqrt (w / var)) vW vdMudEta vVariance

compute_dBetaU
  :: EffectsIO r
  => CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> BetaVec
  -> UVec
  -> CovarianceVec
  -> P.Sem
       r
       (BetaU, CholeskySolutions, SLA.SpMatrix Double)
compute_dBetaU cf mm@(GeneralizedLinearMixedModel glmmSpec) reCalc vBeta svU vTh
  = P.wrapPrefix "compute_dBetaU" $ do
    P.logLE P.Diagnostic "Starting..."
    let mX        = rmsFixedPredictors $ regressionModelSpec mm
        smZS      = makeZS reCalc vTh
        vW        = glmmsWeights glmmSpec
        lf        = GLM.linkFunction $ linkFunctionType mm
        vEta      = denseLinearPredictor mX smZS (BetaU vBeta svU)
        (smU, mV) = spUV mm reCalc vTh vEta
        vMu       = VS.map (GLM.invLink lf) vEta
        vVarW = GLM.varianceScaledWeights (observationDistribution mm) vW vMu
        neqs      = NormalEquationsGLMM vVarW lf vEta svU
    (chol, dBetaU) <- cholmodCholeskySolutions' cf
                                                smU
                                                mV
                                                neqs
                                                (mixedModelSpec mm)
    P.logLE P.Diagnostic "Finished."
    return (dBetaU, chol, smU)

compute_dBetaU _ (LinearMixedModel _) _ _ _ _ = P.throw
  $ OtherGLMError "Called compute_dBetaU given an LMM instead of a GLMM."
{-
pirlsConvergence
  :: EffectsIO r
  => GeneralizedLinearMixedModel b g
  -> ZStarMatrix -- Z*Lambda(theta)
  -> PMatrix
  -> LMatrix
  -> UMatrix
  -> EtaVec
  -> EtaVec
  -> UVec
  -> UVec
  -> P.Sem r Double
pirlsConvergence (LMM _ _) _ _ _ _ _ _ _ _ =
  P.throw $ OtherGLMError "LMM given to pirlsConvergence"

pirlsConvergence glmm@(GLMM mm _ _ glmmC) smZS smP smL smU vEta vEta' svU svdU
  = P.wrapPrefix "pirlsConvergence" $ do
    let (GLMMControls _ _ _ pirlsC) = glmmC
--    P.logLE P.Diagnostic $ "U=" <> (T.pack $ show svU)
--    P.logLE P.Diagnostic $ "dU=" <> (T.pack $ show svdU)
    case pirlsC of
      ConvergeEta _ _ -> do
        let dEta = vEta' - vEta
--        P.logLE P.Diagnostic $ "Eta=" <> (T.pack $ show vEta')
--        P.logLE P.Diagnostic $ "dEta=" <> (T.pack $ show dEta)
        return $ (LA.norm_2 dEta) / (LA.norm_2 vEta')
      ConvergeOrthogonal _ _ -> do
        let MixedModel (RegressionModel _ mX vY) _ = mm
            lft        = linkFunctionType glmm
            vMu        = VS.map (GLM.invLink $ GLM.linkFunction lft) vEta
            svYMinusMu = SD.toSparseVector $ (vY - vMu)
            x1         = SLA.transpose smL SLA.#> svdU
            x2         = smL SLA.#> x1
            x3 = smP SLA.#> (SLA.transpose smU SLA.#> svYMinusMu) SLA.^-^ svU
            q'         = realToFrac $ snd $ SLA.dim smZS
            n'         = realToFrac $ LA.size vY
--        P.logLE P.Diagnostic $ "L'dU=" <> (T.pack $ show x1)
--        P.logLE P.Diagnostic $ "LL'dU=" <> (T.pack $ show x2)
--        P.logLE P.Diagnostic $ "(y - mu)" <> (T.pack $ show svYMinusMu)
--        P.logLE P.Diagnostic $ "P(smU')*(y - mu) - u=" <> (T.pack $ show x3)
        return $ (SLA.norm2 x1 / sqrt q') / (SLA.norm2 x2 / sqrt (n' - q'))
-}
-- compute a starting point for the linear predictor.  Just assume the given observations but fix any out of bounds
-- That is, assume mu = y* and then eta = g (mu) = g (y*) where
-- y* is y, but with some care at the boundaries.  That is, y but with all elements in the domain of g.
computeEta0 :: GLM.LinkFunctionType -> Observations -> EtaVec
computeEta0 lft vY =
  let lF = GLM.link $ GLM.linkFunction lft
      lowerBound lb x = if x < lb then lb else x
      upperBound ub x = if x > ub then ub else x
      bounded lb ub = lowerBound lb . upperBound ub
      f x = case lft of
        GLM.IdentityLink    -> x
        GLM.LogisticLink    -> bounded 0.00001 0.99999 x --minimum (maximum (0.0001, x), 0.99999)
        GLM.ExponentialLink -> lowerBound 0.0001 x
  in  VS.map (lF . f) vY

-- NB: eta = X <*> Beta + Z* <*> u
-- X <*> Beta = eta - Z* <*> u
-- X'X <*> Beta = X' <*> (eta - Z* <*> u)
-- Beta = inverse(X'X) <*> X' <*> (eta - Z* <*> u)
-- and X'X is symmetric positive definite so we solve the above via Cholesky
betaFromXZStarEtaU
  :: FixedPredictors -> SLA.SpMatrix Double -> UVec -> EtaVec -> BetaVec
betaFromXZStarEtaU mX smZS svU vEta =
  let vZSu = SD.toDenseVector $ smZS SLA.#> svU
      vRhs = LA.tr mX LA.#> (vEta - vZSu)
      cXtX = LA.chol $ LA.trustSym $ LA.tr mX LA.<> mX
  in  head $ LA.toColumns $ LA.cholSolve cXtX (LA.asColumn vRhs)

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
  -> MixedModelSpec b g
  -> P.Sem r (CholeskySolutions, BetaU)
cholmodCholeskySolutions' cholmodFactor smU mV nes mixedModelSpec =
  P.wrapPrefix "cholmodCholeskySolutions'" $ do
    let (cholmodC, cholmodF, smP) = cholmodFactor
        vY = rmsObservations $ mmsRegressionModelSpec mixedModelSpec
        n                         = LA.size vY
        (_, p)                    = LA.size mV
        (_, q)                    = SLA.dim smU
    -- Cholesky factorize to get L_theta *and* update factor for solving with it
    P.logLE P.Diagnostic "Starting..."
    -- TODO: Cholmod has built in support for factorizing XtX + a * I.  Use it.
    let cfs = CH.FactorizeAtAPlusBetaI 1
    liftIO $ CH.spMatrixFactorizeP cholmodC
                                   cholmodF
                                   cfs
                                   CH.UnSymmetric
                                   (SLA.transpose smU)
{-
    let cfs = CH.FactorizeA
    liftIO
      $ CH.spMatrixFactorizeP cholmodC cholmodF cfs CH.SquareSymmetricLower
      $ SLA.filterSM lowerTriangular
      $ xTxPlusI
      $ smU
-}
--    P.logLE P.Diagnostic "After factorize..."
    -- compute Rzx
--    P.logLE P.Diagnostic "Calling normal equations."
    (svRhsZ, vRhsX) <- normalEquationsRHS nes smU mV vY
  {-  P.logLE P.Diagnostic
      $  "Normal Equation RHS:\nsmPRhsZ="
      <> (T.pack $ show smPRhsZ)
      <> "\nvRhsX="
      <> (T.pack $ show vRhsX)
    P.logLE P.Diagnostic "Calling solveSparse for svC."
-}
    smPRhsZ         <- liftIO
      $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P (SD.svColumnToSM svRhsZ)
--    P.logLE P.Diagnostic "After smPRhsZ..."
    svC <-
      SLA.toSV
        <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPRhsZ)
    let smUtV = (SLA.transpose smU) SLA.#~# (SD.toSparseMatrix mV)
--    P.logLE P.Diagnostic "Calling solveSparse for smPUtV."
    smPUtV <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smUtV -- PU'V
--    P.logLE P.Diagnostic "Calling solveSparse for smRzx."
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
--    P.logLE P.Diagnostic "Calling solveSparse for smPu."
    smPu <- liftIO
      $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Lt smCMinusRzxBeta
--    P.logLE P.Diagnostic "Calling solveSparse for svu."
    svu <-
      SLA.toSV
        <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Pt $ smPu)
    -- NB: This has to happen after the solves because it unfactors the factor...
    -- NB: It does *not* undo the analysis
--    P.logLE P.Diagnostic "Calling choleskyFactorSM for smLth."
    smLth <- liftIO $ CH.choleskyFactorSM cholmodF cholmodC
    P.logLE P.Diagnostic
      $  "min(u)="
      <> (T.pack $ show $ FL.fold FL.minimum svu)
      <> "; max(u)="
      <> (T.pack $ show $ FL.fold FL.maximum svu)
      <> "; min(beta)="
      <> (T.pack $ show $ VS.minimum vBeta)
      <> "; max(beta)="
      <> (T.pack $ show $ VS.maximum vBeta)
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
