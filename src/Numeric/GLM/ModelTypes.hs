{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}

module Numeric.GLM.ModelTypes
  ( GLMError(..)
  , Effects
  , EffectsIO
  , GLMEffects
  , FixedPredictors
  , Observations
  , GroupFitSpec(..)
  , makeGroupFitSpec
  , FitSpecByGroup
  , fitSpecByGroup
  , RegressionModelSpec(..)
  , MixedModelSpec(..)
  , LMMOptimizer(..)
  , LMMControls(..)
  , defaultLMMControls
  , PIRLSConvergenceType(..)
  , PIRLSConvergenceCriterion(..)
  , GLMMControls(..)
  , defaultGLMMControls
  , LinearMixedModelSpec(..)
  , GeneralizedLinearMixedModelSpec(..)
  , MixedModel(..)
  , mixedModelSpec
  , regressionModelSpec
  , lmmControls
  , generalized
  , fixedPredictors
  , weights
  , observations
  , changeMMObservations
  , RandomEffectModelMatrix
  , CovarianceVec
  , LambdaMatrix
  , WMatrix
  , UMatrix
  , VMatrix
  , LMatrix
  , RzxMatrix
  , RxxMatrix
  , BetaVec
  , UVec
  , MuVec
  , EtaVec
--  , LinearPredictorComponents
  , LinearPredictor(..)
  , denseLinearPredictor
  , diffLinearPredictor
--  , sparseLinearPredictor
--  , SparseEtaVec
  , BetaU(..)
  , scaleBetaU
  , addBetaU
  , makeZS
  , ZStarMatrix(smZS)
  , RandomEffectCalculated(..)
  , PMatrix
  )
where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Numeric.SparseDenseConversions
                                               as SD

import           Control.Monad                  ( when )
import qualified Data.List                     as L
import qualified Data.Map                      as M
import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS
import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA

import qualified Control.Exception             as X
import           Data.Typeable                  ( Typeable )
import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Polysemy.Reader                as P
import qualified Knit.Effect.Logger            as P


data GLMError = OtherGLMError T.Text | NonGLMError T.Text deriving (Show, Typeable)
instance X.Exception GLMError

type Effects r = (P.Member (P.Error GLMError) r, P.LogWithPrefixesLE r)
type EffectsIO r = ( Effects r
                   , P.Member (P.Embed IO) r)


type GLMEffects a = P.Sem '[P.Reader P.LogWithPrefixIO, P.Logger P.LogEntry, P.PrefixLog, P.Error GLMError, P.Embed IO, P.Final IO] a

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
  :: (GLM.PredictorC b, Effects r)
  => Int -- ^ number of items in group
  -> GLM.FixedEffects b
  -> GLM.IndexedEffectSet b
  -> P.Sem r GroupFitSpec
makeGroupFitSpec n fixedEffects indexedGroupEffects = do
  let indexedFixedEffects = GLM.indexedFixedEffectSet fixedEffects
  -- even if we don't have an intercept in the fixed effects, we can have
  -- group intercepts. I think.
  when (not $ IS.subset indexedGroupEffects (IS.add indexedFixedEffects GLM.Intercept))
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
  :: (GLM.PredictorC b, GLM.GroupC g, Effects r)
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

data RegressionModelSpec b = RegressionModelSpec { rmsFixedEffects :: GLM.FixedEffects b
                                                 , rmsFixedPredictors :: FixedPredictors
                                                 , rmsObservations :: Observations
                                                 } deriving (Show, Eq)

changeRMSObservations
  :: Observations -> RegressionModelSpec b -> RegressionModelSpec b
changeRMSObservations os (RegressionModelSpec fes fps _) =
  RegressionModelSpec fes fps os

data MixedModelSpec b g = MixedModelSpec { mmsRegressionModelSpec :: (RegressionModelSpec b),
                                           mmsFitSpecByGroup :: (FitSpecByGroup g)
                                         } deriving (Show, Eq)

changeMMSObservations
  :: Observations -> MixedModelSpec b g -> MixedModelSpec b g
changeMMSObservations os (MixedModelSpec rms fsbg) =
  MixedModelSpec (changeRMSObservations os rms) fsbg

data LMMOptimizer = LMM_BOBYQA | LMM_NELDERMEAD | LMM_SBPLX | LMM_NEWUOA_BOUND deriving (Show, Eq)


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

changeLMMSObservations
  :: Observations -> LinearMixedModelSpec b g -> LinearMixedModelSpec b g
changeLMMSObservations os (LinearMixedModelSpec mms c) =
  LinearMixedModelSpec (changeMMSObservations os mms) c

data GeneralizedLinearMixedModelSpec b g =
  GeneralizedLinearMixedModelSpec { glmmsLinearMixedModelSpec :: LinearMixedModelSpec b g
                                  , glmmsWeights :: WMatrix -- this should prolly be in LinearMixedModelSpec
                                  , glmmsObservationsDistribution :: GLM.ObservationsDistribution
                                  , glmmsControls :: GLMMControls
                                  } deriving (Show)


changeGLMMSObservations
  :: Observations
  -> GeneralizedLinearMixedModelSpec b g
  -> GeneralizedLinearMixedModelSpec b g
changeGLMMSObservations os (GeneralizedLinearMixedModelSpec lmms vW od c) =
  GeneralizedLinearMixedModelSpec (changeLMMSObservations os lmms) vW od c


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

changeMMObservations :: Observations -> MixedModel b g -> MixedModel b g
changeMMObservations os (LinearMixedModel lmms) =
  LinearMixedModel (changeLMMSObservations os lmms)
changeMMObservations os (GeneralizedLinearMixedModel glmms) =
  GeneralizedLinearMixedModel (changeGLMMSObservations os glmms)

weights :: MixedModel b g -> WMatrix
weights mm@(LinearMixedModel _) = VS.replicate (VS.length $ observations mm) 1
weights (GeneralizedLinearMixedModel glmms) = glmmsWeights glmms

type RandomEffectModelMatrix = SLA.SpMatrix Double
type CovarianceVec = LA.Vector Double
type LambdaMatrix = SLA.SpMatrix Double

type WMatrix = LA.Vector Double -- constant weights
type UMatrix = SLA.SpMatrix Double
type VMatrix = LA.Matrix Double
type LMatrix = SLA.SpMatrix Double -- Cholesky Z block
type RzxMatrix = SLA.SpMatrix Double
type RxxMatrix = LA.Matrix Double
type PMatrix = SLA.SpMatrix Double -- should this be (SLA.SpMatrix Int) ??

type BetaVec = LA.Vector Double
type UVec = SLA.SpVector Double
type MuVec = LA.Vector Double
type EtaVec = LA.Vector Double

data BetaU = BetaU { bu_vBeta :: LA.Vector Double, bu_svU :: SLA.SpVector Double } deriving (Show)

scaleBetaU :: Double -> BetaU -> BetaU
scaleBetaU x (BetaU b u) = BetaU (LA.scale x b) (SLA.scale x u)

addBetaU :: BetaU -> BetaU -> BetaU
addBetaU (BetaU bA uA) (BetaU bB uB) = BetaU (LA.add bA bB) (uA SLA.^+^ uB)

diffBetaU :: BetaU -> BetaU -> BetaU
diffBetaU (BetaU bA uA) (BetaU bB uB) = BetaU (bA - bB) (uA SLA.^-^ uB)

newtype ZStarMatrix = ZStar { smZS :: SLA.SpMatrix Double }
data RandomEffectCalculated = RandomEffectCalculated { recModelMatrix:: RandomEffectModelMatrix
                                                     , recMkLambda :: (CovarianceVec -> LambdaMatrix)
                                                     }

makeZS :: RandomEffectCalculated -> CovarianceVec -> ZStarMatrix
makeZS reCalc vTh =
  let smZ      = recModelMatrix reCalc
      mkLambda = recMkLambda reCalc
  in  ZStar (smZ SLA.## (mkLambda vTh))


--data LinearPredictorComponents = LPComponents FixedPredictors ZStarMatrix BetaU

data LinearPredictor = LP_Computed EtaVec | LP_ComputeFrom BetaU deriving (Show)


denseLinearPredictor
  :: FixedPredictors -> ZStarMatrix -> LinearPredictor -> EtaVec
denseLinearPredictor _ _ (LP_Computed x) = x
denseLinearPredictor mX zStar (LP_ComputeFrom betaU) =
  mX
    LA.#> (bu_vBeta betaU)
    +     (SD.toDenseVector $ (smZS zStar) SLA.#> bu_svU betaU)

diffLinearPredictor
  :: FixedPredictors
  -> ZStarMatrix
  -> LinearPredictor
  -> LinearPredictor
  -> LinearPredictor
diffLinearPredictor mX zStar (LP_ComputeFrom buX) (LP_ComputeFrom buY) =
  LP_Computed $ denseLinearPredictor mX zStar $ LP_ComputeFrom $ diffBetaU buX
                                                                           buY
diffLinearPredictor mX zStar lpX lpY = LP_Computed (vEtaX - vEtaY)
 where
  vEtaX = denseLinearPredictor mX zStar lpX
  vEtaY = denseLinearPredictor mX zStar lpY

--updateLinearPredictor :: Double -> BetaU -> FixedPredictors -> ZStarMatrix -> LinearPredictor -> LinearPredictor
--updateLinearPredictor x dBetaU mX zStar lp = go (GLM.scaleBetaU x dBetaU) mX zStar lp where
--  go (BetaU vdBeta svdU) mX zStar (LP_Computed vEta) =
--    vdEta = denseLinearPredictor 


{-
type SparseEtaVec = SLA.SpVector Double

sparseLinearPredictor :: FixedPredictors -> ZStarMatrix -> BetaU -> SparseEtaVec
sparseLinearPredictor mX zStar betaU =
  SD.toSparseVector (mX LA.#> bu_vBeta betaU)
    SLA.^+^ ((smZS zStar) SLA.#> bu_svU betaU)
-}
