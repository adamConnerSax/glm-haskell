{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.GLM.Predict where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.MixedModel        as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Numeric.SparseDenseConversions
                                               as SD

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Knit.Effect.Logger            as P
import qualified Control.Foldl                 as FL
--import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Map                      as M
import           Data.Maybe                     ( catMaybes )
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.List                     as L
--import qualified Colonnade                     as C

import qualified Numeric.LinearAlgebra.Class   as SLA

import qualified Numeric.LinearAlgebra         as LA

import qualified Statistics.Types              as S
import qualified Statistics.Distribution       as S
import qualified Statistics.Distribution.StudentT
                                               as S

import qualified Data.Text                     as T
import qualified Data.Vector.Storable          as VS
import qualified Data.Vector.Split             as VS


{-
n rows
p fixed effects
l effect levels
n_l number of categories in level l
q_l number of (random) effects in level l
q = sum (n_l * q_l)
X is n x p
Z is n x q
zTz is q x q, this is the part that will be sparse.
-}

eitherToSem :: GLM.Effects r => Either T.Text a -> P.Sem r a
eitherToSem = either (P.throw . GLM.OtherGLMError) return

fixedEffectParameters
  :: GLM.PredictorC b
  => GLM.MixedModel b g --GLM.FixedEffects b
  -> GLM.BetaU
  -> GLM.FixedEffectParameters b
fixedEffectParameters mm (GLM.BetaU vBeta _) =
  let means = vBeta
      fe    = GLM.rmsFixedEffects $ GLM.regressionModelSpec mm
  in  GLM.FixedEffectParameters fe means

fixedEffectStatistics
  :: GLM.PredictorC b
  => GLM.MixedModel b g --GLM.FixedEffects b
  -> Double
  -> GLM.CholeskySolutions
  -> GLM.BetaU
  -> GLM.FixedEffectStatistics b
fixedEffectStatistics mm sigma2 (GLM.CholeskySolutions _ _ mRX) betaU =
  let fep       = fixedEffectParameters mm betaU
      effSigma2 = case GLM.observationsDistribution mm of
        GLM.Normal -> sigma2
        _          -> 1
      covs = LA.scale effSigma2 $ (LA.inv mRX) LA.<> (LA.inv $ LA.tr mRX)
  in  GLM.FixedEffectStatistics fep covs

groupOffsets
  :: (Ord g, Show g, Show b, GLM.Effects r)
  => GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> P.Sem r [Int]
groupOffsets rc epg = do
  let groups = IS.members $ GLM.groupIndices rc
      groupOffset group = do
        size     <- GLM.groupSize rc group
        nEffects <- IS.size <$> GLM.groupEffects epg group
        return $ size * nEffects
  eitherToSem $ FL.prescan FL.sum <$> traverse groupOffset groups

effectParametersByGroup
  :: (Ord g, Show g, Show b, GLM.Effects r)
  => GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> LA.Vector Double
  -> P.Sem r (GLM.EffectParametersByGroup g b)
effectParametersByGroup rc epg vb = do
  offsets <- groupOffsets rc epg
  let
    groups = IS.members $ GLM.groupIndices rc
    ep (group, offset) = do
      size    <- GLM.groupSize rc group
      effects <- GLM.groupEffects epg group
      let nEffects = IS.size effects
          parameterMatrix =
            LA.fromRows
              $ VS.chunksOf nEffects
              $ VS.take (size * nEffects)
              $ VS.drop offset vb
      return (group, GLM.EffectParameters effects parameterMatrix)
  eitherToSem $ fmap M.fromList $ traverse ep $ zip groups offsets


effectCovariancesByGroup
  :: ( Ord g
     , Show g
     , Ord (GLM.GroupKey g)
     , Enum b
     , Bounded b
     , Show b
     , GLM.Effects r
     )
  => GLM.EffectsByGroup g b
  -> GLM.MixedModel b g
  -> Double
  -> LA.Vector Double
  -> P.Sem r (GLM.EffectCovariancesByGroup g b)
effectCovariancesByGroup ebg mm sigma2 vTh = do
  let
    effSigma2 = case GLM.observationsDistribution mm of
      GLM.Normal -> sigma2
      _          -> 1
    groups = M.keys ebg
    nElements n = n * (n + 1) `div` 2
    groupOffset group = nElements . IS.size <$> GLM.groupEffects ebg group
    ltIndices n = [ (r, c) | c <- [0 .. (n - 1)], r <- [c .. (n - 1)] ]
    ltAssocs n v = zip (ltIndices n) (VS.toList v)
    tF ((r, c), x) = if r > c then Just ((c, r), x) else Nothing
    tAssocs l = catMaybes $ fmap tF l
    makeLT n v = let l = ltAssocs n v in LA.assoc (n, n) 0 $ l ++ (tAssocs l)
  offsets <- eitherToSem $ FL.prescan FL.sum <$> traverse groupOffset groups
  let ecs (group, offset) = do
        effects <- GLM.groupEffects ebg group
        let nEffects = IS.size effects
            cm =
              LA.scale effSigma2
                $ makeLT nEffects
                $ VS.take (nElements nEffects)
                $ VS.drop offset vTh
        return (group, GLM.GroupEffectCovariances effects cm)
  eitherToSem $ fmap M.fromList $ traverse ecs $ zip groups offsets

predictFromBetaUB
  :: (GLM.PredictorC b, GLM.GroupC g, GLM.Effects r)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe (GLM.GroupKey g))
  -> GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> GLM.BetaU
  -> VS.Vector Double -- b
  -> P.Sem r Double
predictFromBetaUB mm getPredictorM getLabelM rc ebg betaU vb = do
  let fep = fixedEffectParameters mm betaU
  epg <- effectParametersByGroup rc ebg vb
  fst <$> predict mm getPredictorM getLabelM fep epg rc

predict
  :: (GLM.PredictorC b, GLM.GroupC g, GLM.Effects r)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe (GLM.GroupKey g))
  -> GLM.FixedEffectParameters b
  -> GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> P.Sem r (Double, Double) -- (prediction, link (prediction))
predict mm getPredictorM getLabelM (GLM.FixedEffectParameters fe vFE) ebg rowClassifier
  = do
    let
      invLinkF = GLM.invLink $ GLM.linkFunction $ GLM.linkFunctionType mm
      applyEffectE msgF b e = case b of
        GLM.Intercept -> Right e
        GLM.Predictor b ->
          maybe (Left $ msgF b) Right $ fmap (e *) $ getPredictorM b
      applyEffectsE msgF ies v =
        fmap (FL.fold FL.sum)
          $ traverse
              (\(index, effect) ->
                applyEffectE msgF effect (v `LA.atIndex` index)
              )
          $ IS.toIndexedList ies
      groups = M.keys ebg
      grpMsg g =
        "Failed to find group label for \"" <> (T.pack $ show g) <> "\""
      effMsgGroup g b =
        "Failed to find predictor for effect \""
          <> (T.pack $ show b)
          <> "\" which is modeled in group \""
          <> (T.pack $ show g)
          <> "\".  This is a very weird error and should not happen. There should never be group effects which are not also fixed effects."
      applyGroupEffects group = do
        case getLabelM group of
          Nothing    -> return 0 -- if we don't put in a group label, we assume 0 group effect
          Just label -> do
            labelIndex <- GLM.labelIndex rowClassifier group label
            (GLM.EffectParameters ies mEP) <- GLM.effectParameters group ebg
            applyEffectsE (effMsgGroup group)
                          ies
                          (LA.toRows mEP L.!! labelIndex)
      effMsgFixed b =
        "No predictor for fixed effect \"" <> (T.pack $ show b) <> "\""
    fixedEffects <- eitherToSem
      $ applyEffectsE effMsgFixed (GLM.indexedFixedEffectSet fe) vFE
    groupEffects <- eitherToSem $ traverse applyGroupEffects $ M.keys ebg
    let totalEffects = fixedEffects + FL.fold FL.sum groupEffects
    return (invLinkF totalEffects, totalEffects)

-- predict by constructing vCoeffBeta and vCoeffb
fixedEffectsCoefficientVector
  :: (GLM.Effects r, GLM.PredictorC b)
  => (b -> Maybe Double)
  -> GLM.FixedEffects b
  -> P.Sem r (VS.Vector Double)
fixedEffectsCoefficientVector getPredictorM fes = case fes of
  GLM.InterceptOnly    -> return $ VS.fromList [1]
  GLM.FixedEffects ies -> effectsCoefficientVector getPredictorM ies

effectsCoefficientVector
  :: (GLM.Effects r, GLM.PredictorC b)
  => (b -> Maybe Double)
  -> GLM.IndexedEffectSet b
  -> P.Sem r (VS.Vector Double)
effectsCoefficientVector getPredictorM ies = eitherToSem $ do
  let
    indices = [0 .. (IS.size ies - 1)]
    getCol n =
      maybe
          (  Left
          $  "couldn't find index="
          <> (T.pack $ show n)
          <> " in fixed effect index."
          )
          Right
        $            ies
        `IS.atIndex` n
    getCoeff x = case x of
      GLM.Intercept -> Right 1
      GLM.Predictor y ->
        maybe (Left $ "Couldn't get coefficient for " <> (T.pack $ show y))
              Right
          $ getPredictorM y
    getIndexedCoeff n = getCol n >>= getCoeff
  VS.fromList <$> traverse getIndexedCoeff indices


indexedGroupCoefficients
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => (b -> Maybe Double)
  -> (g -> Maybe (GLM.GroupKey g))
  -> GLM.EffectsByGroup g b
  -> GLM.RowClassifier g
  -> P.Sem r [(Int, Double)]
indexedGroupCoefficients getPredictorM getLabelM ebg rc = do
  gOffsets <- groupOffsets rc ebg
  --P.logLE P.Info $ "offsets=" <> (T.pack $ show gOffsets)
  let vecLength = FL.fold FL.sum gOffsets
  -- get a list of groups    
  let
    gIS    = GLM.groupIndices rc
    groups = IS.members gIS
    groupIndex g =
      eitherToSem
        $ maybe
            (  Left
            $  "Couldn't find "
            <> (T.pack $ show g)
            <> " in "
            <> (T.pack $ show gIS)
            )
            Right
        $ IS.index gIS g

    getIndexedCoeffs g = do
      groupEffectsIndex <-
        eitherToSem
        $ maybe
            (  Left
            $  "Couldn't find \""
            <> (T.pack $ show g)
            <> " in "
            <> (T.pack $ show ebg)
            )
            Right
        $ M.lookup g ebg
      let nEffects = IS.size groupEffectsIndex
      gIndex <- groupIndex g
      let groupOffset = gOffsets !! gIndex
      case getLabelM g of
        Nothing -> P.logLE P.Info ("No entry for group=" <> (T.pack $ show g))
          >> return []
        Just key -> case GLM.categoryNumberFromKey rc key g of
          Nothing ->
            P.throw
              $  GLM.OtherGLMError
              $  "Failed to find \""
              <> (T.pack $ show key)
              <> "\" in group=\""
              <> (T.pack $ show g)
              <> "\". This is BAD."
          Just catOffset -> do
            effectCoeffs <- effectsCoefficientVector getPredictorM
                                                     groupEffectsIndex
            let vecIndices = fmap
                  (\n -> groupOffset + (catOffset * nEffects) + n)
                  [0 .. (VS.length effectCoeffs - 1)]
            return $ zip vecIndices (VS.toList effectCoeffs)
  indexedSparseVecEntries <- concat <$> traverse getIndexedCoeffs groups
  --P.logLE P.Diagnostic $ T.pack $ show indexedSparseVecEntries
  return indexedSparseVecEntries

predict'
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe (GLM.GroupKey g))
  -> GLM.EffectsByGroup g b
  -> GLM.RowClassifier g
  -> GLM.BetaU
  -> LA.Vector Double -- vb
  -> P.Sem r (Double, Double) -- (prediction, link (prediction))
predict' mm getPredictorM getLabelM ebg rowClassifier betaU vb = do
  let fes      = GLM.rmsFixedEffects $ GLM.regressionModelSpec mm
      invLinkF = GLM.invLink $ GLM.linkFunction $ GLM.linkFunctionType mm
  vFixedCoeffs  <- fixedEffectsCoefficientVector getPredictorM fes
--  P.logLE P.Diagnostic $ "vFixed=" <> (T.pack $ show vFixedCoeffs)
  svGroupCoeffs <-
    SLA.fromListSV (VS.length vb)
      <$> indexedGroupCoefficients getPredictorM getLabelM ebg rowClassifier
--  P.logLE P.Diagnostic $ "vGroup=" <> (T.pack $ show svGroupCoeffs)
  let svb = SD.toSparseVector vb
      totalEffects =
        (vFixedCoeffs LA.<.> (GLM.bu_vBeta betaU)) + (svGroupCoeffs SLA.<.> svb)
  return (invLinkF totalEffects, totalEffects)

-- We can compute a very crude CI by **assuming** that the likelihood function is approximately quadratic.
-- We compute the covariance matrices for the beta and the bs
-- The betas and bs are uncorrelated (by construction??)
-- So we compute the variance, v, for the linear predictor given predictors and group(s).
-- From that we can compute a studentized CI via: eta +/- t(alpha/2)*sqrt(v)
-- if the inverse-link is F: mu = F(eta) with CI [F(eta - t(alpha/2) * sqrt(v)), F(eta + t(alpha/2) * sqrt(v))]
-- It's *deeply* approximate.  For better answers, Bootstrap!
-- TODO: Wald Test??


-- Here we compute the conditional covariances of the random/group effects
-- This is expensive, but not as expensive as bootstrapping
type ConditionalCovarianceMatrix = SLA.SpMatrix Double

conditionalCovariances
  :: GLM.EffectsIO r
  => GLM.MixedModel b g
  -> GLM.CholmodFactor
  -> GLM.RandomEffectCalculated
  -> GLM.CovarianceVec
  -> GLM.BetaU
  -> P.Sem r ConditionalCovarianceMatrix
conditionalCovariances mm cf reCalc vTh betaU = do
  let (cholmodC, cholmodF, _) = cf
      mX = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      zStar                   = GLM.makeZS reCalc vTh
      smLambda                = (GLM.recMkLambda reCalc) vTh
      smLambdat               = SLA.transpose smLambda
      vEta = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom betaU)
      lf        = GLM.linkFunction $ GLM.linkFunctionType mm
      vMu       = VS.map (GLM.invLink lf) vEta
      vM        = VS.map (GLM.derivInv lf) vEta
      vW        = GLM.weights mm
      vVSW = GLM.varianceScaledWeights (GLM.observationsDistribution mm) vW vMu
  (smU, _) <- GLM.spUV mm zStar vVSW vM 
  let lf                      = GLM.linkFunction $ GLM.linkFunctionType mm
      sigma2                  = GLM.devScale (GLM.observationsDistribution mm)
                                             (GLM.weights mm)
                                             (GLM.observations mm)
                                             (VS.map (GLM.invLink lf) vEta)
  -- reFactorize because the factor gets trashed when we extract L at the end of chomod
      cfs = CH.FactorizeAtAPlusBetaI 1
  liftIO $ CH.spMatrixFactorizeP cholmodC
                                 cholmodF
                                 cfs
                                 CH.UnSymmetric
                                 (SLA.transpose smU)
  smLLtInvLambda <- liftIO
    $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_A smLambdat
  let smCondVar = smLLtInvLambda SLA.## smLambda
  return $ sigma2 SLA..* smCondVar


predictWithCondVarCI
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe (GLM.GroupKey g))
  -> GLM.EffectsByGroup g b
  -> GLM.RowClassifier g
  -> GLM.BetaU
  -> LA.Vector Double -- vb
  -> S.CL Double
  -> LA.Matrix Double -- Cov_Beta  
  -> ConditionalCovarianceMatrix --Cov_U
  -> P.Sem r (Double, (Double, Double)) -- (prediction, lower, upper)
predictWithCondVarCI mm getPredictorM getLabelM ebg rowClassifier betaU vb cl mCovBeta smCondVar
  = do
    let fes      = GLM.rmsFixedEffects $ GLM.regressionModelSpec mm
        invLinkF = GLM.invLink $ GLM.linkFunction $ GLM.linkFunctionType mm
        n        = VS.length $ GLM.observations mm
    vFixedCoeffs  <- fixedEffectsCoefficientVector getPredictorM fes
  --  P.logLE P.Diagnostic $ "vFixed=" <> (T.pack $ show vFixedCoeffs)
    svGroupCoeffs <-
      SLA.fromListSV (VS.length vb)
        <$> indexedGroupCoefficients getPredictorM getLabelM ebg rowClassifier
  --  P.logLE P.Diagnostic $ "vGroup=" <> (T.pack $ show svGroupCoeffs)
    let svb = SD.toSparseVector vb
        totalEffects =
          (vFixedCoeffs LA.<.> (GLM.bu_vBeta betaU))
            + (svGroupCoeffs SLA.<.> svb)
        varianceTE =
          ((mCovBeta LA.#> vFixedCoeffs) LA.<.> vFixedCoeffs)
            + ((smCondVar SLA.#> svGroupCoeffs) SLA.<.> svGroupCoeffs)
        tDist   = S.studentT (realToFrac $ n - 1)
        tFactor = S.quantile tDist (S.significanceLevel cl / 2)  -- This is negative
        lower   = totalEffects + (tFactor * sqrt varianceTE)
        upper   = totalEffects - (tFactor * sqrt varianceTE)
    return (invLinkF totalEffects, (invLinkF lower, invLinkF upper))


