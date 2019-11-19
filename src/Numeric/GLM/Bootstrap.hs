{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.Bootstrap
  ( module Numeric.GLM.Bootstrap
  , module Polysemy.RandomFu
  , module Polysemy.Async
  , module Statistics.Types
  )
where

import qualified Numeric.LinearAlgebra         as LA
import qualified Data.List                     as L
import           Data.Maybe                     ( catMaybes )
import qualified Data.Text                     as T
import qualified Data.Vector.Storable          as VS
import qualified Data.Vector                   as VB

import qualified Control.Foldl                 as FL

import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )
import qualified Control.Concurrent.KazuraQueue
                                               as KQ
import           GHC.Conc                       ( getNumCapabilities )


import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.MixedModel        as GLM
import qualified Numeric.GLM.Report            as GLM
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH

import qualified Statistics.Types              as S
import           Statistics.Types               ( mkCL )
import qualified Statistics.Distribution       as S
import qualified Statistics.Distribution.StudentT
                                               as S
import qualified Statistics.Distribution.Normal
                                               as S
import qualified Numeric.SpecFunctions         as SF

import qualified Data.Random                   as DR
import qualified Data.Random.Distribution.Normal
                                               as DR
import qualified Data.Random.Distribution.Bernoulli
                                               as DR
import qualified Data.Random.Distribution.Binomial
                                               as DR
import qualified Data.Random.Distribution.Poisson
                                               as DR
import qualified Data.Random.Distribution.Gamma
                                               as DR

import qualified Knit.Effect.Logger            as P
import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Polysemy.RandomFu             as P
import           Polysemy.RandomFu              ( RandomFu
                                                , runRandomIO
                                                )
import qualified Polysemy.Async                as P
import           Polysemy.Async                 ( Async
                                                , asyncToIO
                                                , asyncToIOFinal
                                                )
--import qualified Polysemy.ConstraintAbsorber.MonadRandom
--                                               as P

import           Control.Concurrent            as CC
import           Control.Concurrent.MVar       as CC

-- TODO: handle the sometimes requiring dev parameter in a cleaner way.
generateSample
  :: P.Member P.RandomFu r
  => GLM.ObservationDistribution
  -> Double
  -> Double
  -> P.Sem r Double
generateSample od mu dev = P.sampleRVar d
 where
  d = case od of
    GLM.DNormal    -> DR.normal mu dev
    GLM.DBernoulli -> DR.bernoulli mu
    GLM.DBinomial n ->
      fmap (\x -> realToFrac x / realToFrac n) $ DR.binomial n mu -- binomial observations are as fractions
    GLM.DPoisson -> DR.poisson mu
    GLM.DGamma   -> DR.gamma (1 / dev) (mu * dev) -- Shape & Scale as in deviance


-- NB: assumes the same deviance across all when deviance matters.  Not sure this is always correct. 
generateSamples
  :: P.Member P.RandomFu r
  => GLM.ObservationsDistribution
  -> VS.Vector Double
  -> Double
  -> P.Sem r (VS.Vector Double)
generateSamples osd vMu dev = do
  let n = VS.length vMu
      gs od mu = generateSample od mu dev
      vOD = case osd of
        GLM.Normal      -> VB.replicate n GLM.DNormal
        GLM.Bernoulli   -> VB.replicate n GLM.DBernoulli
        GLM.Binomial vN -> VB.map GLM.DBinomial (VB.convert vN) -- these are boxed
        GLM.Poisson     -> VB.replicate n GLM.DPoisson
        GLM.Gamma       -> VB.replicate n GLM.DGamma
  VS.convert <$> VB.zipWithM gs vOD (VB.convert vMu)


data CholmodStructure = CS_ToCompute GLM.RandomEffectCalculated | CS_Computed GLM.CholmodFactor

cholmodFactor
  :: GLM.EffectsIO r => CholmodStructure -> P.Sem r (GLM.CholmodFactor)
cholmodFactor (CS_ToCompute reCalc) = GLM.cholmodAnalyzeProblem reCalc
cholmodFactor (CS_Computed  cf    ) = return cf

parametricBootstrap
  :: (GLM.EffectsIO r, P.Member P.RandomFu r, P.Member P.Async r)
  => GLM.MinimizeDevianceVerbosity
  -> GLM.DevianceType
  -> GLM.MixedModel b g -- we'll modify this, substuting bootstrapped data
  -> GLM.RandomEffectCalculated
  -> GLM.CholmodFactor
  -> GLM.CovarianceVec -- should be solution from a minimize run
  -> VS.Vector Double -- mu from solution
  -> Double -- deviance from solution, may be ignored
  -> Int -- number of times to bootstrap
  -> Bool
  -> P.Sem r [(GLM.BetaU, VS.Vector Double)] -- beta, u and b
parametricBootstrap mdv dt mm0 reCalc cf thSol vMuSol devSol n doConcurrently =
  do
    let (fpC, fpF, smP) = cf
    factorQueue <- liftIO $ do
      nFactors <- if doConcurrently then getNumCapabilities else return 1
      let factors =
            CS_Computed cf : (replicate (nFactors - 1) $ CS_ToCompute reCalc)
      queue <- KQ.newQueue
      mapM (KQ.writeQueue queue) factors -- load the concurrent Queue with the copied factors
      return queue

{-    
    commonQueue <- do
      nCommons <- if doConcurrently then liftIO getNumCapabilities else return 1
      queue <- liftIO KQ.newQueue
      extraCommons <- traverse (const GLM.cholmodMakeCommon)
        $ replicate (n - 1) ()
      let commons = fpC : extraCommons
      liftIO $ mapM (KQ.writeQueue queue) commons -- load the concurrent Queue with the copied factors
      return queue
-}
    let generateOne =
          generateSamples (GLM.observationDistribution mm0) vMuSol devSol

        solveOne (n, newMM) =
          P.wrapPrefix ("Bootstrap (" <> (T.pack $ show n) <> ")") $ do
            P.logLE P.Diagnostic "Starting..."
            cs                     <- liftIO $ KQ.readQueue factorQueue -- get a factor when available.  Blocks until then
            cf'                    <- cholmodFactor cs
            (_, _, _, betaU, b, _) <- GLM.minimizeDevianceInner mdv
                                                                dt
                                                                newMM
                                                                reCalc
                                                                cf'
                                                                thSol
            liftIO $ KQ.writeQueue factorQueue (CS_Computed cf') -- put the factor back
            P.logLE P.Diagnostic "Finished."
            return (betaU, b)
    newObservations <- mapM (const generateOne) $ replicate n ()
    let newMixedModels =
          fmap (\os -> GLM.changeMMObservations os mm0) newObservations
        seqF = if doConcurrently
          then fmap catMaybes . P.sequenceConcurrently
          else sequence
    res <- seqF $ fmap solveOne $ zip [1 ..] newMixedModels
    let m = length res
    when (m < n)
      $  P.throw
      $  GLM.OtherGLMError
      $  "parametric bootstrap attempted "
      <> (T.pack $ show n)
      <> " refits but only "
      <> (T.pack $ show m)
      <> " succeeded."
    return res

data BootstrapCIType = BCI_Student | BCI_Percentile | BCI_Pivot | BCI_ExpandedPercentile | BCI_Accelerated

bootstrappedConfidence
  :: (Enum b, Bounded b, Ord b, Ord g, Show g, Show b, GLM.Effects r)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe T.Text)
  -> GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> (GLM.BetaU, VS.Vector Double) -- fit result
  -> [(GLM.BetaU, VS.Vector Double)] -- bootstrap fits
  -> BootstrapCIType
  -> S.CL Double
  -> P.Sem r (Double, (Double, Double))
bootstrappedConfidence mm getPredM getLabelM rc ebg (fitBetaU, fitvb) bootstraps ciType cl
  = do
    let predict = GLM.predictFromBetaUB mm getPredM getLabelM rc ebg
    fitX        <- predict fitBetaU fitvb
    bootstrapXs <- traverse (uncurry predict) bootstraps
    ci          <- bootstrapCI ciType fitX bootstrapXs cl
    return (fitX, ci)

bootstrapCI
  :: GLM.Effects r
  => BootstrapCIType
  -> Double
  -> [Double]
  -> S.CL Double
  -> P.Sem r (Double, Double)
bootstrapCI BCI_Student x bxs cl = do
  let sigmaB  = sqrt $ FL.fold FL.variance bxs
      n       = length bxs
      tDist   = S.studentT (realToFrac $ n - 1)
      tFactor = S.quantile tDist (S.significanceLevel cl / 2) -- NB: this is negative 
  return (x + tFactor * sigmaB, x - tFactor * sigmaB)

bootstrapCI BCI_Percentile x bxs cl = do
  let sortedBxs  = L.sort bxs
      n          = length bxs
      indexLower = round $ realToFrac (n - 1) * (S.significanceLevel cl / 2)
      indexUpper = n - 1 - indexLower
  return (sortedBxs !! indexLower, sortedBxs !! indexUpper)

bootstrapCI BCI_Pivot x bxs cl = do
  (l, u) <- bootstrapCI BCI_Percentile x bxs cl
  return (2 * x - u, 2 * x - l)


bootstrapCI BCI_ExpandedPercentile x bxs cl = bootstrapCI BCI_Percentile
                                                          x
                                                          bxs
                                                          cl'
 where
  n        = length bxs
  tDist    = S.studentT (realToFrac $ n - 1)
  tFactor  = S.quantile tDist (S.significanceLevel cl / 2)
  normDist = S.standard
  alpha'   = 2 * S.cumulative
    normDist
    (sqrt (realToFrac n / realToFrac (n - 1)) * tFactor)
  cl' = S.mkCLFromSignificance alpha'

bootstrapCI BCI_Accelerated x bxs cl = do
  let
    n        = length bxs
    nDist    = S.standard
    zAlpha   = S.quantile nDist $ S.significanceLevel cl
    z1mAlpha = S.quantile nDist $ S.confidenceLevel cl

    jxs =
      let s = FL.fold FL.sum bxs
      in  fmap (\x -> (s - x) / (realToFrac $ n - 1)) bxs
    mjx    = FL.fold FL.mean jxs
    numF   = FL.premap (\x -> (mjx - x) ^^ (3 :: Int)) FL.sum
    denomF = fmap (\s -> 6 * (s ** 1.5))
      $ FL.premap (\x -> (mjx - x) ^^ (2 :: Int)) FL.sum
    a            = FL.fold ((/) <$> numF <*> denomF) jxs
    countSmaller = length $ filter (< x) bxs
    alphaZ x = case countSmaller of
      0 -> 0
      m ->
        let fracSmaller = realToFrac m / realToFrac n
            z0          = S.quantile nDist fracSmaller
        in  S.cumulative nDist (z0 + ((z0 + x) / (1 - a * (z0 + x))))
    indexScale = realToFrac (n - 1) -- index varies from 0 to (n-1)
    indexL     = round $ indexScale * (alphaZ zAlpha)
    indexU     = round $ indexScale * (alphaZ z1mAlpha)
    sortedBxs  = L.sort bxs
  P.logLE P.Info
    $  "indexL="
    <> (T.pack $ show indexL)
    <> "; indexU="
    <> (T.pack $ show indexU)
  P.logLE P.Info $ "a=" <> (T.pack $ show a)
  P.logLE P.Info $ "zAlpha=" <> (T.pack $ show zAlpha)
  P.logLE P.Info $ "z1mAlpha=" <> (T.pack $ show z1mAlpha)
  return (sortedBxs !! indexL, sortedBxs !! indexU)

