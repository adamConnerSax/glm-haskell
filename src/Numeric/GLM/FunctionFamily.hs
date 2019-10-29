{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.FunctionFamily where

import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Vector.Storable          as VS
import qualified Numeric.MathFunctions.Constants
                                               as MC

import qualified Statistics.Distribution       as S
import qualified Statistics.Distribution.Binomial
                                               as S
import qualified Statistics.Distribution.Gamma as S
import qualified Statistics.Distribution.Poisson
                                               as S
import qualified Statistics.Distribution.Normal
                                               as S

--data FittedScale = FittedScale Double | Unfitted deriving (Show, Eq)-- some distributions have scale parameters determined by the fit

data BinomialCount = CountOne Int | CountEach (VS.Vector Int)

data ObservationDistribution = DNormal
                             | DBinomial Int -- vector of total counts, scale free
                             | DBernoulli
                             | DPoisson
                             | DGamma deriving (Show, Eq)

data ObservationsDistribution = Normal | Binomial (VS.Vector Int) | Bernoulli | Poisson | Gamma

data LinkFunctionType = IdentityLink | LogisticLink | ExponentialLink deriving (Show, Eq)

data LinkFunction = LinkFunction { link :: Double -> Double -- map from observation to linear predictor
                                 , invLink :: Double -> Double -- map from linear predictor to obervations
                                 , derivInv :: Double -> Double -- useful in PIRLS algo
                                 }

canonicalLink :: ObservationsDistribution -> LinkFunctionType
canonicalLink Normal       = IdentityLink
canonicalLink (Binomial _) = LogisticLink
canonicalLink Bernoulli    = LogisticLink
canonicalLink Poisson      = ExponentialLink

linkFunction :: LinkFunctionType -> LinkFunction
linkFunction IdentityLink = LinkFunction id id (const 1)
linkFunction LogisticLink = LinkFunction
  (\x -> log (x / (1 - x)))
  (\x -> let y = exp (-x) in 1 / (1 + y))
  (\x ->
    let y = exp (-x)
        z = 1 + y
    in  y / (z * z)
  )
linkFunction ExponentialLink = LinkFunction log exp exp


data UseLink = UseCanonical | UseOther LinkFunctionType deriving (Show, Eq)

-- notation here:
-- y are the observations
-- mu is the conditional mean of the linear predictor after mapping via the link function
familyWeights
  :: ObservationsDistribution -> LA.Vector Double -> LA.Vector Double
familyWeights (Binomial vN) vW =
  let f w n = w * realToFrac n in VS.zipWith f vW vN
familyWeights _ vW = vW

varianceScaledWeights
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
varianceScaledWeights od vW vMu =
  let vVar = scaledVariance od vMu in VS.zipWith (/) (familyWeights od vW) vVar


devianceOne :: ObservationDistribution -> Double -> Double -> Double
devianceOne DNormal    y mu = (y - mu)
devianceOne DBernoulli y mu = -2 * (y * log mu + (1 - y) * log (1 - mu))
devianceOne (DBinomial n) y mu =
  let eps = 1e-12
      x   = if y < eps
        then negate $ log (1 - mu)
        else if (1 - y) < eps
          then negate $ log mu
          else (y * log (y / mu) + (1 - y) * log ((1 - y) / (1 - mu)))
  in  2 * (realToFrac n) * x
devianceOne DPoisson y mu =
  let eps = 1e-12
      x   = if y < eps then mu else (y * log (y / mu) - (y - mu))
  in  2 * x
devianceOne DGamma y mu = 2 * ((y - mu) / mu - log (y / mu))

deviance
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
deviance od vW vY vMu = case od of
  Normal ->
    VS.sum $ VS.zipWith3 (\w y mu -> w * devianceOne DNormal y mu) vW vY vMu
  Bernoulli -> VS.sum
    $ VS.zipWith3 (\w y mu -> w * devianceOne DBernoulli y mu) vW vY vMu
  Binomial vN -> VS.sum $ VS.zipWith4
    (\n w y mu -> w * devianceOne (DBinomial n) y mu)
    vN
    vW
    vY
    vMu
  Poisson ->
    VS.sum $ VS.zipWith3 (\w y mu -> w * devianceOne DPoisson y mu) vW vY vMu
  Gamma -> VS.sum $ VS.zipWith3 (\w y mu -> devianceOne DGamma y mu) vW vY vMu

-- don't quite get this
devScale
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
devScale od vW vY vMu =
  let n   = LA.size vW
      dev = deviance od vW vY vMu
  in  case od of
        Normal -> VS.replicate n (dev / realToFrac n)
        Gamma  -> let sumWgts = VS.sum vW in VS.map (* (dev / sumWgts)) vW
        _      -> VS.replicate (LA.size vMu) 1

devianceCondRel
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
devianceCondRel = deviance

devianceCondAbs
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
devianceCondAbs od vW vY vMu =
  let vDev = devScale od vW vY vMu in -2 * logLikelihood od vW vY vMu vDev


logLikelihoodOne
  :: ObservationDistribution
  -> Double
  -> Double
  -> Double -- dev, might be igbored
  -> Double
logLikelihoodOne DNormal y mu dev =
  S.logDensity (S.normalDistr 0 (sqrt dev)) (y - mu) --let z = (y - mu) in -log (2 * pi * dev)/2 - z*z/(2*dev)
logLikelihoodOne DBernoulli y mu _ =
  S.logProbability (S.binomial 1 mu) (round y)
logLikelihoodOne (DBinomial count) y mu _ =
  S.logProbability (S.binomial count mu) (round $ realToFrac count * y)
logLikelihoodOne DPoisson y mu _ = S.logProbability (S.poisson mu) (round y)
logLikelihoodOne DGamma y mu dev =
  S.logDensity (S.gammaDistr (1 / dev) (mu * dev)) y

logLikelihood
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double -- dev, might be ignored
  -> Double
logLikelihood Normal vW vY vMu vDev = VS.sum $ VS.zipWith4
  (\w y mu d -> w * logLikelihoodOne DNormal y mu d)
  vW
  vY
  vMu
  vDev
logLikelihood Bernoulli vW vY vMu vDev = VS.sum $ VS.zipWith4
  (\w y mu d -> w * logLikelihoodOne DBernoulli y mu d)
  vW
  vY
  vMu
  vDev
logLikelihood (Binomial vN) vW vY vMu vDev = VS.sum $ VS.zipWith5
  (\n w y mu d -> w * logLikelihoodOne (DBinomial n) y mu d)
  vN
  vW
  vY
  vMu
  vDev
logLikelihood Poisson vW vY vMu vDev = VS.sum $ VS.zipWith4
  (\w y mu d -> w * logLikelihoodOne DPoisson y mu d)
  vW
  vY
  vMu
  vDev
logLikelihood Gamma vW vY vMu vDev = VS.sum $ VS.zipWith4
  (\w y mu d -> w * logLikelihoodOne DGamma y mu d)
  vW
  vY
  vMu
  vDev

aicOne
  :: ObservationDistribution
  -> Double
  -> Double
  -> Double -- for scale free dists, we ignore this argument
  -> Double
aicOne od y mu dev =
  let nParam = case od of
        DNormal -> 1
        DGamma  -> 1
        _       -> 0
  in  2 * (nParam - logLikelihoodOne od y mu dev)

aic
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double -- may be ignored
  -> Double
aic od vW vY vMu vDev =
  let nParam = case od of
        Normal -> 1
        Gamma  -> 1
        _      -> 0
  in  2 * (1 - logLikelihood od vW vY vMu vDev)

-- this needs a better name.  Which means I need to understand it.
-- scaled variance?
scaledVarianceOne :: ObservationDistribution -> Double -> Double
scaledVarianceOne DNormal       _ = 1
scaledVarianceOne DBernoulli    x = x * (1 - x)
scaledVarianceOne (DBinomial n) x = realToFrac n * x * (1 - x) -- not sure about n here
scaledVarianceOne DPoisson      x = x
scaledVarianceOne DGamma        x = x * x

scaledVariance
  :: ObservationsDistribution -> LA.Vector Double -> LA.Vector Double
scaledVariance Normal    vMu = VS.map (scaledVarianceOne DNormal) vMu
scaledVariance Bernoulli vMu = VS.map (scaledVarianceOne DBernoulli) vMu
scaledVariance (Binomial vN) vMu =
  VS.zipWith (\n mu -> scaledVarianceOne (DBinomial n) mu) vN vMu
scaledVariance Poisson vMu = VS.map (scaledVarianceOne DPoisson) vMu
scaledVariance Gamma   vMu = VS.map (scaledVarianceOne DGamma) vMu


-- Numeric helpers 
{-
yLny :: Double -> Double
yLny x =
  let eps = 1e-12 -- do something type/machine appropriate here
  in case x < eps of
    True -> eps
-}
