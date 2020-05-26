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
--import qualified Numeric.MathFunctions.Constants
--                                               as MC

import qualified Statistics.Distribution       as S
import qualified Statistics.Distribution.Binomial
                                               as S
import qualified Statistics.Distribution.Gamma as S
import qualified Statistics.Distribution.Poisson
                                               as S
import qualified Statistics.Distribution.Normal
                                               as S


import qualified Numeric.IEEE as IEEE

data BinomialCount = CountOne Int | CountEach (VS.Vector Int)

data ObservationDistribution = DNormal
                             | DBinomial Int -- vector of total counts, scale free
                             | DBernoulli
                             | DPoisson
                             | DGamma deriving (Show, Eq)

data ObservationsDistribution = Normal | Binomial (VS.Vector Int) | Bernoulli | Poisson | Gamma deriving (Show, Eq)

data LinkFunctionType = IdentityLink | LogisticLink | ExponentialLink | ReciprocalLink deriving (Show, Eq)

data LinkFunction = LinkFunction { link :: Double -> Double -- map from observation to linear predictor, mu -> eta
                                 , invLink :: Double -> Double -- map from linear predictor to observations, eta -> mu
                                 , derivInv :: Double -> Double -- useful in PIRLS algo
                                 }


-- One property of the canonical link (does this define it?), is that
-- it has the property that d(eta)/dmu * (1/var (eta)) = 1
-- this quantity comes up in the PIRLS algo so, whent he link is canonical
-- we simplify.
canonicalLink :: ObservationsDistribution -> LinkFunctionType
canonicalLink Normal       = IdentityLink
canonicalLink (Binomial _) = LogisticLink
canonicalLink Bernoulli    = LogisticLink
canonicalLink Poisson      = ExponentialLink
canonicalLink Gamma        = ReciprocalLink

bounded :: Ord a => a -> a -> a -> a
bounded lb ub = max lb . min ub

unitEpsilonBounded :: (Ord a, IEEE.IEEE a) => a -> a
unitEpsilonBounded = let eps = IEEE.epsilon in bounded eps (1 - eps)

epsilonAwayFrom0 :: (Num a, Ord a, IEEE.IEEE a) => a -> a
epsilonAwayFrom0 x =
  let eps = IEEE.epsilon
  in if x > eps then x
     else if x < -eps then x
          else eps * signum x
               
-- inv is clamped within [epsilon, 1 - epsilon]
linkFunction :: LinkFunctionType -> LinkFunction
linkFunction IdentityLink = LinkFunction id id (const 1)
linkFunction LogisticLink = LinkFunction
  (\mu -> let mu' = unitEpsilonBounded mu in  log (mu' / ( 1 - mu'))) -- eta (mu)
  (\eta -> let y = exp (-eta) in unitEpsilonBounded (1/(1+y))) -- mu (eta) 
  (\eta ->
    let eps :: Double = IEEE.epsilon
        y = exp (-eta)        
        z = 1 + y
    in  max eps (y / (z * z))
  )

-- TODO; do something about x == 0  
linkFunction ExponentialLink = LinkFunction log exp exp
linkFunction ReciprocalLink =
  LinkFunction
  (\x -> -1 / epsilonAwayFrom0 x)
  (\x -> -1 / epsilonAwayFrom0 x)
  (\x -> 1 / epsilonAwayFrom0 (x * x))

data UseLink = UseCanonical | UseOther LinkFunctionType deriving (Show, Eq)
-- notation here:
-- y are the observations
-- mu is the conditional mean of the linear predictor after mapping via the link function
{-
familyWeights
  :: ObservationsDistribution -> LA.Vector Double -> LA.Vector Double
--familyWeights (Binomial vN) vW = VS.zipWith (\w n -> w * realToFrac n) vW vN
familyWeights _ vW = vW
-}
varianceScaledWeights
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
varianceScaledWeights od vW vMu =
  let vVar = scaledVariance od vMu in VS.zipWith (/) vW vVar
--  let vVar = scaledVariance od vMu in VS.zipWith (/) vW vVar


devianceOne :: ObservationDistribution -> Double -> Double -> Double
devianceOne DNormal    y mu = let z = (y - mu) in z * z
devianceOne DBernoulli y mu = -2 * (y * log mu + (1 - y) * log (1 - mu))
devianceOne (DBinomial n) y mu =
  let eps = IEEE.epsilon
      x   = if y < eps
        then negate $ log (1 - (min (1 - eps) mu))
        else if y > (1 - eps)
          then negate $ log $ max eps mu
          else (y * log (y / (max eps mu) + (1 - y) * log ((1 - y) / (1 - (min (1 - eps) mu)))))
  in  2 * (realToFrac n) * x -- should this be multiplied by n?
devianceOne DPoisson y mu =
  let eps = IEEE.epsilon
      x   = if y < eps then mu else (y * log (y / (min eps mu)) - (y - mu))
  in  2 * x
devianceOne DGamma y mu = let mu' = min IEEE.epsilon mu in 2 * ((y - mu) / mu' - log (y / mu'))

deviance
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
deviance od vW vY vMu =
  VS.sum $ zipDist3 od (\od w y mu -> w * devianceOne od y mu) vW vY vMu

-- don't quite get this
devScale
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
devScale od vW vY vMu =
  let vDev = devScaleVec od vW vY vMu
  in  (VS.sum vDev) / (realToFrac $ VS.length vMu)


devScaleVec
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
devScaleVec od vW vY vMu =
  let n   = LA.size vW
      dev = deviance od vW vY vMu
  in  case od of
        Normal -> VS.replicate n (dev / realToFrac n)
        Gamma  -> let sumWgts = VS.sum vW in VS.map (* (dev / sumWgts)) vW
        _      -> VS.replicate n 1

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
  let vDev = devScaleVec od vW vY vMu in -2 * logLikelihood od vW vY vMu vDev


logLikelihoodOne
  :: ObservationDistribution
  -> Double
  -> Double
  -> Double -- dev, might be ignored
  -> Double
logLikelihoodOne DNormal y mu dev =
  S.logDensity (S.normalDistr 0 (sqrt dev)) (y - mu) --let z = (y - mu) in -log (2 * pi * dev)/2 - z*z/(2*dev)
logLikelihoodOne DBernoulli y mu _ =
  S.logProbability (S.binomial 1 mu) (round y)
logLikelihoodOne (DBinomial count) y mu _ =
  S.logProbability (S.binomial count mu) (round $ (realToFrac count) * y)
logLikelihoodOne DPoisson y mu _ = S.logProbability (S.poisson mu) (round y)
logLikelihoodOne DGamma y mu dev =
  S.logDensity (S.gammaDistr (1 / dev) (mu * dev)) y -- NB: this is shape & scale.  mean = mu, var = dev

logLikelihood
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double -- dev, scale parameter, might be ignored
  -> Double
logLikelihood od vW vY vMu vDev = VS.sum $ zipDist4
  od
  (\od w y mu dev -> w * logLikelihoodOne od y mu dev)
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
  -> LA.Vector Double -- dev, may be ignored
  -> Double
aic od vW vY vMu vDev =
  let nParam = case od of
        Normal -> 1
        Gamma  -> 1
        _      -> 0
  in  realToFrac (2 * nParam) - 2 * logLikelihood od vW vY vMu vDev

aicR
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double -- may be ignored. When not ignored is assumed to be the sum of deviance residuals.  Like R.
  -> Double
aicR od vW vY vMu devResid =
  let n      = VS.length vY
      vDev   = VS.replicate n (devResid / realToFrac n)
      nParam = case od of
        Normal -> 1
        Gamma  -> 1
        _      -> 0
  in  realToFrac (2 * nParam) - 2 * logLikelihood od vW vY vMu vDev


-- this needs a better name.  Which means I need to understand it.
-- scaled variance?
scaledVarianceOne :: ObservationDistribution -> Double -> Double
scaledVarianceOne DNormal       _ = 1
scaledVarianceOne DBernoulli    x = x * (1 - x)
scaledVarianceOne (DBinomial n) x = x * (1 - x) / realToFrac n -- not sure about n here
scaledVarianceOne DPoisson      x = x
scaledVarianceOne DGamma        x = x * x

scaledVariance
  :: ObservationsDistribution -> LA.Vector Double -> LA.Vector Double
scaledVariance (Binomial vN) vMu =
  VS.zipWith (\n mu -> scaledVarianceOne (DBinomial n) mu) vN vMu
scaledVariance Normal    vMu = VS.map (scaledVarianceOne DNormal) vMu
scaledVariance Bernoulli vMu = VS.map (scaledVarianceOne DBernoulli) vMu
scaledVariance Poisson   vMu = VS.map (scaledVarianceOne DPoisson) vMu
scaledVariance Gamma     vMu = VS.map (scaledVarianceOne DGamma) vMu


canonicaldEtadMuOverVarOne :: ObservationDistribution -> Double -> Double
canonicaldEtadMuOverVarOne DNormal _ = 0
canonicaldEtadMuOverVarOne DBernoulli eta = let y = exp(-eta) in sqrt (y) / (1+y)
canonicaldEtadMuOverVarOne (DBinomial n) eta = let y = exp(-eta) in sqrt (realToFrac n * y) / (1+y)
canonicaldEtadMuOverVarOne DPoisson eta = exp (eta / 2)
canonicaldEtadMuOverVarOne DGamma eta = 1 / eta

canonicaldEtadMuOverVar :: ObservationsDistribution -> LA.Vector Double -> LA.Vector Double
canonicaldEtadMuOverVar Normal vEta = VS.map (canonicaldEtadMuOverVarOne DNormal) vEta
canonicaldEtadMuOverVar Bernoulli vEta = VS.map (canonicaldEtadMuOverVarOne DBernoulli) vEta
canonicaldEtadMuOverVar (Binomial vN) vEta = VS.zipWith (\n eta -> canonicaldEtadMuOverVarOne (DBinomial n) eta) vN vEta
canonicaldEtadMuOverVar Poisson vEta = VS.map (canonicaldEtadMuOverVarOne DPoisson) vEta
canonicaldEtadMuOverVar Gamma vEta = VS.map (canonicaldEtadMuOverVarOne DGamma) vEta



-- Numeric helpers 

zipDist3
  :: (VS.Storable x1, VS.Storable x2, VS.Storable x3, VS.Storable x4)
  => ObservationsDistribution
  -> (ObservationDistribution -> x1 -> x2 -> x3 -> x4)
  -> VS.Vector x1
  -> VS.Vector x2
  -> VS.Vector x3
  -> VS.Vector x4
zipDist3 Normal        f = VS.zipWith3 (f DNormal)
zipDist3 Bernoulli     f = VS.zipWith3 (f DBernoulli)
zipDist3 (Binomial vN) f = VS.zipWith4 (\n -> f (DBinomial n)) vN
zipDist3 Poisson       f = VS.zipWith3 (f DPoisson)
zipDist3 Gamma         f = VS.zipWith3 (f DGamma)

zipDist4
  :: ( VS.Storable x1
     , VS.Storable x2
     , VS.Storable x3
     , VS.Storable x4
     , VS.Storable x5
     )
  => ObservationsDistribution
  -> (ObservationDistribution -> x1 -> x2 -> x3 -> x4 -> x5)
  -> LA.Vector x1
  -> LA.Vector x2
  -> LA.Vector x3
  -> LA.Vector x4
  -> LA.Vector x5
zipDist4 Normal        f = VS.zipWith4 (f DNormal)
zipDist4 Bernoulli     f = VS.zipWith4 (f DBernoulli)
zipDist4 (Binomial vN) f = VS.zipWith5 (\n -> f (DBinomial n)) vN
zipDist4 Poisson       f = VS.zipWith4 (f DPoisson)
zipDist4 Gamma         f = VS.zipWith4 (f DGamma)



{-
yLny :: Double -> Double
yLny x =
  let eps = 1e-12 -- do something type/machine appropriate here
  in case x < eps of
    True -> eps
-}
