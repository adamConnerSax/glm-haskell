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

--data FittedScale = FittedScale Double | Unfitted deriving (Show, Eq)-- some distributions have scale parameters determined by the fit

data BinomialCount = CountOne Int | CountEach (VS.Vector Int)

data ObservationDistr = DNormal
                      | DBinomial Int -- vector of total counts, scale free
                      | DBernoulli
                      | DPoisson
                      | DGamma

data ObservationsDistribution = Normal | Binomial (VS.Vector Int) | Bernoulli | Poisson | Gamma

data LinkFunctionType = IdentityLink | LogisticLink | ExponentialLink deriving (Show, Eq)

data LinkFunction = LinkFunction { link :: Double -> Double -- map from observation to linear predictor
                                 , invLink :: Double -> Double -- map from linear predictor to obervations
                                 , derivInv :: Double -> Double -- useful in PIRLS algo
                                 }

canonicalLink :: ObservationsDistribution -> LinkFunctionType
canonicalLink (Normal   _) = IdentityLink
canonicalLink (Binomial _) = LogisticLink
canonicalLink Bernoulli    = LogisticLink
canonicalLink (Poisson _)  = ExponentialLink

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
  :: ObservationsDistribution -> LA.Vector Double -> Maybe LA.Vector Double
familyWeights (Binomial vN) vW =
  let f w n = w * realToFrac n in Just (VS.zipWith f vW vN)
familyWeights _ vW = vW

varianceScaledWeights
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
varianceScaledWeights od vW vMu =
  let vVar = variance od vMu in VS.zipWith (/) (familyWeights od vW) vVar


devianceOne :: ObservationDistribution -> Double -> Double -> Double -> Double
devianceOne DNormal    w y mu = w * (y - mu)
devianceOne DBernoulli w y mu = -w * 2 * (y * log mu + (1 - y) * log (1 - mu))
devianceOne (DBinomial n) w y mu =
  let eps = 1e-12
      x   = if y < eps
        then negate $ log (1 - mu)
        else if (1 - y) < eps
          then negate $ log mu
          else (y * log (y / mu) + (1 - y) * log ((1 - y) / (1 - mu)))
  in  2 * n * w * x
devianceOne DPoisson w y mu =
  let eps = 1e-12
      x   = if y < eps then mu else (y * log (y / mu) - (y - mu))
  in  2 * x
devianceOne DGamma w y mu = 2 * w * ((y - mu) / mu - log (y / mu))

deviance
  :: ObservationDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
deviance od vW vY vMu = case od of
  Normal      -> VS.zip3 (devianceOne DNormal) vW vY vMu
  Bernoulli   -> VS.zip3 (devianceOne DBernoulli) vW vY vMu
  Binomial vN -> VS.zip4 (\n -> devianceOne (DBinomial n)) vN vW vY vMu
  Poisson     -> VS.zip3 (devianceOne DPoisson) vW vY vMu
  Gamma       -> VS.zip3 (devianceOne DGamma) vW vY vMu
{-    
  = let
      eps     = 1e-12
      weights = familyWeights od vW
      f y mu = case od of
        Normal   -> (y - mu) ** 2
        Bernoulli    -> -2 * (y * log mu + (1 - y) * log (1 - mu))
        (Binomial _) -> if y < eps
          then -2 * log (1 - mu) -- y = 0
          else if (1 - y) < eps
            then -2 * log mu  -- y = 1
            else 2 * (y * log (y / mu) + (1 - y) * log ((1 - y) / (1 - mu)))
--            else -2 * (y * log mu + (1 - y) * log (1 - mu)) -- 
        (Poisson _) ->
          if y > eps then 2 * (y * log (y / mu) - (y - mu)) else 2 * mu
        (Gamma _) -> 2 * ((y - mu) / mu - log (y / mu))
      g h w y mu = w * h y mu
    in
      VS.sum $ VS.zipWith3 (g f) weights vY vMu
-}


-- don't quite get this
devScale
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Maybe (LA.Vector Double)
devScale od vW vY vMu =
  let n   = LA.dim vW
      dev = deviance od vW vY vMu
  in  case od of
        Normal -> Just $ VS.fromList $ L.replicate n (dev / realToFrac n)
        Gamma ->
          let sumWgts = VS.sum vW in Just $ VS.map (* (dev / sumWgts)) vW
        _ -> Nothing

devianceAbsOne :: ObservationDistribution

logLikelihood
  :: ObservationsDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Maybe Double
logLikelihood od vW vY vMu =
  let dof = case od of
        (Normal _) -> 1
        (Gamma  _) -> 1
        _          -> 0
  in  2 * dof - aic od vW vY vMu


aicOne
  :: ObservationDistribution
  -> Double
  -> Double
  -> Double
  -> Double -- for scale free dists, we ignore this argument
  -> Double
aicOne DNormal w y mu dev = 2 - (log w - (log (2 * pi * dev) + 1)) -- FIX FIX FIX
aicOne DBernouilli w y mu _ =
  -2 * w * S.logProbability (S.Binomial 1 mu) (round y)
aicOne (DBinomial n) w y mu _ =
  -2 * w * S.logProbability (S.Binomial counts mu) (round $ y * n)
aicOne DPoisson w y mu _ = -2 * w * S.logProbability (S.Poisson mu) (round y)
aicOne DGamma w y mu dev =
  2 - 2 * w * S.logDensity (S.gammaDistr (w / dev) (mu * dev / w)) y


aic
  :: ObservationDistribution
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double -- may be ignored
  -> Double
aic od vW vY vMu vDev = case od of
  Normal -> 2 + VS.sum $ VS.zipWith4
    (\w y mu dev -> (aicOne DNormal w y mu dev) - 2)
    vW
    vY
    vMu
    vDev
  Bernoulli -> VS.sum $ VS.zipWith4 (aicOne DBernoulli) vW vY vMu vDev
  Binomial vN ->
    VS.sum $ VS.zipWith5 (\n -> aicOne (DBinomial n)) vN vW vY vMu vDev
  Poisson -> VS.sum $ VS.zipWith4 (aicOne DPoisson) vW vY vMu vDev
  Gamma   -> 2 + VS.sum $ VS.zipWith4
    (\w y mu dev -> (aicOne DGamma w y mu dev) - 2)
    vW
    vY
    vMu
    vDev
{-    
  Normal -> case vDevM of
    Nothing -> Nothing
    (Just vDev ->
      let n       = realToFrac $ LA.size vMu
          sumLnWt = VS.sum $ VS.map log vW
      in  Just (n * (log (2 * pi * dev / n) + 1) + 2 - sumLnWt)
  Bernoulli ->
    let vSuccess = VS.map round vY
--        logProbBernoulli y mu  = y* log mu + (1-y)*log(1-mu) 
        logProbBernoulli succ prob = S.logProbability (S.binomial 1 prob) succ
        wgtdLPB wgt succ prob = wgt * logProbBernoulli succ prob
        ans = VS.sum $ VS.zipWith3 wgtdLPB vW vSuccess vMu
    in  Just $ -2 * ans
  Binomial counts
    -> let
         vNSuccess :: VS.Vector Int =
           VS.zipWith (\n x -> round $ realToFrac n * x) counts vY
         logProbBinom nSucc nTrials prob =
           S.logProbability (S.binomial nTrials prob) nSucc
         wgtdLPB wgt cnt succ prob =
           (wgt / realToFrac cnt) * logProbBinom succ cnt prob
         ans = VS.sum $ VS.zipWith4 wgtdLPB vW counts vNSuccess vMu
       in
         Just $ -2 * ans
  Poisson ->
    let vCount :: VS.Vector Int = VS.map round vY
        logProbPoisson count prob = S.logProbability (S.poisson prob) count
        wgtdLPP wgt count prob = wgt * logProbPoisson count prob
        ans = VS.sum $ VS.zipWith3 wgtdLPP vW vCount vMu
    in  Just $ -2 * ans
  Gamma fs -> case fs of
    Unfitted -> Nothing
    FittedScale dev ->
      let
        sumWt   = VS.sum vW
        disp    = dev / sumWt
        invDisp = 1 / disp
        logDensityGamma y mu =
          S.logDensity (S.gammaDistr invDisp (mu * disp)) y
        wgtdLDG wgt y mu = wgt * logDensityGamma y mu
        ans = VS.sum $ VS.zipWith3 wgtdLDG vW vY vMu
      in
        Just (2 - 2 * ans)
-}

-- this needs a better name.  Which means I need to understand it.
-- scaled variance?
variance :: ObservationDistribution -> LA.Vector Double -> LA.Vector Double
variance od =
  let eps = 1e-12
      f x = case od of
        (Normal   _) -> 1 -- this is wrong so maybe this is a scaled variance for distributions with an overall scale??
        (Binomial _) -> x * (1 - x)
        Bernoulli    -> x * (1 - x)  -- ??
        (Poisson _)  -> x
        (Gamma   _)  -> x * x
      g x = if x > eps then f x else eps
  in  VS.map g


-- Numeric helpers 
{-
yLny :: Double -> Double
yLny x =
  let eps = 1e-12 -- do something type/machine appropriate here
  in case x < eps of
    True -> eps
-}
