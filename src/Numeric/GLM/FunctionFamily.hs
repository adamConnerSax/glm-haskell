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
import qualified Data.Vector.Storable          as LA
import qualified Numeric.MathFunctions.Constants
                                               as MC


data ObservationDistribution = Normal
                             | Binomial (LA.Vector Int) -- vector of total counts
                             | Poisson
                             | Gamma deriving (Show, Eq)



data LinkFunctionType = IdentityLink | LogisticLink | ExponentialLink deriving (Show, Eq)

data LinkFunction = LinkFunction { link :: Double -> Double -- map from observation to linear predictor
                                 , invLink :: Double -> Double -- map from linear predictor to obervations
                                 , derivInv :: Double -> Double -- useful in PIRLS algo
                                 }




canonicalLink :: ObservationDistribution -> LinkFunctionType
canonicalLink Normal       = IdentityLink
canonicalLink (Binomial _) = LogisticLink
canonicalLink Poisson      = ExponentialLink

linkFunction :: LinkFunctionType -> LinkFunction
linkFunction IdentityLink = LinkFunction id id (const 1)

linkFunction LogisticLink = LinkFunction
  (\x -> log (x / (1 - x)))
  (\x -> let y = exp x in y / (1 + y))
  (\x -> let y = exp x in y / ((1 + y) ** 2))

linkFunction ExponentialLink = LinkFunction log exp exp


data UseLink = UseCanonical | UseOther LinkFunctionType deriving (Show, Eq)

-- notation here:
-- y are the observations
-- mu is the conditional mean of the linear predictor after mapping via the link function

familyWeights :: ObservationDistribution -> LA.Vector Double -> LA.Vector Double
familyWeights (Binomial vN) vW =
  let f w n = w * realToFrac n in LA.zipWith f vW vN
familyWeights _ vW = vW

deviance
  :: ObservationDistribution
  -> UseLink
  -> LA.Vector Double
  -> LA.Vector Double
  -> LA.Vector Double
  -> Double
deviance od ul vW vY vEta
  = let
      eps     = 1e-12
      weights = familyWeights od vW
      iLink   = invLink . linkFunction $ case ul of
        UseCanonical -> canonicalLink od
        UseOther lft -> lft
      vMu = LA.map iLink vEta
      f y mu = case od of
        Normal       -> (y - mu) ** 2
        (Binomial _) -> if y < eps
          then -2 * log (1 - mu) -- y = 0
          else if (1 - y) < eps
            then -2 * log mu  -- y = 1
            else 2 * (y * log (y / mu) + (1 - y) * log ((1 - y) / (1 - mu)))
--            else -2 * (y * log mu + (1 - y) * log (1 - mu)) -- 
        Poisson ->
          if y > eps then 2 * (y * log (y / mu) - (y - mu)) else 2 * mu
        Gamma -> 2 * ((y - mu) / mu - log (y / mu))
      g h w y mu = w * h y mu
    in
      LA.sumElements $ LA.zipWith3 (g f) weights vY vMu
{-
  go Normal =
    let f y mu = (y - mu) ^^ 2 in LA.sumElements $ LA.zipWith f vY vMu'
  go (Binomial vN) =
    let f n y mu =
          let m = realToFrac n
          in  m * (y * log (y / mu) - (1 - y) * log ((1 - y) / (1 - mu)))
    in  2 * (LA.sumElements $ LA.zipWith3 f vN vY vMu')
  go Poisson =
    let f y mu = y * log (y / mu) - (y - mu)
    in  2 * (LA.sumElements $ LA.zipWith f vY vMu')
  go Gamma =
    let f y mu = (y - mu) / mu - log (y / mu)
    in  2 * (LA.sumElements $ LA.zipWith f vY vMu')
-}

