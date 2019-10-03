{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.LinkFunction where


data LinkFunctionType = IdentityLink | LogisticLink | PoissonLink deriving (Show, Eq)
data LinkFunction a = LinkFunction { link :: a -> a -- map from observation to linear predictor
                                   , invLink :: a -> a -- map from linear predictor to obervations
                                   , derivInv :: a -> a -- useful in PIRLS algo
                                   }

linkFunction :: Floating a => LinkFunctionType -> LinkFunction a
linkFunction IdentityLink = LinkFunction id id (const 1)

linkFunction LogisticLink = LinkFunction (\x -> log $ x / (1 - x))
                                         (\x -> exp x / (1 + exp x))
                                         (\x -> exp x / (1 + exp x) ^^ 2)

linkFunction PoissonLink = LinkFunction log exp exp






