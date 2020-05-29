{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.Confidence
  ( module Numeric.GLM.Bootstrap
  , module Numeric.GLM.Confidence
  )
where

import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Vector.Storable          as VS
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.ProblemTypes      as GLM

import qualified Numeric.GLM.Predict           as GLM
import qualified Numeric.GLM.Bootstrap         as GLM
import           Numeric.GLM.Bootstrap          ( BootstrapCIType(..) )

import qualified Statistics.Types              as S

import qualified Knit.Effect.Logger            as P



-- TODO: I should add a likelihood test.  Once I figure out how to do it.
data ConfidenceType = NaiveCondVarCI (LA.Matrix Double) GLM.ConditionalCovarianceMatrix
                    | BootstrapCI GLM.BootstrapCIType [(GLM.BetaVec, VS.Vector Double)]

predictWithCI
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe (GLM.GroupKey g))
  -> GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> GLM.BetaVec
  -> VS.Vector Double --vb
  -> S.CL Double
  -> ConfidenceType
  -> P.Sem r (Double, (Double, Double))
predictWithCI mm getPredictorM getLabelM rc ebg vBeta vb cl (NaiveCondVarCI mBetaCov smCondVar)
  = GLM.predictWithCondVarCI mm
                             getPredictorM
                             getLabelM
                             ebg
                             rc
                             vBeta
                             vb
                             cl
                             mBetaCov
                             smCondVar

predictWithCI mm getPredictorM getLabelM rc ebg vBeta vb cl (BootstrapCI bootCIType bootBetabs)
  = GLM.bootstrappedConfidence mm
                               getPredictorM
                               getLabelM
                               rc
                               ebg
                               vBeta
                               vb
                               bootBetabs
                               bootCIType
                               cl

