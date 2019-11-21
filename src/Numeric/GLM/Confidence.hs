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
import qualified Numeric.GLM.Predict           as GLM
import qualified Numeric.GLM.Bootstrap         as GLM
import           Numeric.GLM.Bootstrap          ( BootstrapCIType(..) )
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



-- TODO: I should add a likelihood test.  Once I figure out how to do it.
data ConfidenceType = NaiveCondVarCI (LA.Matrix Double) GLM.ConditionalCovarianceMatrix
                    | BootstrapCI GLM.BootstrapCIType [(GLM.BetaU, VS.Vector Double)]

predictWithCI
  :: (Enum b, Bounded b, Ord b, Ord g, Show g, Show b, GLM.Effects r)
  => GLM.MixedModel b g
  -> (b -> Maybe Double)
  -> (g -> Maybe T.Text)
  -> GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> GLM.BetaU
  -> VS.Vector Double --vb
  -> S.CL Double
  -> ConfidenceType
  -> P.Sem r (Double, (Double, Double))
predictWithCI mm getPredictorM getLabelM rc ebg betaU vb cl (NaiveCondVarCI mBetaCov smCondVar)
  = GLM.predictWithCondVarCI mm
                             getPredictorM
                             getLabelM
                             ebg
                             rc
                             betaU
                             vb
                             cl
                             mBetaCov
                             smCondVar

predictWithCI mm getPredictorM getLabelM rc ebg betaU vb cl (BootstrapCI bootCIType bootBetaUbs)
  = GLM.bootstrappedConfidence mm
                               getPredictorM
                               getLabelM
                               rc
                               ebg
                               betaU
                               vb
                               bootBetaUbs
                               bootCIType
                               cl

