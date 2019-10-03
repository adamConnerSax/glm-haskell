{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.GeneralizedMixedModel
  ( module Numeric.GLM.GeneralizedMixedModel
  , module Numeric.GLM.MixedModel
  )
where
import qualified Numeric.GLM.Types             as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.SparseDenseConversions
                                               as SD
import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Numeric.LinearAlgebra.Class   as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import           Numeric.LinearAlgebra.Sparse   ( (##)
                                                , (#^#)
                                                , (<#)
                                                , (#>)
                                                , (-=-)
                                                , (-||-)
                                                )
import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL



--import           Data.Either                    ( partitionEithers )
import qualified Data.List                     as L
--import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS




data LinkFunctionType = IdentityLink | LogisticLink | PoissonLink deriving (Show, Eq)
data LinkFunction a = LinkFunction { link :: a -> a -- map from observation to linear predictor
                                   , invLink :: a -> a -- map from linear predictor to obervations
                                   , derivInv :: a -> a -- useful in PIRLS algo
                                   }



data GeneralizedMixedModel b g = GeneralizedMixedModel { mixedModel :: MixedModel b g, weights :: WMatrix, linkType :: LinkFunctionType }

linkFunction :: Floating a => LinkFunctionType -> LinkFunction a
linkFunction IdentityLink = LinkFunction id id (const 1)

linkFunction LogisticLink = LinkFunction (\x -> log $ x / (1 - x))
                                         (\x -> exp x / (1 + exp x))
                                         (\x -> exp x / (1 + exp x) ^^ 2)

linkFunction PoissonLink = LinkFunction log exp exp

--type WMatrix = LA.Vector Double -- constant weights
--type UMatrix = SLA.SpMatrix Double
--type VMatrix = LA.Matrix Double


spUV
  :: GeneralizedMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVector
  -> SLA.SpVector Double
  -> SLA.SpVector Double
  -> (UMatrix, VMatrix)
spUV (GeneralizedMixedModel (MixedModel (RegressionModel _ mX vY) _) vW lft) (RandomEffectCalculated smZ mkLambda) vTh svBeta svu
  = let
      (LinkFunction l inv dinv) = linkFunction lft
      spLambda                  = mkLambda vTh
      vEta = SD.toDenseVector (smZ #> svu) + mX LA.#> (SD.toDenseVector svBeta)
      vdMudEta                  = VS.map dinv vEta
      mWdMudEta                 = LA.diag $ VS.zipWith (*) vW vdMudEta
      smU = SD.toSparseMatrix mWdMudEta SLA.#~# smZ SLA.#~# spLambda
      mV                        = mWdMudEta LA.<> mX
    in
      (smU, mV)

{-
cholmodDeltas
  :: CholmodFactor
  -> GeneralizedMixedModel b g
  -> MatrixU
  -> MatrixV
  -> Lambda
  -> IO (du,dBeta)
cholmodDeltas cholmodFactor gmm vU vV smL =
  let
    (cholmodC, cholmodF, _) = cholmodFactor
    
-}
