{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.GLM.LinearMixedModel
  ( module Numeric.GLM.LinearMixedModel
  , module Numeric.GLM.MixedModel
  )
where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.Types             as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.SparseDenseConversions
                                               as SD

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Map                      as M
import           Data.Maybe                     ( catMaybes )
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA

import qualified Numeric.LinearAlgebra.Class   as SLA
{-
import qualified Data.Sparse.SpVector          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import           Numeric.LinearAlgebra.Sparse   ( (##)
                                                , (#^#)
                                                , (<#)
                                                , (#>)
                                                , (-=-)
                                                , (-||-)
                                                )
-}
import qualified Numeric.LinearAlgebra         as LA





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






--data SolutionComponents = SolutionComponents { mRX :: LA.Matrix Double, vBeta :: LA.Vector Double, vTheta :: LA.Vector Double, vb :: LA.Vector Double }

