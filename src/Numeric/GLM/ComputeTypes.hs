{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Numeric.GLM.ComputeTypes
  ()
where

import qualified Data.List                     as L
import qualified Data.Map                      as M
import qualified Data.Text                     as T
import qualified Data.Vector.Storable          as VB
import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA

