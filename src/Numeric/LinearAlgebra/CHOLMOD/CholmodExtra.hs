{-# LANGUAGE CPP                      #-}
{-# LANGUAGE ForeignFunctionInterface #-} 
{-# LANGUAGE EmptyDataDecls           #-}

module Numeric.LinearAlgebra.CHOLMOD.CholmodExtra where

import Foreign hiding (free)
import Foreign.C.Types

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow as CH
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFace as CH

-- | the array of indices representing the permutation matrix in the factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_permutation"
  permutation :: Ptr CH.Factor -> IO (Ptr CInt) 

-- | the factor L
-- NB: for LDLt decompositions, the diagonal of L is all 1s so CHOLMOD puts
-- D on the diagonal here instead.  I think.
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_get_nrow"
  factor_sparse :: Ptr CH.Common -> Ptr CH.Factor -> IO (CH.Matrix CH.Sparse)

-- TODO
-- SpMatrix to cholmod_sparse or cholmod_triplet
-- cholmod_sparse to SpMatrix
-- cholmod permutation ([Int] ?) to permutationSM
-- permutationSM * (dense) vector
-- driver for LLt decomposition, holding factor constant.
-- possibly some way of simplifying all the copying to
-- do this in an iterative solver.  


