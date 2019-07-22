{-# LANGUAGE CPP                      #-}
{-# LANGUAGE ForeignFunctionInterface #-} 
{-# LANGUAGE EmptyDataDecls           #-}

module Numeric.LinearAlgebra.CHOLMOD.CholmodExtras where

import Foreign hiding (free)
import Foreign.C.Types
import Foreign.Storable (peek)

import qualified Data.Sparse.SpMatrix          as SLA

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow as CH
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFace as CH

import qualified Data.Vector.Storable.Mutable as V

permutation :: ForeignPtr CH.Factor -> IO (V.MVector s CInt)
permutation ffp = withForeignPtr ffp $ \fp -> do
  n <- factorN_L fp :: IO CSize
  putStrLn $ "n=" ++ show n
  pi <- factorPermutationL fp :: IO (Ptr CInt)
  pifp <- newForeignPtr_ pi :: IO (ForeignPtr CInt)
  return $ V.unsafeFromForeignPtr0 pifp (fromIntegral n)
  
choleskyFactorSM :: ForeignPtr CH.Factor -> ForeignPtr CH.Common -> IO (SLA.SpMatrix Double)
choleskyFactorSM ffp cfp = withForeignPtr ffp $ \fp -> do
  withForeignPtr cfp $ \cp -> do
    pSparse <- factorToSparseL fp cp :: IO (Ptr CH.Sparse)
    pTriplet <- sparseToTriplet pSparse cp :: IO (Ptr CH.Triplet)
    nRows <- CH.triplet_get_nrow pTriplet
    nCols <- CH.triplet_get_ncol pTriplet
    nZmax <- CH.triplet_get_nzmax pTriplet
    rowIndicesP <- CH.triplet_get_row_indices pTriplet :: IO (Ptr CInt)
    rowIndicesFP <- newForeignPtr_ rowIndicesP    
    colIndicesP <- CH.triplet_get_col_indices pTriplet :: IO (Ptr CInt)
    colIndicesFP <- newForeignPtr_ colIndicesP
    valsP <- CH.triplet_get_x pTriplet  :: IO (Ptr CDouble)
    valsFP <- newForeignPtr_ valsP
    rowIndices <- readv (V.unsafeFromForeignPtr0 rowIndicesFP (fromIntegral nZmax))
    colIndices <- readv (V.unsafeFromForeignPtr0 colIndicesFP (fromIntegral nZmax))
    vals <- readv (V.unsafeFromForeignPtr0 valsFP (fromIntegral nZmax))
    let triplets = zip3 (fmap fromIntegral rowIndices) (fmap fromIntegral colIndices) (fmap realToFrac vals)
    return $ SLA.fromListSM (fromIntegral nRows, fromIntegral nCols) triplets
  
readv :: V.Storable a => V.IOVector a -> IO [a]
readv v = sequence [ V.read v i | i <- [0 .. (V.length v) - 1] ]


-- | the N of the N x N factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_n"
  factorN_L :: Ptr CH.Factor -> IO CSize

-- | the array of indices representing the permutation matrix in the factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_permutation"
  factorPermutationL :: Ptr CH.Factor -> IO (Ptr CInt)

foreign import ccall unsafe "cholmod.h cholmod_factor_to_sparse"
  factorToSparseL :: Ptr CH.Factor -> Ptr CH.Common -> IO (Ptr CH.Sparse) 

foreign import ccall unsafe "cholmod.h cholmod_sparse_to_triplet"
  sparseToTriplet :: Ptr CH.Sparse -> Ptr CH.Common -> IO (Ptr CH.Triplet)




{-
-- | the factor L
-- NB: for LDLt decompositions, the diagonal of L is all 1s so CHOLMOD puts
-- D on the diagonal here instead.  I think.
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_get_nrow"
  factor_sparse :: Ptr CH.Common -> Ptr CH.Factor -> IO (CH.Matrix CH.Sparse)
-}

-- TODO
-- SpMatrix to cholmod_sparse or cholmod_triplet
-- cholmod_sparse to SpMatrix
-- cholmod permutation ([Int] ?) to permutationSM
-- permutationSM * (dense) vector
-- driver for LLt decomposition, holding factor constant.
-- possibly some way of simplifying all the copying to
-- do this in an iterative solver.  


