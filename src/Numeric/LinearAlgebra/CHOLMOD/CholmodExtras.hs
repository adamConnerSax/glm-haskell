{-# LANGUAGE CPP                      #-}
{-# LANGUAGE EmptyDataDecls           #-}
{-# LANGUAGE ForeignFunctionInterface #-} 
{-# LANGUAGE OverloadedStrings #-}

module Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
  (
    spMatrixAnalyze
  , spMatrixCholesky
  , unsafeSpMatrixCholesky
  , factorPermutation
  , factorPermutationSM
  , choleskyFactorSM
  , spMatrixToTriplet
  , tripletToSpMatrix
  , readv
  , writev
  -- * low-level
  , printCommon
  , printFactor
  , printTriplet
  , printSparse
  , setFinalLL
  -- * re-exports
  , Factor
  , Common
  , ForeignPtr
  , allocCommon
  , startC
  )
where

import Foreign hiding (free)
import Foreign.C.Types
import Foreign.C.String
import Foreign.Storable (peek)
import System.IO.Unsafe (unsafePerformIO)

import qualified Data.Text as T

import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse          as SLA

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow as CH
import           Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow  (Common, Factor)
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFace as CH
import           Numeric.LinearAlgebra.CHOLMOD.CholmodXFace  (allocCommon, startC)

import qualified Data.Vector.Storable.Mutable as V

-- | make an SpMatrix into a CHOLMOD triplet
-- assumes the SpMatrix is square and symmetric and only uses lower triangle
spMatrixToTriplet :: RealFrac a => ForeignPtr CH.Common -> SLA.SpMatrix a -> IO (CH.Matrix CH.Triplet)
spMatrixToTriplet fpc smX = do
  let (nrows, ncols) = SLA.dimSM smX
      nzMax = SLA.nzSM smX
  triplet <- CH.allocTriplet
             (fromIntegral nrows)
             (fromIntegral ncols)
             (fromIntegral nzMax)
             CH.stUnsymmetric
             CH.xtReal
             fpc
  rowIs <- CH.tripletGetRowIndices triplet
  colIs <- CH.tripletGetColIndices triplet
  vals  <- CH.tripletGetX triplet
  let smTriplets = SLA.toListSM smX
  writev rowIs $ fmap (\(rI,_,_) -> fromIntegral rI) smTriplets
  writev colIs $ fmap (\(_,cI,_) -> fromIntegral cI) smTriplets
  writev vals $ fmap (\(_,_,x) -> realToFrac x) smTriplets
  CH.tripletSetNNZ triplet (fromIntegral nzMax)
  return triplet

-- | Compute fill-reducing permutation, etc. for a symmetric positive-definite matrix
-- This only requires that the lower triangle be filled in
spMatrixAnalyze :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                -> SLA.SpMatrix Double -- ^ matrix to analyze
                -> IO (ForeignPtr CH.Factor, SLA.SpMatrix Double) -- ^ analysis and fill-reducing permutation
spMatrixAnalyze fpc smX = do
--  withForeignPtr fpc $ \cp -> setFinalLL 1 cp
  triplet <- spMatrixToTriplet fpc smX
  printTriplet (CH.fPtr triplet) "spMatrixAnalyze" fpc
  sparse <- CH.tripletToSparse triplet fpc
  f <- CH.analyze sparse fpc
  printFactor f "spMatrixAnalyze" fpc
  permSM <- factorPermutationSM f
  return (f, permSM)
  
-- | compute the lower-triangular Cholesky factor using the given analysis stored in Factor
-- the matrix given here must have the same pattern of non-zeroes as the one used for the
-- analysis
spMatrixCholesky :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                 -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                 -> SLA.SpMatrix Double -- ^ matrix to decompose
                 -> IO (SLA.SpMatrix Double) -- ^ lower-triangular Cholesky factor
spMatrixCholesky fpc fpf smX = do  
  triplet <- spMatrixToTriplet fpc smX
  printTriplet (CH.fPtr triplet) "spMatrixCholesky" fpc
  sparse <- CH.tripletToSparse triplet fpc
  printSparse (CH.fPtr sparse) "spMatrixCholesky" fpc
  CH.factorize sparse fpf fpc
  printFactor fpf "spMatrixCholesky (after factorize):" fpc
  choleskyFactorSM fpf fpc

-- |  compute the lower-triangular Cholesky factor using the given analysis stored in Factor
-- the matrix given here must have the same pattern of non-zeroes as the one used for the
-- analysis
-- NB: THis version calls unsafePerformIO which, I think, is okay because the function is referentially transparent.
-- At least in a single threaded environment. Oy.
unsafeSpMatrixCholesky ::  ForeignPtr CH.Common -- ^ CHOLMOD environment
                       -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                       -> SLA.SpMatrix Double -- ^ matrix to decompose
                       -> SLA.SpMatrix Double -- ^ lower-triangular Cholesky factor
unsafeSpMatrixCholesky fpc fpf smX = unsafePerformIO $ spMatrixCholesky fpc fpf smX
{-# NOINLINE unsafeSpMatrixCholesky #-} 
  
factorPermutation :: ForeignPtr CH.Factor -> IO (V.MVector s CInt)
factorPermutation ffp = withForeignPtr ffp $ \fp -> do
  n <- factorN_L fp :: IO CSize
  pi <- factorPermutationL fp :: IO (Ptr CInt)
  pifp <- newForeignPtr_ pi :: IO (ForeignPtr CInt)
  return $ V.unsafeFromForeignPtr0 pifp (fromIntegral n)

-- | retrieve the factor permutation as an SpMatrix
factorPermutationSM :: Num a => ForeignPtr CH.Factor -> IO (SLA.SpMatrix a)
factorPermutationSM ffp = do
  permMV <- factorPermutation ffp
  permL <- readv permMV
  return $ SLA.permutationSM (V.length permMV) (fmap fromIntegral permL)

-- | retrieve the Cholesky Factor as an SLA.SpMatrix
-- | NB: This undoes the factorization in Factor
choleskyFactorSM :: ForeignPtr CH.Factor -> ForeignPtr CH.Common -> IO (SLA.SpMatrix Double)
choleskyFactorSM ffp cfp = withForeignPtr ffp $ \fp -> do
  withForeignPtr cfp $ \cp -> do
    printFactor ffp "Before" cfp
    changeFactorL (fromIntegral $ CH.unXType CH.xtReal) 1 0 0 0 fp cp
    printFactor ffp "After" cfp
    pSparse <- factorToSparseL fp cp :: IO (Ptr CH.Sparse)
    pTriplet <- sparseToTripletL pSparse cp :: IO (Ptr CH.Triplet)
    tripletToSpMatrix pTriplet
{-    
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
        sm = SLA.fromListSM (fromIntegral nRows, fromIntegral nCols) triplets
    SLA.prd sm    
    return sm
-}

tripletToSpMatrix :: Ptr CH.Triplet -> IO (SLA.SpMatrix Double)
tripletToSpMatrix pTriplet = do
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
      sm = SLA.fromListSM (fromIntegral nRows, fromIntegral nCols) triplets
  SLA.prd sm    
  return sm
  

printCommon :: T.Text -> ForeignPtr CH.Common -> IO ()
printCommon t fpc = withForeignPtr fpc $ \fp -> do
  tCS <- newCString (T.unpack t)
  printCommonL tCS fp

printFactor :: ForeignPtr CH.Factor -> T.Text -> ForeignPtr CH.Common -> IO ()
printFactor fpf t fpc = withForeignPtr fpf $ \pf -> do
  withForeignPtr fpc $ \pc -> do
    tSC <- newCString (T.unpack t)
    printFactorL pf tSC pc

printTriplet :: ForeignPtr CH.Triplet -> T.Text -> ForeignPtr CH.Common -> IO ()
printTriplet fpt t fpc = withForeignPtr fpt $ \pt -> do
  withForeignPtr fpc $ \pc -> do
    tSC <- newCString (T.unpack t)
    printTripletL pt tSC pc

printSparse :: ForeignPtr CH.Sparse -> T.Text -> ForeignPtr CH.Common -> IO ()
printSparse fps t fpc = withForeignPtr fps $ \ps -> do
  withForeignPtr fpc $ \pc -> do
    tSC <- newCString (T.unpack t)
    printSparseL ps tSC pc

setFinalLL :: Int -> ForeignPtr CH.Common -> IO ()
setFinalLL n fpc = withForeignPtr fpc $ \pc -> setFinalLLL n pc

readv :: V.Storable a => V.IOVector a -> IO [a]
readv v = sequence [ V.read v i | i <- [0 .. (V.length v) - 1] ]

writev :: (Storable a) => V.IOVector a -> [a] -> IO ()
writev v xs =
  sequence_ [ V.write v i x | (i, x) <- zip [0 .. (Prelude.length xs - 1)] xs ]

-- | the N of the N x N factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_n"
  factorN_L :: Ptr CH.Factor -> IO CSize

-- | the array of indices representing the permutation matrix in the factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_permutation"
  factorPermutationL :: Ptr CH.Factor -> IO (Ptr CInt)

foreign import ccall unsafe "cholmod.h cholmod_factor_to_sparse"
  factorToSparseL :: Ptr CH.Factor -> Ptr CH.Common -> IO (Ptr CH.Sparse) 

foreign import ccall unsafe "cholmod.h cholmod_sparse_to_triplet"
  sparseToTripletL :: Ptr CH.Sparse -> Ptr CH.Common -> IO (Ptr CH.Triplet)

foreign import ccall unsafe "cholmod.h cholmod_change_factor"
  changeFactorL :: Int -> Int -> Int -> Int -> Int -> Ptr CH.Factor -> Ptr CH.Common -> IO CInt 

foreign import ccall unsafe "cholmod_extras.h cholmod_set_final_ll"
  setFinalLLL :: Int -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_common"
  printCommonL :: CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_factor"
  printFactorL :: Ptr CH.Factor -> CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_triplet"
  printTripletL :: Ptr CH.Triplet -> CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_sparse"
  printSparseL :: Ptr CH.Sparse -> CString -> Ptr CH.Common -> IO ()



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


