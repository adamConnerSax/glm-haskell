{-# LANGUAGE CPP                      #-}
{-# LANGUAGE EmptyDataDecls           #-}
{-# LANGUAGE ForeignFunctionInterface #-} 
{-# LANGUAGE OverloadedStrings #-}

module Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
  (
    MatrixSymmetry (..)
  , SolveSystem (..)
  , solveSparse
  , spMatrixAnalyze
  , spMatrixAnalyzeWP
  , spMatrixCholesky
  , unsafeSpMatrixCholesky
  , spMatrixFactorize
  , CholmodFactorizationStyle (..)
  , spMatrixFactorizeP
  , unsafeSpMatrixFactorize
  , factorPermutation
  , factorPermutationSM
  , unsafeFactorPermutationSM
  , choleskyFactorSM
  , unsafeCholeskyFactorSM
  , spMatrixToTriplet
  , tripletToSpMatrix
  , readv
  , writev
  , setFinalLL
  -- * low-level
  , printCommon
  , printFactor
  , printTriplet
  , printSparse
  , copyFactor
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
--import Foreign.Storable (peek)
--import Foreign.Marshal.Alloc (alloca)
import System.IO.Unsafe (unsafePerformIO)


import           Control.Monad (when)
import qualified Data.Text as T

import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse          as SLA

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow as CH
import           Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow  (Common, Factor)
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodXFace as CH
import           Numeric.LinearAlgebra.CHOLMOD.CholmodXFace  (allocCommon, startC)

import qualified Data.Vector.Storable.Mutable as MVS
import qualified Data.Vector.Storable as VS

data MatrixSymmetry = UnSymmetric | SquareSymmetricUpper | SquareSymmetricLower

data SolveSystem =
  CHOLMOD_A      -- Ax    = b
  | CHOLMOD_LDLt -- LDL'x = b 
  | CHOLMOD_LD   -- LDx   = b
  | CHOLMOD_DLt  -- DL'x  = b
  | CHOLMOD_L    -- Lx    = b
  | CHOLMOD_Lt   -- L'x   = b
  | CHOLMOD_D    -- Dx    = b
  | CHOLMOD_P    -- x     = Pb
  | CHOLMOD_Pt   -- x     = P'b
    deriving (Enum, Show, Eq)

cholmodSystem :: SolveSystem -> CInt
cholmodSystem CHOLMOD_A = 0
cholmodSystem CHOLMOD_LDLt = 1
cholmodSystem CHOLMOD_LD = 2
cholmodSystem CHOLMOD_DLt = 3
cholmodSystem CHOLMOD_L = 4
cholmodSystem CHOLMOD_Lt = 5
cholmodSystem CHOLMOD_D = 6
cholmodSystem CHOLMOD_P = 7
cholmodSystem CHOLMOD_Pt = 8

hcholmodSType :: MatrixSymmetry -> CH.SType
hcholmodSType UnSymmetric = CH.stUnsymmetric
hcholmodSType SquareSymmetricUpper = CH.stSquareSymmetricUpper
hcholmodSType SquareSymmetricLower = CH.stSquareSymmetricLower

debug = False
debugSolve = False

solveSparse :: ForeignPtr CH.Common -- ^ CHOLMOD environment
            -> ForeignPtr CH.Factor -- ^ contains P and LL' (or LDL')  
            -> SolveSystem -- ^ which system to solve
            -> SLA.SpMatrix Double -- ^ RHS
            -> IO (SLA.SpMatrix Double) -- ^ solutions
solveSparse fpc fpf ss smB = do
  mTripletB <- spMatrixToTriplet fpc UnSymmetric smB
  when debugSolve $ printTriplet (CH.fPtr mTripletB) "solveSparse input triplet" fpc
  mSparseB <- CH.tripletToSparse mTripletB fpc
  when debugSolve $ printSparse (CH.fPtr mSparseB) "solveSparse input sparse" fpc
  withForeignPtr fpc $ \pc -> do
    withForeignPtr fpf $ \pf -> do
      withForeignPtr (CH.fPtr mSparseB) $ \pSparseB -> do
        pSparseX <- sparseSolveL (cholmodSystem ss) pf pSparseB pc
        when debugSolve $ printSparse' pSparseX "solveSparse output sparse" fpc
        pTripletX <- sparseToTripletL pSparseX pc
        when debugSolve $ printTriplet' pTripletX "solveSparse output triplet" fpc
        smX <- tripletToSpMatrix pTripletX
--        when debugSolve $ SLA.prd smX
        CH.sparse_free pSparseB pc
        withForeignPtr (CH.fPtr mTripletB) $ \pt -> CH.triplet_free pt pc       
        CH.triplet_free pTripletX pc
        CH.sparse_free pSparseX pc
        return smX
        
getDims :: MatrixSymmetry ->  SLA.SpMatrix a -> (CSize, CSize, CSize)
getDims ms smX =
  let (nrows, ncols) = SLA.dimSM smX
      nzMax = SLA.nzSM smX
      nnz = case ms of
        UnSymmetric -> nzMax
        _ -> min nzMax (nrows * (nrows + 1) `div` 2)
  in (fromIntegral nrows, fromIntegral ncols, fromIntegral nnz)
    
-- | make an SpMatrix into a CHOLMOD triplet
-- this thing should, but doesn't, free itself via ForeignPtr.
spMatrixToTriplet :: RealFrac a => ForeignPtr CH.Common -> MatrixSymmetry -> SLA.SpMatrix a -> IO (CH.Matrix CH.Triplet)
spMatrixToTriplet fpc ms smX = do
  let (nrows, ncols, nnz) = getDims ms smX
  triplet <- CH.allocTriplet
             nrows
             ncols
             nnz
             (hcholmodSType ms)
             CH.xtReal
             fpc
  rowIs <- CH.tripletGetRowIndices triplet
  colIs <- CH.tripletGetColIndices triplet
  vals  <- CH.tripletGetX triplet
  let symmetryFilter (r,c,_) = case ms of
        UnSymmetric -> True
        SquareSymmetricUpper -> (c >= r)
        SquareSymmetricLower -> (c <= r)
      smTriplets = filter symmetryFilter $ SLA.toListSM smX
  writev rowIs $ fmap (\(rI,_,_) -> fromIntegral rI) smTriplets
  writev colIs $ fmap (\(_,cI,_) -> fromIntegral cI) smTriplets
  writev vals $ fmap (\(_,_,x) -> realToFrac x) smTriplets
  CH.tripletSetNNZ triplet (fromIntegral $ length smTriplets)
  return triplet

-- | Compute fill-reducing permutation, etc. for a symmetric positive-definite matrix
-- This only requires that the lower triangle be filled in
spMatrixAnalyzeWP :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                -> MatrixSymmetry
                -> SLA.SpMatrix Double -- ^ matrix to analyze
                -> IO (ForeignPtr CH.Factor, SLA.SpMatrix Double) -- ^ analysis and fill-reducing permutation
spMatrixAnalyzeWP fpc ms smX = do
  triplet <- spMatrixToTriplet fpc ms smX
  when debug $ printTriplet (CH.fPtr triplet) "spMatrixAnalyze" fpc
  sparse <- CH.tripletToSparse triplet fpc
  f <- CH.analyze sparse fpc
  when debug $ printFactor f "spMatrixAnalyze" fpc
  permSM <- factorPermutationSM f
  withForeignPtr fpc $ \pc -> do
    withForeignPtr (CH.fPtr triplet) $ \pt -> CH.triplet_free pt pc    
  return (f, permSM)

-- | Compute fill-reducing permutation, etc. for a symmetric positive-definite matrix
-- This only requires that the lower triangle be filled in
spMatrixAnalyze :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                -> MatrixSymmetry
                -> SLA.SpMatrix Double -- ^ matrix to analyze
                -> IO (ForeignPtr CH.Factor)
spMatrixAnalyze fpc ms smX = do
  triplet <- spMatrixToTriplet fpc ms smX
  when debug $ printTriplet (CH.fPtr triplet) "spMatrixAnalyze" fpc
  sparse <- CH.tripletToSparse triplet fpc
  f <- CH.analyze sparse fpc
  when debug $ printFactor f "spMatrixAnalyze" fpc
  withForeignPtr fpc $ \pc -> do
    withForeignPtr (CH.fPtr triplet) $ \pt -> CH.triplet_free pt pc    
  return f

  
-- | compute the lower-triangular Cholesky factor using the given analysis stored in Factor
-- the matrix given here must have the same pattern of non-zeroes as the one used for the
-- analysis
spMatrixCholesky :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                 -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                 -> MatrixSymmetry
                 -> SLA.SpMatrix Double -- ^ matrix to decompose
                 -> IO (SLA.SpMatrix Double) -- ^ lower-triangular Cholesky factor
spMatrixCholesky fpc fpf ms smX = do
  spMatrixFactorize fpc fpf ms smX
  choleskyFactorSM fpf fpc

spMatrixFactorize :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                  -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                  -> MatrixSymmetry
                  -> SLA.SpMatrix Double -- ^ matrix to decompose
                  -> IO () -- ^ The factor in ForeignPtr Factor will be updated :(
spMatrixFactorize fpc fpf ms smX = do
  triplet <- spMatrixToTriplet fpc ms smX
  when debug $ printTriplet (CH.fPtr triplet) "spMatrixCholesky" fpc
  sparse <- CH.tripletToSparse triplet fpc
  when debug $ printSparse (CH.fPtr sparse) "spMatrixCholesky" fpc
  CH.factorize sparse fpf fpc
  withForeignPtr fpc $ \pc -> do
    withForeignPtr (CH.fPtr triplet) $ \pt -> CH.triplet_free pt pc
    withForeignPtr (CH.fPtr sparse) $ \ps -> CH.sparse_free ps pc
  return ()

data CholmodFactorizationStyle = FactorizeA | FactorizeAtAPlusBetaI Double

spMatrixFactorizeP :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                   -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                   -> CholmodFactorizationStyle
                   -> MatrixSymmetry
                   -> SLA.SpMatrix Double -- ^ matrix to decompose
                   -> IO () -- ^ The factor in ForeignPtr Factor will be updated :(
spMatrixFactorizeP fpc fpf cfs ms smX = do
  triplet <- spMatrixToTriplet fpc ms smX
  when debug $ printTriplet (CH.fPtr triplet) "spMatrixCholesky" fpc
  sparse <- CH.tripletToSparse triplet fpc
  when debug $ printSparse (CH.fPtr sparse) "spMatrixCholesky" fpc
  case cfs of
    FactorizeA -> CH.factorize sparse fpf fpc
    FactorizeAtAPlusBetaI x -> do
      let (_, ncols) = SLA.dimSM smX
          beta = [x,0]
      withForeignPtr (CH.fPtr sparse) $ \pm -> do
        alloca $ \pbeta -> do
          pokeArray pbeta beta
          factorize_p pm pbeta nullPtr 0 fpf fpc
  withForeignPtr fpc $ \pc -> do
    withForeignPtr (CH.fPtr triplet) $ \pt -> CH.triplet_free pt pc
    withForeignPtr (CH.fPtr sparse) $ \ps -> CH.sparse_free ps pc
  return ()


  
-- |  compute the lower-triangular Cholesky factor using the given analysis stored in Factor
-- the matrix given here must have the same pattern of non-zeroes as the one used for the
-- analysis
-- NB: THis version calls unsafePerformIO which, I think, is okay because the function is referentially transparent.
-- At least in a single threaded environment. Oy.
-- TODO: But I am not sure all the allocation and freeing is really correct.
unsafeSpMatrixCholesky ::  ForeignPtr CH.Common -- ^ CHOLMOD environment
                       -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                       -> MatrixSymmetry
                       -> SLA.SpMatrix Double -- ^ matrix to decompose                       
                       -> SLA.SpMatrix Double -- ^ lower-triangular Cholesky factor
unsafeSpMatrixCholesky fpc fpf ms smX = unsafePerformIO $ spMatrixCholesky fpc fpf ms smX
{-# NOINLINE unsafeSpMatrixCholesky #-} 

unsafeSpMatrixFactorize ::  ForeignPtr CH.Common -- ^ CHOLMOD environment
                        -> ForeignPtr CH.Factor -- ^ result of CHOLMOD analysis
                        -> MatrixSymmetry
                        -> SLA.SpMatrix Double -- ^ matrix to decompose
                        -> () -- ^ The factor in ForeignPtr Factor will be updated :(
unsafeSpMatrixFactorize fpc fpf sm smX = unsafePerformIO $ spMatrixFactorize fpc fpf sm smX
{-# NOINLINE unsafeSpMatrixFactorize #-}

unsafeSolveSparse :: ForeignPtr CH.Common -- ^ CHOLMOD environment
                  -> ForeignPtr CH.Factor -- ^ contains P and LL' (or LDL')  
                  -> SolveSystem -- ^ which system to solve
                  -> SLA.SpMatrix Double -- ^ RHS
                  -> SLA.SpMatrix Double -- ^ solutions
unsafeSolveSparse fpc fpf ss smB = unsafePerformIO $ solveSparse fpc fpf ss smB
{-# NOINLINE unsafeSolveSparse #-}
  
factorPermutation :: ForeignPtr CH.Factor -> IO (VS.MVector s CInt)
factorPermutation ffp = withForeignPtr ffp $ \fp -> do
  n <- factorN_L fp :: IO CSize
  pi <- factorPermutationL fp :: IO (Ptr CInt)
  pifp <- newForeignPtr_ pi :: IO (ForeignPtr CInt)
  return $ MVS.unsafeFromForeignPtr0 pifp (fromIntegral n)

-- | retrieve the factor permutation as an SpMatrix
factorPermutationSM :: Num a => ForeignPtr CH.Factor -> IO (SLA.SpMatrix a)
factorPermutationSM ffp = do
  permMV <- factorPermutation ffp
  permL <- readv permMV
  return $ SLA.permutationSM (MVS.length permMV) (fmap fromIntegral permL)

unsafeFactorPermutationSM :: Num a => ForeignPtr CH.Factor -> SLA.SpMatrix a
unsafeFactorPermutationSM = unsafePerformIO . factorPermutationSM
{-# NOINLINE unsafeFactorPermutationSM #-}

-- | retrieve the Cholesky Factor as an SLA.SpMatrix
-- | NB: This undoes the factorization in Factor
choleskyFactorSM :: ForeignPtr CH.Factor -> ForeignPtr CH.Common -> IO (SLA.SpMatrix Double)
choleskyFactorSM ffp cfp = withForeignPtr ffp $ \fp -> do
  withForeignPtr cfp $ \cp -> do
    when debug $ printFactor ffp "Before" cfp
    changeFactorL (fromIntegral $ CH.unXType CH.xtReal) 1 0 0 0 fp cp
    when debug $ printFactor ffp "After" cfp
    pSparse <- factorToSparseL fp cp :: IO (Ptr CH.Sparse)
    pTriplet <- sparseToTripletL pSparse cp :: IO (Ptr CH.Triplet)
    cf <- tripletToSpMatrix pTriplet
    CH.triplet_free pTriplet cp
    CH.sparse_free pSparse cp
    return cf

unsafeCholeskyFactorSM :: ForeignPtr CH.Factor -> ForeignPtr CH.Common -> SLA.SpMatrix Double
unsafeCholeskyFactorSM fpf fpc = unsafePerformIO $ choleskyFactorSM fpf fpc
                       
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
  rowIndices <- readv (MVS.unsafeFromForeignPtr0 rowIndicesFP (fromIntegral nZmax))
  colIndices <- readv (MVS.unsafeFromForeignPtr0 colIndicesFP (fromIntegral nZmax))
  vals <- readv (MVS.unsafeFromForeignPtr0 valsFP (fromIntegral nZmax))
  let triplets = zip3 (fmap fromIntegral rowIndices) (fmap fromIntegral colIndices) (fmap realToFrac vals)
      sm = SLA.fromListSM (fromIntegral nRows, fromIntegral nCols) triplets
  when debug $ SLA.prd sm    
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

printCommon' :: T.Text -> Ptr CH.Common -> IO ()
printCommon' t pc = do
  tCS <- newCString (T.unpack t)
  printCommonL tCS pc

printFactor' :: Ptr CH.Factor -> T.Text -> ForeignPtr CH.Common -> IO ()
printFactor' pf t fpc = do
  withForeignPtr fpc $ \pc -> do
    tSC <- newCString (T.unpack t)
    printFactorL pf tSC pc

printTriplet' :: Ptr CH.Triplet -> T.Text -> ForeignPtr CH.Common -> IO ()
printTriplet' pt t fpc = do
  withForeignPtr fpc $ \pc -> do
    tSC <- newCString (T.unpack t)
    printTripletL pt tSC pc

printSparse' :: Ptr CH.Sparse -> T.Text -> ForeignPtr CH.Common -> IO ()
printSparse' ps t fpc = do
  withForeignPtr fpc $ \pc -> do
    tSC <- newCString (T.unpack t)
    printSparseL ps tSC pc

setFinalLL :: Int -> ForeignPtr CH.Common -> IO ()
setFinalLL n fpc = withForeignPtr fpc $ \pc -> setFinalLLL n pc

factorize_p :: Ptr CH.Sparse -> Ptr Double -> Ptr CInt -> CSize -> ForeignPtr Factor -> ForeignPtr Common -> IO ()
factorize_p ps pb i s fpf fpc = do
  withForeignPtr fpc $ \pc -> do
    withForeignPtr fpf $ \pf -> factorize_pL ps pb i s pf pc
      
readv :: MVS.Storable a => MVS.IOVector a -> IO [a]
readv v = sequence [ MVS.read v i | i <- [0 .. (MVS.length v) - 1] ]

writev :: MVS.Storable a => MVS.IOVector a -> [a] -> IO ()
writev v xs =
  sequence_ [ MVS.write v i x | (i, x) <- zip [0 .. (Prelude.length xs - 1)] xs ]

copyFactor :: ForeignPtr CH.Factor -> ForeignPtr CH.Common -> IO (ForeignPtr CH.Factor)
copyFactor fpF fpC = do
  withForeignPtr fpC $ \pC -> do
    withForeignPtr fpF $ \pF -> do
      pFactorCopy <- copyFactorL pF pC
      newForeignPtrEnv CH.factor_free_ptr pC pFactorCopy

-- | the N of the N x N factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_n"
  factorN_L :: Ptr CH.Factor -> IO CSize

-- | the array of indices representing the permutation matrix in the factor
foreign import ccall unsafe "cholmod_extras.h cholmod_factor_permutation"
  factorPermutationL :: Ptr CH.Factor -> IO (Ptr CInt)

foreign import ccall unsafe "cholmod_extras.h cholmod_set_final_ll"
   setFinalLLL :: Int -> Ptr CH.Common -> IO () 

foreign import ccall unsafe "cholmod.h cholmod_factor_to_sparse"
  factorToSparseL :: Ptr CH.Factor -> Ptr CH.Common -> IO (Ptr CH.Sparse) 

foreign import ccall unsafe "cholmod.h cholmod_sparse_to_triplet"
  sparseToTripletL :: Ptr CH.Sparse -> Ptr CH.Common -> IO (Ptr CH.Triplet)

foreign import ccall unsafe "cholmod.h cholmod_change_factor"
  changeFactorL :: Int -> Int -> Int -> Int -> Int -> Ptr CH.Factor -> Ptr CH.Common -> IO CInt 

foreign import ccall unsafe "cholmod.h cholmod_solve"
  denseSolveL :: CInt -> Ptr CH.Factor -> Ptr CH.Dense -> Ptr CH.Common -> IO (Ptr CH.Dense)

foreign import ccall unsafe "cholmod.h cholmod_spsolve"
  sparseSolveL :: CInt -> Ptr CH.Factor -> Ptr CH.Sparse -> Ptr CH.Common -> IO (Ptr CH.Sparse)

foreign import ccall unsafe "cholmod.h cholmod_print_common"
  printCommonL :: CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_factor"
  printFactorL :: Ptr CH.Factor -> CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_triplet"
  printTripletL :: Ptr CH.Triplet -> CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_print_sparse"
  printSparseL :: Ptr CH.Sparse -> CString -> Ptr CH.Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_factorize_p"
  factorize_pL :: Ptr CH.Sparse -> Ptr Double -> Ptr CInt -> CSize -> Ptr Factor -> Ptr Common -> IO ()

foreign import ccall unsafe "cholmod.h cholmod_copy_factor"
  copyFactorL :: Ptr CH.Factor -> Ptr CH.Common -> IO (Ptr CH.Factor)



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


