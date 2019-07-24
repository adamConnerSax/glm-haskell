{-# LANGUAGE  FlexibleInstances    #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverloadedStrings #-}
-- * base
import           Foreign                        ( ForeignPtr
                                                , Storable
                                                , withForeignPtr
                                                )
import           Foreign.C.Types                ( CDouble
                                                , CInt
                                                , CSize
                                                )

-- * vector
import qualified Data.Vector.Storable.Mutable  as V
import qualified Data.Vector.Storable          as V
                                                ( freeze )

-- * hmatrix
import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.LinearAlgebra.HMatrix as LA

-- * sparse-linear-algebra
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
-- * cholmod
import           Numeric.LinearAlgebra.CHOLMOD.CholmodXFaceLow
import           Numeric.LinearAlgebra.CHOLMOD.CholmodXFace



-- * glm additions to cholmod
import           Numeric.LinearAlgebra.CHOLMOD.CholmodExtras

-- * glm
import           Numeric.MixedModel
--------------------------------------------------------------------------------

main :: IO ()
main = do

  c <- allocCommon :: IO (ForeignPtr Common)
  startC c
  setFinalLL 1 c
  printCommon "" c
{-
  let n  = 5 :: CSize
      nz = (n + 1) * n `div` 2 :: CSize
      nn = fromIntegral n :: CSize


  at <-
    allocTriplet n n nz stSquareSymmetricLower xtReal c :: IO (Matrix Triplet)

  atNRow  <- getNRow at :: IO CSize
  atNCol  <- getNCol at :: IO CSize
  atNZMax <- getNZMax at :: IO CSize

  putStrLn ""
  putStrLn $ "at: (" ++ show atNRow ++ " x " ++ show atNCol ++ ")"
  putStrLn $ "nzmax: " ++ show atNZMax

  iv <- tripletGetRowIndices at :: IO (V.IOVector CInt)
  jv <- tripletGetColIndices at :: IO (V.IOVector CInt)
  xt <- tripletGetX at


  let nni = fromIntegral nn

  let
    ij       = [ (i, j) | i <- [0 .. nni - 1], j <- [0 .. i] ] :: [(CInt, CInt)]
    (ia, ja) = unzip ij :: ([CInt], [CInt])

    -- the elements of the matrix are (lower half)
    --
    --  11.0
    --  21.0  22.0
    --  31.0  32.0  33.0
    --  41.0  42.0  43.0  44.0
    --  51.0  52.0  53.0  54.0  55.0

    xp =
      [11, 21, 22, 31, 32, 33, 41, 42, 43, 44, 51, 52, 53, 54, 55] :: [CDouble]

  writev iv ia
  writev jv ja
  writev xt xp

  let nnz = fromIntegral $ Prelude.length xp :: CSize

  tripletSetNNZ at nnz


  as     <- tripletToSparse at c :: IO (Matrix Sparse)
  b      <- ones n 1 xtReal c :: IO (Matrix Dense)
  l      <- analyze as c :: IO (ForeignPtr Factor)

  permM  <- factorPermutation l
  permSM <- factorPermutationSM l
  let permLength = V.length permM
  permV <- V.freeze permM
  putStrLn $ "perm length=" ++ show permLength
  putStrLn $ "perm=" ++ show permV
  putStrLn $ "permSM="
  LA.disp 1 $ asDense permSM

  factorize as l c
  choleskyL_SM <- choleskyFactorSM l c
  LA.disp 2 $ asDense choleskyL_SM
  factorize as l c
  x     <- solve stA l b c :: IO (Matrix Dense)

  xNRow <- getNRow x
  xNCol <- getNCol x

  xv    <- getX x :: IO (V.IOVector CDouble)

  putStrLn ""
  putStrLn $ "xv: (" ++ show xNRow ++ " x " ++ show xNCol ++ ")"

  xl <- readv xv
  mapM_ (putStrLn . show) xl

  r <- denseCopy b c

  let oneL = [1, 0] :: [CDouble]
      m1L  = [-1, 0] :: [CDouble]

  _    <- sdMult as noTranspose m1L oneL x r c

  norm <- denseNorm r infNorm c :: IO CDouble

  putStrLn $ "norm: " ++ show norm


  putStrLn $ "Now via SpMatrix"
  let mX = SLA.fromListSM
        (5, 5)
        [ (0, 0, 11)
        , (1, 0, 21)
        , (2, 0, 31)
        , (3, 0, 41)
        , (4, 0, 51)
        , (1, 1, 22)
        , (2, 1, 32)
        , (3, 1, 42)
        , (4, 1, 52)
        , (2, 2, 33)
        , (3, 2, 43)
        , (4, 2, 53)
        , (3, 3, 44)
        , (4, 3, 54)
        , (4, 4, 55)
        ]
  let mX' = SLA.fromListSM
        (5, 5)
        [ (0, 0, 11)
        , (1, 0, 121)
        , (2, 0, 31)
        , (3, 0, 41)
        , (4, 0, 51)
        , (1, 1, 122)
        , (2, 1, 32)
        , (3, 1, 42)
        , (4, 1, 52)
        , (2, 2, 33)
        , (3, 2, 43)
        , (4, 2, 53)
        , (3, 3, 44)
        , (4, 3, 54)
        , (4, 4, 55)
        ]
-}
  let smY = SLA.fromListSM
        (3, 3)
        [ (0, 0, 2)
        , (0, 1, -1)
        , (0, 2, 0)
        , (1, 0, -1)
        , (1, 1, 2)
        , (1, 2, -1)
        , (2, 0, 0)
        , (2, 1, -1)
        , (2, 2, 2)
        ]

  putStrLn $ "smY="
  LA.disp 1 $ asDense smY
  mt   <- spMatrixToTriplet c smY
  smY' <- withForeignPtr (fPtr mt) tripletToSpMatrix
  putStrLn $ "smY'="
  LA.disp 1 $ asDense smY'
  (cholmodFactor, cholPerm) <- spMatrixAnalyze c smY
  cholL_CM                  <- spMatrixCholesky c cholmodFactor smY

--      cholL'_CM = unsafeSpMatrixCholesky c cholmodFactor mX'
--      cholL_HM  = LA.tr $ LA.chol $ LA.trustSym $ LA.tr $ asDense mX
--      cholL'_HM = LA.tr $ LA.chol $ LA.trustSym $ LA.tr $ asDense mX'
  putStrLn $ "perm="
  LA.disp 0 $ asDense cholPerm

  putStrLn $ "L (CHOLMOD)="
  LA.disp 2 $ asDense cholL_CM
  putStrLn $ "LLt="
  LA.disp 2 $ asDense (cholL_CM SLA.## (SLA.transposeSM cholL_CM))

  return ()
{-
writev :: (Storable a) => V.IOVector a -> [a] -> IO ()
writev v xs =
  sequence_ [ V.write v i x | (i, x) <- zip [0 .. (Prelude.length xs - 1)] xs ]


readv :: (Storable a) => V.IOVector a -> IO [a]
readv v = sequence [ V.read v i | i <- [0 .. (V.length v) - 1] ]
-}
