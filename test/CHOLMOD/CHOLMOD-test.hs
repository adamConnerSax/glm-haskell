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

  let smX = SLA.fromListSM
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
  let ms = SquareSymmetricLower
  putStrLn $ "smX="
  LA.disp 1 $ asDense smX
--  mt                        <- spMatrixToTriplet c ms smX
  (cholmodFactor, cholPerm) <- spMatrixAnalyze c ms smX
  cholL_CM                  <- spMatrixCholesky c cholmodFactor ms smX
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
