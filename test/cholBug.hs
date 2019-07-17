{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import qualified Data.Sparse.Common            as SLA
import qualified Data.List                     as L

toSparseMatrix
  :: (LA.Container LA.Vector a, RealFrac a) => LA.Matrix a -> SLA.SpMatrix a
toSparseMatrix x
  = let
      (r, c)     = LA.size x
      colIndices = take c $ L.iterate (+ 1) 0
      rowIndices = take r $ L.iterate (+ 1) 0
      indexRow rI =
        fmap (\(cI, val) -> (rI, cI, val)) . zip colIndices . LA.toList
      indexedRows =
        concat . fmap (\(rI, rV) -> indexRow rI rV) . zip rowIndices . LA.toRows
    in
      SLA.fromListSM (r, c) $ indexedRows x

asDense sp =
  LA.matrix (snd $ SLA.dimSM sp) $ fmap (\(_, _, x) -> x) $ SLA.toDenseListSM sp


cholBug :: IO ()
cholBug = do
  let rowMajor :: [Double] =
        [ 95.08
        , 0
        , 0
        , 0
        , 0
        , 0
        , 16.80
        , -908
        , 0
        , 95.08
        , 0
        , 0
        , 0
        , 0
        , 16.80
        , -532
        , 0
        , 0
        , 95.08
        , 0
        , 0
        , 0
        , 16.80
        , -1422
        , 0
        , 0
        , 0
        , 95.08
        , 0
        , 0
        , 16.80
        , -1612
        , 0
        , 0
        , 0
        , 0
        , 95.08
        , 0
        , 16.80
        , -840
        , 0
        , 0
        , 0
        , 0
        , 0
        , 95.08
        , 16.80
        , -1389
        , 16.80
        , 16.80
        , 16.80
        , 16.80
        , 16.80
        , 16.80
        , 18
        , -1197
        , -908
        , -532
        , -1422
        , -1612
        , -840
        , -1389
        , -1197
        , 89105
        ]
      dense   = LA.matrix 8 rowMajor
      sparse1 = toSparseMatrix dense
      sparse2 = SLA.sparsifySM sparse1
      denseCL = LA.tr $ LA.chol $ LA.trustSym $ dense
  sparseCL1 <- SLA.chol sparse1
  sparseCL2 <- SLA.chol sparse2
  putStrLn $ "X="
  LA.disp 2 dense
  putStrLn $ "hmatrix L="
  LA.disp 2 denseCL
  putStrLn $ "sparse-linear-algebra (unsparsified)="
  LA.disp 2 $ asDense sparseCL1
  putStrLn $ "sparse-linear-algebra (sparsified)="
  LA.disp 2 $ asDense sparseCL2
  return ()




