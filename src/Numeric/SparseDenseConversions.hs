{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.SparseDenseConversions where

import qualified Control.Foldl                 as FL

import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Numeric.LinearAlgebra.Class   as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA

import qualified Numeric.LinearAlgebra         as LA

import qualified Data.List                     as L
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

toDenseMatrix :: SLA.SpMatrix Double -> LA.Matrix Double
toDenseMatrix sp =
  LA.matrix (snd $ SLA.dimSM sp) $ fmap (\(_, _, x) -> x) $ SLA.toDenseListSM sp

toDenseVector :: (Num a, VS.Storable a) => SLA.SpVector a -> LA.Vector a
toDenseVector = VB.convert . SLA.toVectorDense

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

toSparseVector
  :: (LA.Container LA.Vector a, RealFrac a) => LA.Vector a -> SLA.SpVector a
toSparseVector v = SLA.fromListDenseSV (LA.size v) $ LA.toList v

svColumnToSM :: SLA.SpVector a -> SLA.SpMatrix a
svColumnToSM svX =
  SLA.fromListSM (SLA.dim svX, 1) $ fmap (\(r, x) -> (r, 0, x)) $ SLA.toListSV
    svX

svRowToSM :: SLA.SpVector a -> SLA.SpMatrix a
svRowToSM svX =
  SLA.fromListSM (1, SLA.dim svX) $ fmap (\(c, x) -> (0, c, x)) $ SLA.toListSV
    svX
