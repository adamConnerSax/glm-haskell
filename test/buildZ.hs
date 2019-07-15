{-# LANGUAGE ScopedTypeVariables #-}
module Main where

import           Numeric.PIRLS

import qualified Polysemy                      as P
import qualified Control.Foldl                 as FL
import           Control.Monad.IO.Class         ( MonadIO )
import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import qualified Data.Sparse.Common            as SLA
import qualified Data.Sequence                 as Seq
import qualified Data.Vector                   as VB


xM :: LA.Matrix Double = LA.matrix 2 [1, 1, 1, 2, 1, 3, 1, 4]
levels =
  VB.fromList [(2, True, Just $ VB.fromList [False, True]), (2, True, Nothing)]
rowClassifier n = fmap VB.fromList [[0, 0], [0, 1], [1, 0], [1, 1]] !! n

main :: IO ()
main = do
  z <- P.runM $ makeZ xM levels rowClassifier
  SLA.prd0 z

