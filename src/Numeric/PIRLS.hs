{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.PIRLS where

import qualified Polysemy                      as P
import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )
import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import qualified Data.Sequence                 as Seq
import qualified Data.Vector                   as VB

type Levels = VB.Vector (Int, Bool, Maybe (VB.Vector Bool)) -- level sizes and effects

-- classify row into its levels
type RowClassifier = Int -> VB.Vector Int

-- for group k, what group effects are we modeling?
-- first Bool is for intercept, then (optional) vector,
-- of same length as X has columns, to indicate which
-- predictors in X get random slopes
-- NB:  If X has a constant column, there is redundancy
-- here between the Bool in the tuple and the first bool
-- in the vector.
type LevelEffectsSpec = [(Bool, Maybe (LA.Vector Bool))]

makeZ
  :: (LA.Container LA.Vector a, RealFrac a, MonadIO (P.Sem r), a ~ Double)
  => LA.Matrix a
  -> Levels
  -> RowClassifier
  -> P.Sem r (SLA.SpMatrix a)
makeZ mX levels rc = do
  let (nO, nP) = LA.size mX
      k        = VB.length levels -- number of levels
      levelSize n = let (s, _, _) = levels VB.! n in s
      colsForLevel (qL, b, vbM) =
        qL * ((if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM))
      q = FL.fold FL.sum $ fmap colsForLevel levels -- total number of columns in Z
{-      
  liftIO $ do
    putStrLn "X="
    LA.disp 2 mX
    putStrLn $ "We have " ++ show k ++ " levels"
    putStrLn $ "And should get " ++ show q ++ " columns in Z"
-}
    -- build list of items to put in Z, level by level
  let
    obsIndices = Seq.fromFunction nO id
    intercept startingCol level =
      let spIndex obsIndex = (obsIndex, rc obsIndex VB.! level)
          spEntry obsIndex =
            let spi = spIndex obsIndex in (fst spi, startingCol + snd spi, 1)
      in  fmap spEntry obsIndices
    slope startingCol level predIndex =
      let spIndex obsIndex = (obsIndex, rc obsIndex VB.! level)
          spValue obsIndex = mX `LA.atIndex` (obsIndex, predIndex)
          spEntry obsIndex =
            let spi = spIndex obsIndex
            in  (fst spi, startingCol + snd spi, spValue obsIndex)
      in  fmap spEntry obsIndices
    zEntriesForSlope startingCol level vb =
      let predIndices = fmap fst $ VB.filter snd $ VB.zip
            (VB.generate (VB.length vb) id)
            vb
          startingCols = VB.generate
            (VB.length predIndices)
            (\n -> startingCol + n * (levelSize level))
          newStart = startingCol + (VB.length predIndices * (levelSize level))
          entriesF (pi, sc) = slope sc level pi
      in  ( FL.fold FL.mconcat $ fmap entriesF $ VB.zip predIndices startingCols
          , newStart
          )
    zEntriesForLevel startingCol level (qL, b, vbM) =
      let (interceptZs, newStartI) = if b
            then (intercept startingCol level, startingCol + qL)
            else (Seq.empty, startingCol)
          (slopeZs, newStartS) = case vbM of
            Nothing -> (Seq.empty, newStartI)
            Just vb -> zEntriesForSlope newStartI level vb
      in  (interceptZs <> slopeZs, newStartS)
    zFold = FL.Fold
      (\(zs, sc, l) x ->
        let (newZs, newStart) = zEntriesForLevel sc l x
        in  (zs <> newZs, newStart, l + 1)
      )
      (Seq.empty, 0, 0)
      id
    (zEntries, numCols, _) = FL.fold zFold levels
{-    
  liftIO $ do
    putStrLn $ "numCols=" ++ show numCols
    putStrLn $ "entries=" ++ show zEntries
-}
  return $ SLA.fromListSM (nO, q) zEntries
