{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE PolyKinds           #-}
{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE QuasiQuotes         #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TupleSections       #-}

module DataFrames where

import qualified Frames                        as F
import qualified Frames.CSV                    as F
import qualified Frames.InCore                 as FI
import qualified Data.Vinyl                    as V
--import qualified Data.Vinyl.TypeLevel          as V

import qualified Numeric.LinearAlgebra         as LA
--import qualified Data.Sparse.SpMatrix          as SLA

import qualified Data.Map                      as M
import           Data.Maybe                     ( catMaybes )
import qualified Data.Text                     as T
import           Data.Text                      ( Text )

import qualified Pipes                         as P
import qualified Pipes.Prelude                 as P
import qualified Control.Foldl                 as FL
import qualified Control.Monad.State           as ST
import qualified Data.Vector                   as VB

railCSV = "data/Rail.csv"
F.tableTypes "RailRow" "data/Rail.csv"

sleepStudyCSV = "data/SleepStudy.csv"
F.tableTypes "SleepStudyRow" "data/SleepStudy.csv"

oatsCSV = "data/Oats.csv"
F.tableTypes "OatsRow" "data/Oats.csv"

loadToFrame
  :: forall rs
   . (F.ReadRec rs, FI.RecVec rs, V.RMap rs)
  => F.ParserOptions
  -> FilePath
  -> (F.Record rs -> Bool)
  -> IO (F.FrameRec rs)
loadToFrame po fp filterF = do
  let producer = F.readTableOpt po fp P.>-> P.filter filterF
  frame <- F.inCoreAoS producer
  let reportRows :: Foldable f => f x -> FilePath -> IO ()
      reportRows f fn =
        putStrLn
          .  T.unpack
          $  T.pack (show $ FL.fold FL.length f)
          <> " rows in "
          <> T.pack fn
  reportRows frame fp
  return frame


defaultLoadToFrame
  :: forall rs
   . (F.ReadRec rs, FI.RecVec rs, V.RMap rs)
  => FilePath
  -> (F.Record rs -> Bool)
  -> IO (F.FrameRec rs)
defaultLoadToFrame = loadToFrame F.defaultParser


--type OrdEqSnd t = (Ord (V.Snd t), Eq (V.Snd t))
--categorizeFields :: V.AllConstrained OrdEqSnd rs => (F.Record rs -> [Int])



lmePrepFrameOne
  :: (Ord b, Eq b)
  => (F.Record rs -> Double) -- ^ observations
  -> (F.Record rs -> LA.Vector Double) -- ^ predictors
  -> (F.Record rs -> b) -- ^ classifier
  -> FL.Fold
       (F.Record rs)
       ( LA.Vector Double
       , LA.Matrix Double
       , (VB.Vector (VB.Vector Int), Int)
       ) -- ^ (X,y,(row-classifier, size of class))
lmePrepFrameOne observationF predictorF classF =
  let foldMapF
        :: Ord b
        => ST.State Int (M.Map b Int)
        -> b
        -> ST.State Int (M.Map b Int)
      foldMapF mM b = do
        m <- mM
        case M.lookup b m of
          Nothing -> do
            nextIndex <- ST.get
            let m' = M.insert b nextIndex m
            ST.put (nextIndex + 1)
            return m'
          _ -> return m
      foldClassMapM = FL.Fold foldMapF (return M.empty) id
      makeClassifiers m bs =
        ( (VB.fromList . fmap (VB.fromList . pure) . catMaybes $ fmap
            (\b -> M.lookup b m)
            bs
          )
        , M.size m
        )
      foldObs   = fmap LA.fromList $ FL.premap observationF FL.list
      foldPred  = fmap LA.fromRows $ FL.premap predictorF FL.list
      foldClass = FL.premap classF FL.list
      foldMapM  = FL.premap classF foldClassMapM
      g (vY, mX, mM, bs) = (vY, mX, makeClassifiers (ST.evalState mM 0) bs)
  in  fmap g $ ((,,,) <$> foldObs <*> foldPred <*> foldMapM <*> foldClass)


lmePrepFrameTwo
  :: (Ord b, Eq b, Ord c, Eq c)
  => (F.Record rs -> Double) -- ^ observations
  -> (F.Record rs -> LA.Vector Double) -- ^ predictors
  -> (F.Record rs -> b) -- ^ classifier1
  -> (F.Record rs -> c) -- ^ classifier2
  -> FL.Fold
       (F.Record rs)
       ( LA.Vector Double
       , LA.Matrix Double
       , (VB.Vector (VB.Vector Int), Int, Int)
       ) -- ^ (X,y,(row-classifier, size of class 1, size of class 2))
lmePrepFrameTwo observationF predictorF classbF classcF =
  let
    foldMapF
      :: Ord z => ST.State Int (M.Map z Int) -> z -> ST.State Int (M.Map z Int)
    foldMapF mM z = do
      m <- mM
      case M.lookup z m of
        Nothing -> do
          nextIndex <- ST.get
          let m' = M.insert z nextIndex m
          ST.put (nextIndex + 1)
          return m'
        _ -> return m
    foldClassMapM = FL.Fold foldMapF (return M.empty) id
    makeClassifiers mb mc bs cs =
      let x = zip bs cs
          f (b, c) = do
            bCat <- M.lookup b mb
            cCat <- M.lookup c mc
            return $ VB.fromList [bCat, cCat]
      in  ((VB.fromList . catMaybes $ fmap f x), M.size mb, M.size mc)
    foldObs    = fmap LA.fromList $ FL.premap observationF FL.list
    foldPred   = fmap LA.fromRows $ FL.premap predictorF FL.list
    foldClassb = FL.premap classbF FL.list
    foldMapbM  = FL.premap classbF foldClassMapM
    foldClassc = FL.premap classcF FL.list
    foldMapcM  = FL.premap classcF foldClassMapM
    g (vY, mX, mMb, bs, mMc, cs) =
      (vY, mX, makeClassifiers (ST.evalState mMb 0) (ST.evalState mMc 0) bs cs)
  in
    fmap g
      $ (   (,,,,,)
        <$> foldObs
        <*> foldPred
        <*> foldMapbM
        <*> foldClassb
        <*> foldMapcM
        <*> foldClassc
        )


