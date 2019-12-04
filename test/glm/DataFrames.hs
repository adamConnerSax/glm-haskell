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

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.ProblemTypes      as GLM

import qualified Frames                        as F
import qualified Frames.CSV                    as F
import qualified Frames.InCore                 as FI
import qualified Data.Vinyl                    as V

import qualified Numeric.LinearAlgebra         as LA

import qualified Data.Array                    as A
import qualified Data.List                     as L
import qualified Data.Map                      as M
import           Data.Maybe                     ( catMaybes
                                                , fromMaybe
                                                )
import qualified Data.Text                     as T
import           Data.Text                      ( Text )

import qualified Pipes                         as P
import qualified Pipes.Prelude                 as P
import qualified Control.Foldl                 as FL
import qualified Control.Monad.State           as ST
import qualified Data.Vector                   as VB

railCSV = "data/Rail.csv"
F.tableTypes "RailRow" "data/Rail.csv"

type RailEffect = GLM.WithIntercept ()
railPredictor :: RailRow -> () -> Double
railPredictor _ _ = undefined -- we need something of this type but if this gets called something has gone wrong!

data RailGroup = RG_Rail deriving (Show, Enum, Bounded, Eq, Ord, A.Ix)
type instance GLM.GroupKey RailGroup = T.Text

railGroupLabels :: RailRow -> RailGroup -> T.Text
railGroupLabels row group = case group of
  RG_Rail -> F.rgetField @Rail row



--

sleepStudyCSV = "data/SleepStudy.csv"
F.tableTypes "SleepStudyRow" "data/SleepStudy.csv"

data SleepStudyPredictor = SleepStudyDays deriving (Eq, Ord, Enum, Bounded, Show, A.Ix)
type SleepStudyEffect = GLM.WithIntercept SleepStudyPredictor

sleepStudyPredictor :: SleepStudyRow -> SleepStudyPredictor -> Double
sleepStudyPredictor r _ = realToFrac $ F.rgetField @Days r

data SleepStudyGroup = SSG_Subject deriving (Show, Enum, Bounded, Eq, Ord, A.Ix)
type instance GLM.GroupKey SleepStudyGroup = T.Text

sleepStudyGroupLabels :: SleepStudyRow -> SleepStudyGroup -> T.Text
sleepStudyGroupLabels row group = case group of
  SSG_Subject -> T.pack . show $ F.rgetField @Subject row


--

oatsCSV = "data/Oats.csv"
F.tableTypes "OatsRow" "data/Oats.csv"

data OatsPredictor = OatsNitro deriving (Eq, Ord, Enum, Bounded, Show, A.Ix)
type OatsEffect = GLM.WithIntercept OatsPredictor

getOatsPredictor :: OatsRow -> OatsPredictor -> Double
getOatsPredictor r _ = F.rgetField @Nitro r

data OatsGroup = OG_Block | OG_VarietyBlock deriving (Show, Enum, Bounded, Eq, Ord, A.Ix)
type instance GLM.GroupKey OatsGroup = T.Text

oatsGroupLabels :: OatsRow -> OatsGroup -> T.Text
oatsGroupLabels row group = case group of
  OG_Block        -> F.rgetField @Block row
  OG_VarietyBlock -> F.rgetField @Variety row <> "_" <> F.rgetField @Block row



--

cbppCSV = "data/cbpp.csv"
F.tableTypes "CbppRow" "data/cbpp.csv"

data CbppPredictor = CbppPeriod deriving (Eq, Ord, Enum, Bounded, Show, A.Ix)
type CbppEffect = GLM.WithIntercept CbppPredictor

getCbppPredictor :: CbppRow -> CbppPredictor -> Double
getCbppPredictor r _ = realToFrac $ F.rgetField @Period r

data CbppGroup = CbppHerd deriving (Show, Enum, Bounded, Eq, Ord, A.Ix)
type instance GLM.GroupKey CbppGroup = T.Text

cbppGroupLabels :: CbppRow -> CbppGroup -> T.Text
cbppGroupLabels r _ = T.pack $ show $ F.rgetField @Herd r
--

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


lmePrepFrame
  :: forall p g rs
   . (GLM.PredictorC p, GLM.GroupC g)
  => (F.Record rs -> Double) -- ^ observations
  -> GLM.FixedEffects p
  -> IS.IndexedSet g
  -> (F.Record rs -> p -> Double) -- ^ predictors
  -> (F.Record rs -> g -> GLM.GroupKey g)  -- ^ classifiers
  -> FL.Fold
       (F.Record rs)
       ( LA.Vector Double
       , LA.Matrix Double
       , Either Text (GLM.RowClassifier g)
       ) -- ^ (X,y,(row-classifier, size of class))
lmePrepFrame observationF fe groupIndices getPredictorF classifierLabelF
  = let
      makeInfoVector
        :: M.Map g (M.Map (GLM.GroupKey g) Int)
        -> M.Map g (GLM.GroupKey g)
        -> Either Text (VB.Vector (GLM.ItemInfo g))
      makeInfoVector indexMaps groupKeys =
        let
          g (grp, groupKey) =
            GLM.ItemInfo
              <$> (   maybe
                      (Left $ "Failed on " <> (T.pack $ show (grp, groupKey)))
                      Right
                  $   M.lookup grp indexMaps
                  >>= M.lookup groupKey
                  )
              <*> pure groupKey
        in  fmap VB.fromList $ traverse g $ M.toList groupKeys
      makeRowClassifier
        :: Traversable f
        => M.Map g (M.Map (GLM.GroupKey g) Int)
        -> f (M.Map g (GLM.GroupKey g))
        -> Either Text (GLM.RowClassifier g)
      makeRowClassifier indexMaps labels = do
        let sizes = fmap M.size indexMaps
        indexed <- traverse (makeInfoVector indexMaps) labels
        return $ GLM.RowClassifier groupIndices
                                   sizes
                                   (VB.fromList $ FL.fold FL.list indexed)
                                   indexMaps
      getPredictorF' _   GLM.Intercept     = 1
      getPredictorF' row (GLM.Predictor x) = getPredictorF row x
      predictorF row = LA.fromList $ case fe of
        GLM.FixedEffects indexedFixedEffects ->
          fmap (getPredictorF' row) $ IS.members indexedFixedEffects
        GLM.InterceptOnly -> [1]
      getClassifierLabels :: F.Record rs -> M.Map g (GLM.GroupKey g)
      getClassifierLabels r =
        M.fromList $ fmap (\g -> (g, classifierLabelF r g)) $ IS.members
          groupIndices
      foldObs   = fmap LA.fromList $ FL.premap observationF FL.list
      foldPred  = fmap LA.fromRows $ FL.premap predictorF FL.list
      foldClass = FL.premap getClassifierLabels FL.list
      g (vY, mX, ls) =
        ( vY
        , mX
        , makeRowClassifier
          (snd $ ST.execState (addAll ls) (M.empty, M.empty))
          ls
        )
    in
      fmap g $ ((,,) <$> foldObs <*> foldPred <*> foldClass)



addOne
  :: GLM.GroupC g
  => (g, GLM.GroupKey g)
  -> ST.State (M.Map g Int, M.Map g (M.Map (GLM.GroupKey g) Int)) ()
addOne (grp, label) = do
  (nextIndexMap, groupIndexMaps) <- ST.get
  let groupIndexMap = fromMaybe M.empty $ M.lookup grp groupIndexMaps
  case M.lookup label groupIndexMap of
    Nothing -> do
      let index         = fromMaybe 0 $ M.lookup grp nextIndexMap
          nextIndexMap' = M.insert grp (index + 1) nextIndexMap
          groupIndexMaps' =
            M.insert grp (M.insert label index groupIndexMap) groupIndexMaps
      ST.put (nextIndexMap', groupIndexMaps')
      return ()
    _ -> return ()

addMany
  :: (GLM.GroupC g, Traversable h)
  => h (g, GLM.GroupKey g)
  -> ST.State (M.Map g Int, M.Map g (M.Map (GLM.GroupKey g) Int)) ()
addMany x = traverse addOne x >> return ()

addAll
  :: GLM.GroupC g
  => [M.Map g (GLM.GroupKey g)]
  -> ST.State (M.Map g Int, M.Map g (M.Map (GLM.GroupKey g) Int)) ()
addAll x = traverse (addMany . M.toList) x >> return ()

