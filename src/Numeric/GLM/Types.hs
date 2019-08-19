{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Numeric.GLM.Types
  ( WithIntercept(..)
  , IndexedEffectSet
  , FixedEffects(..)
  , indexedFixedEffectSet
  , allFixedEffects
  , ItemInfo(..)
  , RowClassifier(..)
  , groupIndices
  , groupSizes
  , rowInfos
  )
where

import qualified Data.IndexedSet               as IS

import qualified Control.Foldl                 as FL
import qualified Data.Array                    as A
import qualified Data.List                     as L
import qualified Data.Map                      as M
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB

data WithIntercept b where
  Intercept :: WithIntercept b
  Predictor :: b -> WithIntercept b

instance Show b => Show (WithIntercept b) where
  show Intercept = "Intercept"
  show (Predictor x) = "Predictor " ++ show x

instance Eq b => Eq (WithIntercept b) where
  Intercept == Intercept = True
  (Predictor x) == (Predictor y) = x == y
  _ == _ = False

instance Ord b => Ord (WithIntercept b) where
  compare Intercept Intercept = EQ
  compare Intercept (Predictor _) = LT
  compare (Predictor _) Intercept = GT
  compare (Predictor x) (Predictor y) = compare x y

instance (Bounded b, Enum b) => Enum (WithIntercept b) where
  fromEnum Intercept = 0
  fromEnum (Predictor x) = 1 + fromEnum x
  toEnum 0 = Intercept
  toEnum n = Predictor $ toEnum (n - 1)
  succ Intercept = Predictor $ toEnum 0
  succ (Predictor x) = Predictor $ succ x
  pred Intercept = error "pred{WithIntercept b} called on \"Intercept\""
  pred x = toEnum $ fromEnum x - 1
  enumFrom Intercept = Intercept : fmap Predictor [minBound..]
  enumFrom (Predictor x) = fmap Predictor [x..]
  enumFromThen = error "enumFromThen{WithIntercept b} called.  Likely using [a,b..] notation.  Not supported for \"WithIntercept\""
  enumFromTo Intercept Intercept = [Intercept]
  enumFromTo Intercept (Predictor x) = Intercept : fmap Predictor [minBound..x]
  enumFromTo (Predictor x) (Predictor y) = fmap Predictor [x..y]
  enumFromTo (Predictor _) Intercept = []

instance Bounded b => Bounded (WithIntercept b) where
  minBound = Intercept
  maxBound = Predictor (maxBound :: b)

instance (Bounded b, Enum b, A.Ix b) => A.Ix (WithIntercept b) where
{-  range (Intercept, Intercept) = [Intercept]
  range (Intercept, (Predictor x)) = Intercept : fmap Predictor (A.range (minBound, x))
  range ((Predictor x), (Predictor y)) = fmap Predictor $ A.range (x, y)
  range (_, _) = [] -}
  range (a,b) = [a..b]
  index (Intercept, _) Intercept = 0
  index (Intercept, Intercept) (Predictor _) = error "Ix{WithIntercept b}.index: Index out of range.  \"index Intercept Intercept (Predictor _)\" called."
  index (Intercept, (Predictor x)) (Predictor y) = 1 + A.index (minBound, x) y
  index ((Predictor _), (Predictor _)) Intercept = error "Ix{WithIntercept b}.index: Index out of range. Intercept is not between any two Predictors."
  index ((Predictor x) ,(Predictor y)) (Predictor z) = A.index (x, y) z
  index (_, _) _ = error "Ix{WithIntercept b}.index: Index out of range."

  inRange (Intercept, _) Intercept = True
  inRange (Intercept, Intercept) (Predictor _) = False
  inRange (Intercept, (Predictor x)) (Predictor y) = A.inRange (minBound, x) y
  inRange ((Predictor _), (Predictor _)) Intercept = False
  inRange ((Predictor x), (Predictor y)) (Predictor z) = A.inRange (x, y) z
  inRange (_, _) _ = False

type IndexedEffectSet b = IS.IndexedSet (WithIntercept b)

{-
data IndexedEffectSet b where
  IndexedEffectSet :: M.Map (WithIntercept b) Int -> M.Map Int (WithIntercept b) -> IndexedEffectSet b
    deriving (Show)

effectSetMembers :: IndexedEffectSet b -> [WithIntercept b]
effectSetMembers (IndexedEffectSet indexByEffect _) = M.keys indexByEffect

indexedEffectSetFromList :: Ord b => [WithIntercept b] -> IndexedEffectSet b
indexedEffectSetFromList es =
  let pairs         = zip [0 ..] $ L.sort es
      indexByEffect = M.fromList $ fmap (\(a, b) -> (b, a)) pairs
      effectByIndex = M.fromList pairs
  in  IndexedEffectSet indexByEffect effectByIndex

makeIndexedEffectSet
  :: (Ord b, Foldable f) => f (WithIntercept b) -> IndexedEffectSet b
makeIndexedEffectSet = indexedEffectSetFromList . FL.fold FL.list

addIndexedEffect
  :: Ord b => IndexedEffectSet b -> WithIntercept b -> IndexedEffectSet b
addIndexedEffect (IndexedEffectSet m _) x =
  indexedEffectSetFromList . (x :) . fmap fst $ M.toList m

effectIndex :: Ord b => IndexedEffectSet b -> WithIntercept b -> Maybe Int
effectIndex (IndexedEffectSet indexByEffect _) x = M.lookup x indexByEffect

effectAtIndex :: Ord b => IndexedEffectSet b -> Int -> Maybe (WithIntercept b)
effectAtIndex (IndexedEffectSet _ effectByIndex) n = M.lookup n effectByIndex
-}

-- serves as a holder for the intercept Choice and a proxy for the type b
data FixedEffects b where
  FixedEffects :: IndexedEffectSet b -> FixedEffects b
  InterceptOnly :: FixedEffects b
  deriving (Show)

allFixedEffects :: (Enum b, Ord b, Bounded b) => Bool -> FixedEffects b
allFixedEffects True = FixedEffects $ IS.fromList [minBound ..]
allFixedEffects False =
  FixedEffects $ IS.fromList $ fmap Predictor [minBound ..]

modelsIntercept :: (Ord b, Enum b, Bounded b) => FixedEffects b -> Bool
modelsIntercept (FixedEffects ef) =
  maybe False (const True) $ IS.index ef Intercept
modelsIntercept InterceptOnly = True

indexedFixedEffectSet
  :: (Ord b, Enum b, Bounded b)
  => FixedEffects b
  -> IS.IndexedSet (WithIntercept b)
indexedFixedEffectSet InterceptOnly    = IS.fromList [Intercept]
indexedFixedEffectSet (FixedEffects x) = x


data ItemInfo = ItemInfo { itemIndex :: Int, itemName :: T.Text } deriving (Show)
{-
g is the group and we need to map it to an Int, representing which group. E.g., "state" -> 0, "county" -> 1
The IndexedSet has our map to and from Int, i.e., the position in the vectors
The Vector of Vectors has our membership info where the vector for each row is indexed as in the indexed set
-}
data RowClassifier g where
  RowClassifier :: IS.IndexedSet g -> M.Map g Int -> VB.Vector (VB.Vector ItemInfo) -> RowClassifier g

instance Show g => Show (RowClassifier g) where
  show (RowClassifier indices sizes infos) = "RowClassifier " ++ show sizes ++ " " ++ show infos

groupSizes :: RowClassifier g -> M.Map g Int
groupSizes (RowClassifier _ sizes _) = sizes

groupIndices :: RowClassifier g -> IS.IndexedSet g
groupIndices (RowClassifier groupIndices _ _) = groupIndices

rowInfos :: RowClassifier g -> VB.Vector (VB.Vector ItemInfo)
rowInfos (RowClassifier _ _ infos) = infos


