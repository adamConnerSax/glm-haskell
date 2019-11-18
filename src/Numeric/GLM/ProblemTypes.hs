{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Numeric.GLM.ProblemTypes
  ( WithIntercept(..)
  , IndexedEffectSet
  , FixedEffects(..)
  , indexedFixedEffectSet
  , allFixedEffects
  , ItemInfo(..)
  , RowClassifier(..)
  , groupIndices
  , groupSizes
  , groupSize
  , rowInfos
  , labelMap
  , labelIndex
  , EffectsByGroup
  , groupEffects
  , categoryNumberFromRowIndex
  , categoryNumberFromLabel
  , FixedEffectParameters(..)
  , FixedEffectStatistics(..)
  , GroupEffectCovariances(..)
  , EffectCovariancesByGroup
  , EffectParameters(..)
  , EffectParametersByGroup
  , effectParameters
  )
where

import qualified Data.IndexedSet               as IS

import qualified Control.Foldl                 as FL
import qualified Data.Array                    as A
import qualified Data.List                     as L
import qualified Data.Map                      as M
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Numeric.LinearAlgebra         as LA

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

eitherLookup :: (Show k, Ord k, Show x) => k -> M.Map k x -> Either T.Text x
eitherLookup k m =
  maybe
      (Left $ "\"" <> (T.pack $ show k) <> " not found in " <> (T.pack $ show m)
      )
      Right
    $ M.lookup k m

-- serves as a holder for the intercept Choice and a proxy for the type b
data FixedEffects b where
  FixedEffects :: IndexedEffectSet b -> FixedEffects b
  InterceptOnly :: FixedEffects b
  deriving (Show, Eq)

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

type EffectsByGroup g b = M.Map g (IndexedEffectSet b)

groupEffects
  :: (Show g, Ord g, Show b)
  => EffectsByGroup g b
  -> g
  -> Either T.Text (IndexedEffectSet b)
groupEffects ebg group = eitherLookup group ebg

groupEffectIndex
  :: (Ord g, Ord b) => EffectsByGroup g b -> g -> (WithIntercept b) -> Maybe Int
groupEffectIndex ebg group effect = M.lookup group ebg >>= flip IS.index effect

data ItemInfo = ItemInfo { itemIndex :: Int, itemName :: T.Text } deriving (Show)

{-
g is the group and we need to map it to an Int, representing which group. E.g., "state" -> 0, "county" -> 1
The IndexedSet has our map to and from Int, i.e., the position in the vectors
The Vector of Vectors has our membership info where the vector for each row is indexed as in the indexed set
Finally, the @Map g (Map Text Int)@ contains the info required to map
a specific category in a group to its category index  
-}
data RowClassifier g where
  RowClassifier :: IS.IndexedSet g -> M.Map g Int -> VB.Vector (VB.Vector ItemInfo) -> M.Map g (M.Map T.Text Int) -> RowClassifier g

-- get the row category for a given row and group
categoryNumberFromRowIndex :: Ord g => RowClassifier g -> Int -> g -> Maybe Int
categoryNumberFromRowIndex (RowClassifier groupIndices _ rowClassifierV _) rowIndex group
  = do
    vectorIndex <- groupIndices `IS.index` group
    return $ itemIndex $ (rowClassifierV VB.! rowIndex) VB.! vectorIndex --rowIndex VB.! levelIndex  

categoryNumberFromLabel :: Ord g => RowClassifier g -> T.Text -> g -> Maybe Int
categoryNumberFromLabel (RowClassifier _ _ _ categoryIndexByGroup) label group
  = M.lookup group categoryIndexByGroup >>= M.lookup label

instance Show g => Show (RowClassifier g) where
  show (RowClassifier indices sizes infos labelMaps) = "RowClassifier " ++ show indices ++ " " ++ show sizes ++ " " ++ show infos ++ " " ++ show labelMaps

groupSizes :: RowClassifier g -> M.Map g Int
groupSizes (RowClassifier _ sizes _ _) = sizes

groupSize :: (Show g, Ord g) => RowClassifier g -> g -> Either T.Text Int
groupSize rc group = eitherLookup group (groupSizes rc)

groupIndices :: RowClassifier g -> IS.IndexedSet g
groupIndices (RowClassifier groupIndices _ _ _) = groupIndices

rowInfos :: RowClassifier g -> VB.Vector (VB.Vector ItemInfo)
rowInfos (RowClassifier _ _ infos _) = infos

labelMaps :: RowClassifier g -> M.Map g (M.Map T.Text Int)
labelMaps (RowClassifier _ _ _ labelMaps) = labelMaps

labelMap
  :: (Show g, Ord g) => RowClassifier g -> g -> Either T.Text (M.Map T.Text Int)
labelMap rc group = eitherLookup group $ labelMaps rc

labelIndex
  :: (Show g, Ord g) => RowClassifier g -> g -> T.Text -> Either T.Text Int
labelIndex rc group label = do
  labelMap <- eitherLookup group (labelMaps rc)
  eitherLookup label labelMap

-- types for storing results
-- we fill this in at the end with the (mean) estimates and covariances
data FixedEffectParameters b where
  FixedEffectParameters :: FixedEffects b -> LA.Vector Double -> FixedEffectParameters b
  deriving (Show)

data FixedEffectStatistics b where
  FixedEffectStatistics :: FixedEffectParameters b -> LA.Matrix Double -> FixedEffectStatistics b
  deriving (Show)

data GroupEffectCovariances b where
  GroupEffectCovariances :: IndexedEffectSet b -> LA.Matrix Double -> GroupEffectCovariances b
  deriving (Show)

type EffectCovariancesByGroup g b = M.Map g (GroupEffectCovariances b)

-- One column per effect, estimates for each category. 
data EffectParameters b where
  EffectParameters :: IndexedEffectSet b -> LA.Matrix Double -> EffectParameters b
  deriving (Show)

type EffectParametersByGroup g b = M.Map g (EffectParameters b)

effectParameters
  :: (Show g, Ord g, Show b)
  => g
  -> EffectParametersByGroup g b
  -> Either T.Text (EffectParameters b)
effectParameters = eitherLookup

