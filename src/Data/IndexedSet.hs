{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Data.IndexedSet
  (
    -- * Types
    IndexedSet
    -- * Constuction
  , empty
  , fromList
  , fromFoldable
  , add
    -- * Conversion
  , size
  , toList
  , toSet
  , members
    -- * Access
  , index
  , atIndex
  -- * Subset
  , subset
  -- * Map/Traverse
  , mapIS
  , traverseIS
  )
where

import qualified Control.Foldl                 as FL
import qualified Data.List.Extra               as LE
import qualified Data.IntMap                   as IM
import qualified Data.Map                      as M
import qualified Data.Set                      as S
import qualified Data.Text                     as T


-- don't export this constructor !!
data IndexedSet a where
  IndexedSet :: M.Map a Int -> IM.IntMap a -> IndexedSet a
    deriving (Show)

-- construction
empty :: IndexedSet a
empty = IndexedSet M.empty IM.empty

fromList :: Ord a => [a] -> IndexedSet a
fromList as =
  let swap (a, b) = (b, a)
      pairs   = zip [0 ..] $ LE.nubSort as
      indexBy = M.fromAscList $ fmap swap pairs
      byIndex = IM.fromAscList pairs
  in  IndexedSet indexBy byIndex

fromFoldable :: (Ord a, Foldable f) => f a -> IndexedSet a
fromFoldable = fromList . FL.fold FL.list

add :: Ord a => IndexedSet a -> a -> IndexedSet a
add is x = fromList . (x :) $ members is

-- use as index
index :: Ord a => IndexedSet a -> a -> Maybe Int
index (IndexedSet indexBy _) x = M.lookup x indexBy

atIndex :: Ord a => IndexedSet a -> Int -> Maybe a
atIndex (IndexedSet _ byIndex) n = IM.lookup n byIndex

-- conversion
members :: IndexedSet a -> [a]
members = toList

size :: IndexedSet a -> Int
size (IndexedSet indexBy _) = M.size indexBy

toList :: IndexedSet a -> [a]
toList (IndexedSet indexBy _) = M.keys indexBy

toSet :: Ord a => IndexedSet a -> S.Set a
toSet = S.fromList . toList

-- mapping/traversing
mapIS :: Ord b => (a -> b) -> IndexedSet a -> IndexedSet b
mapIS f = fromList . fmap f . toList

traverseIS
  :: (Ord b, Applicative g) => (a -> g b) -> IndexedSet a -> g (IndexedSet b)
traverseIS f = fmap fromList . traverse f . toList

-- subset
subset :: Ord a => IndexedSet a -> IndexedSet a -> Bool
subset (IndexedSet sub _) (IndexedSet super _) =
  let f _ _ = True in M.isSubmapOfBy f sub super


