{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.GLM.Report where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.Types             as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.SparseDenseConversions
                                               as SD

import qualified Polysemy                      as P
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Map                      as M
import           Data.Maybe                     ( catMaybes )
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.List                     as L
import qualified Colonnade                     as C

import qualified Numeric.LinearAlgebra.Class   as SLA

import qualified Numeric.LinearAlgebra         as LA


import qualified Data.Text                     as T
import qualified Data.Vector.Storable          as VS
import qualified Data.Vector.Split             as VS


{-
n rows
p fixed effects
l effect levels
n_l number of categories in level l
q_l number of (random) effects in level l
q = sum (n_l * q_l)
X is n x p
Z is n x q
zTz is q x q, this is the part that will be sparse.
-}

fixedEffectStatistics
  :: (Ord b, Enum b, Bounded b)
  => GLM.FixedEffects b
  -> Double
  -> CholeskySolutions
  -> SparseBetaU
  -> GLM.FixedEffectStatistics b
fixedEffectStatistics fe sigma2 (CholeskySolutions _ _ mRX) (BetaU svBeta _) =
  let means = SD.toDenseVector svBeta
      covs  = LA.scale sigma2 $ (LA.inv mRX) LA.<> (LA.inv $ LA.tr mRX)
  in  GLM.FixedEffectStatistics fe means covs

effectParametersByGroup
  :: (Ord g, Show g, Show b)
  => GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> LA.Vector Double
  -> Either T.Text (GLM.EffectParametersByGroup g b)
effectParametersByGroup rc ebg vb = do
  let groups = IS.members $ GLM.groupIndices rc
      groupOffset group = do
        size     <- GLM.groupSize rc group
        nEffects <- IS.size <$> GLM.groupEffects ebg group
        return $ size * nEffects
  offsets <- FL.prescan FL.sum <$> traverse groupOffset groups
  let
    ep (group, offset) = do
      size    <- GLM.groupSize rc group
      effects <- GLM.groupEffects ebg group
      let nEffects = IS.size effects
          parameterMatrix =
            LA.fromRows
              $ VS.chunksOf nEffects
              $ VS.take (size * nEffects)
              $ VS.drop offset vb
      return (group, GLM.EffectParameters effects parameterMatrix)
  fmap M.fromList $ traverse ep $ zip groups offsets


effectCovariancesByGroup
  :: (Ord g, Show g, Enum b, Bounded b, Show b)
  => GLM.EffectsByGroup g b
  -> Double
  -> LA.Vector Double
  -> Either T.Text (GLM.EffectCovariancesByGroup g b)
effectCovariancesByGroup ebg sigma2 vTh = do
  let
    groups = M.keys ebg
    nElements n = n * (n + 1) `div` 2
    groupOffset group = nElements . IS.size <$> GLM.groupEffects ebg group
    ltIndices n = [ (r, c) | c <- [0 .. (n - 1)], r <- [c .. (n - 1)] ]
    ltAssocs n v = zip (ltIndices n) (VS.toList v)
    tF ((r, c), x) = if r > c then Just ((c, r), x) else Nothing
    tAssocs l = catMaybes $ fmap tF l
    makeLT n v = let l = ltAssocs n v in LA.assoc (n, n) 0 $ l ++ (tAssocs l)
  offsets <- FL.prescan FL.sum <$> traverse groupOffset groups
  let ecs (group, offset) = do
        effects <- GLM.groupEffects ebg group
        let nEffects = IS.size effects
            cm =
              LA.scale sigma2
                $ makeLT nEffects
                $ VS.take (nElements nEffects)
                $ VS.drop offset vTh
        return (group, GLM.GroupEffectCovariances effects cm)
  fmap M.fromList $ traverse ecs $ zip groups offsets


report
  :: (LA.Container LA.Vector Double, SemC r)
  => Int -- ^ p
  -> Int -- ^ q
  -> FitSpecByGroup g
  -> LA.Vector Double -- ^ y
  -> LA.Matrix Double -- ^ X
  -> SLA.SpMatrix Double -- ^ Z
  -> SLA.SpVector Double -- ^ beta
  -> SLA.SpVector Double -- ^ b
  -> P.Sem r ()
report p q groupFSM vY mX smZ svBeta svb = do
  let
    vBeta = SD.toDenseVector svBeta
    vb    = SD.toDenseVector svb
    reportStats prefix v = do
      let mean = meanV v
          var =
            let v' = LA.cmap (\x -> x - mean) v
            in  (v' LA.<.> v') / (realToFrac $ LA.size v')
      putStrLn $ (T.unpack prefix) ++ ": mean (should be 0)=" ++ show mean
      putStrLn $ (T.unpack prefix) ++ ": variance=" ++ show var
      putStrLn $ (T.unpack prefix) ++ ": std. dev=" ++ show (sqrt var)
    vEps = vY - ((mX LA.#> vBeta) + SD.toDenseVector (smZ SLA.#> svb))
    meanV v = LA.sumElements v / realToFrac (LA.size v)
    mEps   = meanV vEps --LA.sumElements vEps / realToFrac (LA.size vEps)
    varEps = vEps LA.<.> vEps -- assumes mEps is 0.  Which it should be!!                
    groupReport l@(GroupFitSpec n b _) b' = do
      putStrLn $ show n ++ " groups"
      when b $ reportStats "Intercept" $ VS.take n b'
      let (numSlopes, bS) = case b of
            True  -> (effectsForGroup l - 1, VS.drop n b')
            False -> (effectsForGroup l, b')
      mapM_
        (\s -> reportStats ("Slope " <> (T.pack $ show s))
                           (VS.take n $ VS.drop (n * s) bS)
        )
        [0 .. (numSlopes - 1)]
  liftIO $ putStrLn $ "p=" ++ show p ++ "; q=" ++ show q
  liftIO $ reportStats "Residual" vEps
  let numberedGroups = zip [0 ..] (FL.fold FL.list groupFSM)
  liftIO $ mapM_
    (\(lN, l) -> putStrLn ("Level " ++ show lN) >> groupReport l vb)
    numberedGroups



-- fitted values of input data
-- like "fitted" in lme4
fitted
  :: (Ord g, Show g, Show b, Ord b, Enum b, Bounded b)
  => (r -> b -> Double)
  -> (r -> g -> T.Text)
  -> GLM.FixedEffectStatistics b
  -> GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> r
  -> Either T.Text Double -- the map in effectParametersByGroup and the map in RowClassifier might not match
fitted getPred getLabel (GLM.FixedEffectStatistics fe vFE _) ebg rowClassifier row
  = do
    let applyEffect b e r = case b of
          GLM.Intercept   -> e
          GLM.Predictor b -> e * getPred r b
        applyEffects ies v r =
          FL.fold FL.sum
            $ fmap
                (\(index, effect) -> applyEffect effect (v `LA.atIndex` index) r
                )
            $ IS.toIndexedList ies
        groups = M.keys ebg
        applyGroupEffects r group = do
          labelIndex <- GLM.labelIndex rowClassifier group (getLabel r group)
          (GLM.EffectParameters ies mEP) <- GLM.effectParameters group ebg
          return $ applyEffects ies (LA.toRows mEP L.!! labelIndex) r
        fixedEffects = applyEffects (GLM.indexedFixedEffectSet fe) vFE row
    groupEffects <- traverse (applyGroupEffects row) $ M.keys ebg
    return $ fixedEffects + FL.fold FL.sum groupEffects


-- Like "ranef" in lme4
randomEffectsByLabel
  :: (Ord g, Show g, Show b, Enum b, Bounded b)
  => GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> Either
       T.Text
       (M.Map g (GLM.IndexedEffectSet b, [(T.Text, LA.Vector Double)]))
randomEffectsByLabel ebg rowClassifier = do
  let byLabel group (GLM.EffectParameters ies mEP) = do
        indexedLabels <- M.toList <$> GLM.labelMap rowClassifier group
        return
          $ (ies, fmap (\(l, r) -> (l, LA.toRows mEP L.!! r)) indexedLabels)
  M.traverseWithKey byLabel ebg


colRandomEffectsByLabel
  :: Show b
  => GLM.IndexedEffectSet b
  -> C.Colonnade C.Headed (T.Text, LA.Vector Double) T.Text
colRandomEffectsByLabel ies =
  let
    indexedEffects = IS.toIndexedList ies
    colLabel       = C.headed "Category" fst
    colEffect (index, effect) = C.headed
      (T.pack $ show effect)
      (T.pack . show . flip LA.atIndex index . snd)
  in
    mconcat $ (colLabel : fmap colEffect indexedEffects)


printRandomEffectsByLabel
  :: (Show g, Show b)
  => M.Map g (GLM.IndexedEffectSet b, [(T.Text, LA.Vector Double)])
  -> T.Text
printRandomEffectsByLabel rebg =
  let printOne (g, (ies, x)) =
        (T.pack $ show g)
          <> ":\n"
          <> (T.pack $ C.ascii (fmap T.unpack $ colRandomEffectsByLabel ies) x)
          <> "\n"
  in  mconcat $ fmap printOne $ M.toList rebg
