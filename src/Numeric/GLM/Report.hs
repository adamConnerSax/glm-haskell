{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.GLM.Report where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.MixedModel        as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Numeric.GLM.Predict           as GLM
import qualified Numeric.SparseDenseConversions
                                               as SD

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
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


report
  :: (LA.Container LA.Vector Double, GLM.EffectsIO r)
  => GLM.MixedModel b g
  -> SLA.SpMatrix Double -- ^ smZ
  -> GLM.BetaVec -- ^ beta
  -> SLA.SpVector Double -- ^ b
  -> P.Sem r ()
report mm smZ vBeta svb = do
  let
    mX       = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
    vY       = GLM.rmsObservations $ GLM.regressionModelSpec mm
    groupFSM = GLM.mmsFitSpecByGroup $ GLM.mixedModelSpec mm
    (_, p)   = LA.size mX
    (_, q)   = SLA.dim smZ
    vb       = SD.toDenseVector svb
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
    groupReport l@(GLM.GroupFitSpec n b _) b' = do
      putStrLn $ show n ++ " groups"
      when b $ reportStats "Intercept" $ VS.take n b'
      let (numSlopes, bS) = case b of
            True  -> (GLM.effectsForGroup l - 1, VS.drop n b')
            False -> (GLM.effectsForGroup l, b')
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

colPredictions
  :: (Show b, GLM.GroupC g)
  => GLM.EffectParametersByGroup g b
  -> (Double -> T.Text)
  -> C.Colonnade
       C.Headed
       (T.Text, M.Map g (GLM.GroupKey g), Double)
       T.Text
colPredictions epgMap formatPrediction =
  let label (l, _, _) = l
      cats (_, c, _) = c
      pred (_, _, p) = p
      gLabel g = maybe "--" (T.pack . show) . M.lookup g
      makeGroupCol g = C.headed (T.pack $ show g) (gLabel g . cats)
  in  C.headed "Category" label
      <> mconcat (fmap makeGroupCol $ M.keys epgMap)
      <> C.headed "Prediction" (formatPrediction . pred)

printPredictions
  :: (Traversable f, GLM.PredictorC b, GLM.GroupC g, GLM.Effects r)
  => GLM.MixedModel b g
  -> GLM.FixedEffectParameters b
  -> GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> f (T.Text, M.Map b Double, M.Map g (GLM.GroupKey g))
  -> P.Sem r T.Text
printPredictions mm fes epg rowClassifier labeledToPredict =
  let f predMap gLabelMap =
        fst
          <$> GLM.predict mm
                          (flip M.lookup predMap)
                          (flip M.lookup gLabelMap)
                          fes
                          epg
                          rowClassifier
      makePrediction (l, predictors, groupLabels) = do
        prediction <- f predictors groupLabels
        return (l, groupLabels, prediction)
  in  fmap
        (T.pack . C.ascii (fmap T.unpack $ colPredictions epg (T.pack . show)))
        (traverse makePrediction labeledToPredict)


-- fitted values of input data
-- like "fitted" in lme4
fitted
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => GLM.MixedModel b g
  -> (q -> b -> Double)
  -> (q -> g -> GLM.GroupKey g)
  -> GLM.FixedEffectParameters b
  -> GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> q
  -> P.Sem r Double -- the map in effectParametersByGroup and the map in RowClassifier might not match
fitted mm getPred getLabel fep epg rowClassifier row =
  fst
    <$> GLM.predict mm
                    (Just . getPred row)
                    (Just . getLabel row)
                    fep
                    epg
                    rowClassifier


fitted'
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => GLM.MixedModel b g
  -> (q -> b -> Double)
  -> (q -> g -> GLM.GroupKey g)
  -> GLM.EffectsByGroup g b
  -> GLM.RowClassifier g
  -> GLM.BetaU
  -> LA.Vector Double
  -> q
  -> P.Sem r Double -- the map in effectParametersByGroup and the map in RowClassifier might not match
fitted' mm getPred getLabel ebg rowClassifier betaU vb row =
  fst
    <$> GLM.predict' mm
                     (Just . getPred row)
                     (Just . getLabel row)
                     ebg
                     rowClassifier
                     betaU
                     vb


-- reports for fixed Effects
-- first a table of intercept/predictor vs. coefficient and variance
colFixedEffects :: Show b => C.Colonnade C.Headed (b, Double, Double) T.Text
colFixedEffects =
  let eff (b, _, _) = b
      mean (_, m, _) = m
      var (_, _, v) = v
  in  C.headed "Effect" (T.pack . show . eff)
      <> C.headed "Parameter" (T.pack . show . mean)
      <> C.headed "Std. Dev" (T.pack . show . sqrt . var)

printFixedEffects
  :: (GLM.Effects r, GLM.PredictorC b)
  => GLM.FixedEffectStatistics b
  -> P.Sem r T.Text
printFixedEffects (GLM.FixedEffectStatistics fep covars) = do
  let
    GLM.FixedEffectParameters fes means = fep
    fromEff fe = GLM.eitherToSem $ do
      index <-
        maybe
            (  Left
            $  (T.pack $ show fe)
            <> " not found in "
            <> (T.pack $ show fes)
            )
            Right
          $ IS.index (GLM.indexedFixedEffectSet fes) fe
      return (fe, means `LA.atIndex` index, covars `LA.atIndex` (index, index))
    effList = IS.toList $ GLM.indexedFixedEffectSet fes
  rows <- traverse fromEff effList
  return $ T.pack $ C.ascii (fmap T.unpack $ colFixedEffects) rows

-- and a function to produce "predictions" from a given set of FE predictors
-- this so we can make a table reflecting the fixed-effects
{-
getIndexed :: Show b => GLM.IndexedEffectSet b -> (Int -> a) -> b -> Maybe a
getIndexed ies coeffF b = fmap coeffF (ies `index` b)

getIndexedE
  :: Show b => GLM.IndexedEffectSet b -> (Int -> a) -> b -> Either T.Text a
getIndexedE ies coeffF b = maybe (Left msg) Right (getIndexed ies coeffF b)
 where
  msg =
    "Failed to find \""
      <> (T.pack $ show b)
      <> "\" in index="
      <> (T.pack $ show ies)
-}

-- Like "ranef" in lme4
randomEffectsByLabel
  :: (GLM.Effects r, GLM.PredictorC b, GLM.GroupC g)
  => GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> P.Sem
       r
       (M.Map g (GLM.IndexedEffectSet b, [(T.Text, LA.Vector Double)]))
randomEffectsByLabel ebg rowClassifier = do
  let byLabel group (GLM.EffectParameters ies mEP) = do
        indexedKeys <- M.toList <$> GLM.labelMap rowClassifier group
        return
          $ ( ies
            , fmap (\(gk, r) -> (T.pack (show gk), LA.toRows mEP L.!! r))
                   indexedKeys
            )
  GLM.eitherToSem $ M.traverseWithKey byLabel ebg


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

