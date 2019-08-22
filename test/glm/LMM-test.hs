{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.Types             as GLM
import           Numeric.GLM.LinearMixedModel
import qualified Numeric.SparseDenseConversions
                                               as SD
import           DataFrames
import qualified Frames                        as F


import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )
import qualified Numeric.LinearAlgebra         as LA
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import qualified Data.Array                    as A
import qualified Data.Map                      as M
import qualified Data.Sparse.Common            as SLA
import qualified Data.List                     as L
import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB


verbose = False

throwEither :: (P.Member (P.Error T.Text) r) => Either T.Text a -> P.Sem r a
throwEither x = case x of
  Left  msg -> P.throw msg
  Right x   -> return x

throwMaybe :: (P.Member (P.Error T.Text) r) => T.Text -> Maybe a -> P.Sem r a
throwMaybe msg x = throwEither $ maybe (Left msg) Right x

main :: IO ()
main = do
{-
  frame <- defaultLoadToFrame @'[Rail, Travel] railCSV (const True)
  let getObservation = realToFrac . F.rgetField @Travel
      fixedEffects :: GLM.FixedEffects ()
      fixedEffects   = GLM.InterceptOnly
      groups         = IS.fromList [RG_Rail]
      getPredictor   = railPredictor
      groupLabels    = railGroupLabels
      effectsByGroup = M.fromList [(RG_Rail, IS.fromList [GLM.Intercept])]
-}
{-
  frame <- defaultLoadToFrame @'[Reaction, Days, Subject] sleepStudyCSV
                                                          (const True)
  let
    getObservation = realToFrac . F.rgetField @Reaction
    fixedEffects :: GLM.FixedEffects SleepStudyPredictor
    fixedEffects   = GLM.allFixedEffects True
    groups         = IS.fromList [SSG_Subject]
    getPredictor   = sleepStudyPredictor
    groupLabels    = sleepStudyGroupLabels
    effectsByGroup = M.fromList
      [(SSG_Subject, IS.fromList [GLM.Intercept, GLM.Predictor SleepStudyDays])]
-}

  frame <- defaultLoadToFrame @'[Block, Variety, Nitro, Yield] oatsCSV
                                                               (const True)
  let getObservation = realToFrac . F.rgetField @Yield
      fixedEffects :: GLM.FixedEffects OatsPredictor
      fixedEffects   = GLM.allFixedEffects True -- model using OatsPredictors and with intercept
      groups         = IS.fromList [OG_Block, OG_VarietyBlock]
      getPredictor   = getOatsPredictor
      groupLabels    = oatsGroupLabels
      effectsByGroup = M.fromList
        [ (OG_Block       , IS.fromList [GLM.Intercept])
        , (OG_VarietyBlock, IS.fromList [GLM.Intercept])
        ]

  resultEither <- runPIRLS_M $ do
    let (vY, mX, rcM) = FL.fold
          (lmePrepFrame getObservation
                        fixedEffects
                        groups
                        getPredictor
                        groupLabels
          )
          frame
    rowClassifier  <- throwEither rcM
    fitSpecByGroup <- throwEither
      $ fitSpecByGroup fixedEffects effectsByGroup rowClassifier
    let (n, p) = LA.size mX
        rcRows = VB.length $ GLM.rowInfos rowClassifier
    when verbose $ liftIO $ do
      putStrLn $ "classifiers=" ++ show rowClassifier
    when (rcRows /= n)
      $  P.throw
      $  "Only "
      <> (T.pack $ show rcRows)
      <> " in vRC but there are "
      <> (T.pack $ show n)
      <> " rows in the data!"
    --let rowClassifier n = vRC VB.! n
    when verbose $ liftIO $ do
      putStrLn $ show $ fmap colsForGroup fitSpecByGroup
      putStrLn $ "y=" ++ show vY
      putStrLn "X="
      LA.disp 2 mX
      putStrLn $ "levels=" ++ show fitSpecByGroup
    smZ <- throwMaybe "Error making Z, the random effect model matrix"
      $ makeZ mX fitSpecByGroup rowClassifier
    let
      mkLambda = makeLambda fitSpecByGroup
      (_, q)   = SLA.dim smZ
      mixedModel =
        MixedModel (RegressionModel fixedEffects mX vY) fitSpecByGroup
      randomEffectCalc = RandomEffectCalculated smZ mkLambda
      th0              = setCovarianceVector fitSpecByGroup 1 0 -- LA.fromList [2, 2]
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ SD.toDenseMatrix smZ
    checkProblem mixedModel randomEffectCalc
    let mdVerbosity = if verbose then MDVSimple else MDVNone
    (th2_ML, pd2_ML, sigma2_ML, vBeta2_ML, vu2_ML, vb2_ML, cs_ML) <-
      minimizeDeviance mdVerbosity ML mixedModel randomEffectCalc th0
    liftIO $ do
      putStrLn $ "ML Via method 2"
      putStrLn $ "deviance=" ++ show pd2_ML
      putStrLn $ "beta=" ++ show vBeta2_ML
      putStrLn $ "u=" ++ show vu2_ML
      putStrLn $ "b=" ++ show vb2_ML
    report p
           q
           fitSpecByGroup
           vY
           mX
           smZ
           (SD.toSparseVector vBeta2_ML)
           (SD.toSparseVector vb2_ML)
    (th2_REML, pd2_REML, sigma2_REML, vBeta2_REML, vu2_REML, vb2_REML, cs_REML) <-
      minimizeDeviance mdVerbosity REML mixedModel randomEffectCalc th0
    liftIO $ do
      putStrLn $ "REML Via method 2"
      putStrLn $ "deviance=" ++ show pd2_REML
      putStrLn $ "sigma=" ++ show (sqrt sigma2_REML)
      putStrLn $ "beta=" ++ show vBeta2_REML
      putStrLn $ "u=" ++ show vu2_REML
      putStrLn $ "b=" ++ show vb2_REML
    report p
           q
           fitSpecByGroup
           vY
           mX
           smZ
           (SD.toSparseVector vBeta2_REML)
           (SD.toSparseVector vb2_REML)
    let fes_REML = fixedEffectStatistics fixedEffects sigma2_REML cs_REML
    liftIO $ putStrLn $ "FixedEffectStatistics: " ++ show fes_REML
    epg <- throwEither
      $ effectParametersByGroup rowClassifier effectsByGroup vb2_REML
    liftIO $ putStrLn $ "EffectParametersByGroup: " ++ show epg
    gec <- throwEither
      $ effectCovariancesByGroup effectsByGroup sigma2_REML th2_REML
    liftIO $ putStrLn $ "EffectCovariancesByGroup: " ++ show gec
    rebl <- throwEither $ randomEffectsByLabel epg rowClassifier
    liftIO $ putStrLn $ "Random Effects: " ++ show rebl
    let f r = do
          let obs = getObservation r
          fitted <- fitted getPredictor groupLabels fes_REML epg rowClassifier r
          return (obs, fitted)
    fitted <- throwEither $ traverse f (FL.fold FL.list frame)
    liftIO $ putStrLn $ "Fitted:\n" ++ show fitted
    when verbose $ do
      cholmodFactor <- cholmodAnalyzeProblem randomEffectCalc
      (pdTest, sigma2Test, betaTest, uTest, bTest, _) <-
        liftIO $ profiledDeviance PDVAll
                                  cholmodFactor
                                  REML
                                  mixedModel
                                  randomEffectCalc
                                  th2_REML
      liftIO $ do
        putStrLn $ "pdTest=" ++ show pdTest
        putStrLn $ "betaTest=" ++ show betaTest
        putStrLn $ "uTest=" ++ show (SD.toDenseVector uTest)
        putStrLn $ "bTest=" ++ show (SD.toDenseVector bTest)
  case resultEither of
    Left  err -> putStrLn $ "Error: " ++ (T.unpack err)
    Right ()  -> putStrLn $ "Success!"




