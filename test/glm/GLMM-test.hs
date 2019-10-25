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
import qualified Numeric.GLM.FunctionFamily    as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.GLM.Report            as GLM

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

import           System.IO                      ( hSetBuffering
                                                , stdout
                                                , BufferMode(..)
                                                )

verbose = False

throwEither :: (P.Member (P.Error GLMError) r) => Either GLMError a -> P.Sem r a
throwEither x = case x of
  Left  e -> P.throw e
  Right x -> return x

throwMaybe :: (P.Member (P.Error GLMError) r) => T.Text -> Maybe a -> P.Sem r a
throwMaybe msg x = throwEither $ maybe (Left $ OtherGLMError $ msg) Right x

main :: IO ()
main = do
  hSetBuffering stdout NoBuffering
  let lmmControls = LMMControls LMM_NELDERMEAD --defaultLMMControls
  let glmmControls =
        GLMMControls lmmControls GLM.UseCanonical 10 (ConvergeSimple 0.05 20)

{-
  frame <- defaultLoadToFrame @'[Rail, Travel] railCSV (const True)
  let getObservation = realToFrac . F.rgetField @Travel
      fixedEffects :: GLM.FixedEffects ()
      fixedEffects   = GLM.InterceptOnly
      groups         = IS.fromList [RG_Rail]
      getPredictor   = railPredictor
      groupLabels    = railGroupLabels
      effectsByGroup = M.fromList [(RG_Rail, IS.fromList [GLM.Intercept])]
      asLMM x = LMM x lmmControls
      vW = LA.fromList $ L.replicate (FL.fold FL.length frame) 1.0
      asGLMM x = GLMM x vW GLM.Normal glmmControls
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
    asLMM x = LMM x lmmControls
    vW = LA.fromList $ L.replicate (FL.fold FL.length frame) 1.0
    asGLMM x = GLMM x vW GLM.Normal glmmControls
-}
{-
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
      asLMM x = LMM x lmmControls
      vW = LA.fromList $ L.replicate (FL.fold FL.length frame) 1.0
      asGLMM x = GLMM x vW GLM.Normal glmmControls
-}

  frame <- defaultLoadToFrame @'[Row, Herd, Incidence, Size, Period, Obs]
    cbppCSV
    (const True)
  let getObservation r =
        (realToFrac $ F.rgetField @Incidence r)
          / (realToFrac $ F.rgetField @Size r)
      fixedEffects :: GLM.FixedEffects CbppPredictor
      fixedEffects   = GLM.InterceptOnly --GLM.allFixedEffects True -- model using CbppPredictor and with intercept
      groups         = IS.fromList [CbppHerd]
      getPredictor   = getCbppPredictor
      groupLabels    = cbppGroupLabels
      effectsByGroup = M.fromList [(CbppHerd, IS.fromList [GLM.Intercept])]
      vW             = LA.fromList $ L.replicate (FL.fold FL.length frame) 1.0
      vN = LA.fromList $ fmap (F.rgetField @Size) $ FL.fold FL.list frame
      asGLMM mm = GLMM mm vW (GLM.Binomial vN) glmmControls

  resultEither <- runEffectsIO $ do
    let (vY, mX, rcM) = FL.fold
          (lmePrepFrame getObservation
                        fixedEffects
                        groups
                        getPredictor
                        groupLabels
          )
          frame
    rowClassifier  <- throwEither $ either (Left . OtherGLMError) Right rcM
    fitSpecByGroup <- fitSpecByGroup fixedEffects effectsByGroup rowClassifier
    let (n, p) = LA.size mX
        rcRows = VB.length $ GLM.rowInfos rowClassifier
    when verbose $ liftIO $ do
      putStrLn $ "classifiers=" ++ show rowClassifier
    when (rcRows /= n)
      $  P.throw
      $  OtherGLMError
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
    smZ <- makeZ mX fitSpecByGroup rowClassifier
    let mkLambda         = makeLambda fitSpecByGroup
        (_, q)           = SLA.dim smZ
        mm = MixedModel (RegressionModel fixedEffects mX vY) fitSpecByGroup
        randomEffectCalc = RandomEffectCalculated smZ mkLambda
        th0              = setCovarianceVector fitSpecByGroup 1 0 -- LA.fromList [2, 2]
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ SD.toDenseMatrix smZ
    checkProblem mm randomEffectCalc
    let mdVerbosity = if verbose then MDVSimple else MDVNone
-- compare LMM and GLMM with ObservationDistribution set to Normal
{-    liftIO $ putStrLn "LMM"
    (th2_LMM, pd2_LMM, sigma2_LMM, vBetaU2_LMM, vb2_LMM, cs_LMM) <-
      minimizeDeviance mdVerbosity ML (asLMM mm) randomEffectCalc th0
-}
    liftIO $ putStrLn $ "GLMM"
    (th2_GLMM, pd2_GLMM, sigma2_GLMM, vBetaU2_GLMM, vb2_GLMM, cs_GLMM) <-
      minimizeDeviance mdVerbosity ML (asGLMM mm) randomEffectCalc th0
    liftIO $ do
      putStrLn $ "deviance=" ++ show pd2_GLMM
      putStrLn $ "beta=" ++ show (bu_vBeta vBetaU2_GLMM)
      putStrLn $ "u=" ++ show (bu_svU vBetaU2_GLMM)
      putStrLn $ "b=" ++ show vb2_GLMM
    GLM.report (asGLMM mm)
               smZ
               (bu_vBeta vBetaU2_GLMM)
               (SD.toSparseVector vb2_GLMM)
    let fes_GLMM = GLM.fixedEffectStatistics (asGLMM mm)
                                             sigma2_GLMM
                                             cs_GLMM
                                             vBetaU2_GLMM
    liftIO $ putStrLn $ "FixedEffectStatistics: " ++ show fes_GLMM
    epg <- GLM.effectParametersByGroup rowClassifier effectsByGroup vb2_GLMM
    liftIO $ putStrLn $ "EffectParametersByGroup: " ++ show epg
    gec <- GLM.effectCovariancesByGroup effectsByGroup
                                        (asGLMM mm)
                                        sigma2_GLMM
                                        th2_GLMM
    liftIO $ putStrLn $ "EffectCovariancesByGroup: " ++ show gec
    rebl <- GLM.randomEffectsByLabel epg rowClassifier
    liftIO
      $  putStrLn
      $  "Random Effects:\n"
      ++ (T.unpack $ GLM.printRandomEffectsByLabel rebl)
    let f r = do
          let obs = getObservation r
          fitted <- GLM.fitted (asGLMM mm)
                               getPredictor
                               groupLabels
                               fes_GLMM
                               epg
                               rowClassifier
                               r
          return (obs, fitted)
    fitted <- traverse f (FL.fold FL.list frame)
    liftIO $ putStrLn $ "Fitted:\n" ++ show fitted
    liftIO $ putStrLn "Done"

{-
    (th2_ML, pd2_ML, sigma2_ML, vBetaU2_ML, vb2_ML, cs_ML) <- minimizeDeviance
      mdVerbosity
      ML
      (glmm mm)
      randomEffectCalc
      th0
    liftIO $ do
      putStrLn $ "ML Via method 2"
      putStrLn $ "deviance=" ++ show pd2_ML
      putStrLn $ "beta=" ++ show (vBeta vBetaU2_ML)
      putStrLn $ "u=" ++ show (svU vBetaU2_ML)
      putStrLn $ "b=" ++ show vb2_ML
    GLM.report (glmm mm) smZ (vBeta vBetaU2_ML) (SD.toSparseVector vb2_ML)
    let fes_ML = GLM.fixedEffectStatistics (glmm mm) sigma2_ML cs_ML vBetaU2_ML
    liftIO $ putStrLn $ "FixedEffectStatistics: " ++ show fes_ML
    epg <- GLM.effectParametersByGroup rowClassifier effectsByGroup vb2_ML
    liftIO $ putStrLn $ "EffectParametersByGroup: " ++ show epg
    gec <- GLM.effectCovariancesByGroup effectsByGroup
                                        (glmm mm)
                                        sigma2_ML
                                        th2_ML
    liftIO $ putStrLn $ "EffectCovariancesByGroup: " ++ show gec
    rebl <- GLM.randomEffectsByLabel epg rowClassifier
    liftIO
      $  putStrLn
      $  "Random Effects:\n"
      ++ (T.unpack $ GLM.printRandomEffectsByLabel rebl)
    let f r = do
          let obs = getObservation r
          fitted <- GLM.fitted (glmm mm)
                               getPredictor
                               groupLabels
                               fes_ML
                               epg
                               rowClassifier
                               r
          return (obs, fitted)
    fitted <- traverse f (FL.fold FL.list frame)
    liftIO $ putStrLn $ "Fitted:\n" ++ show fitted
    when (not $ generalized $ glmm mm) $ do
      (th2_REML, pd2_REML, sigma2_REML, vBetaU2_REML, vb2_REML, cs_REML) <-
        minimizeDeviance mdVerbosity REML (glmm mm) randomEffectCalc th0
      liftIO $ do
        putStrLn $ "REML Via method 2"
        putStrLn $ "deviance=" ++ show pd2_REML
        putStrLn $ "sigma=" ++ show (sqrt sigma2_REML)
        putStrLn $ "beta=" ++ show (vBeta vBetaU2_REML)
        putStrLn $ "u=" ++ show (svU vBetaU2_REML)
        putStrLn $ "b=" ++ show vb2_REML
      GLM.report (glmm mm) smZ (vBeta vBetaU2_REML) (SD.toSparseVector vb2_REML)

    when verbose $ do
      cholmodFactor <- cholmodAnalyzeProblem randomEffectCalc
      (pdTest, sigma2Test, svBetaUTest, bTest, _) <- profiledDeviance
        PDVAll
        cholmodFactor
        ML
        (glmm mm)
        randomEffectCalc
        th2_ML
      liftIO $ do
        putStrLn $ "pdTest=" ++ show pdTest
        putStrLn $ "betaTest=" ++ show (vBeta $ svBetaUTest)
        putStrLn $ "uTest=" ++ show (SD.toDenseVector $ svU svBetaUTest)
        putStrLn $ "bTest=" ++ show (SD.toDenseVector $ bTest)
-}
  case resultEither of
    Left  e  -> putStrLn $ "Error: " ++ (show e)
    Right () -> putStrLn $ "Success!"



