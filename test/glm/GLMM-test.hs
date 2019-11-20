{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.GLM.Bootstrap         as GLM
import qualified Numeric.GLM.Predict           as GLM
import qualified Numeric.GLM.Report            as GLM

import qualified Numeric.SparseDenseConversions
                                               as SD
import           DataFrames
import qualified Frames                        as F

import qualified Statistics.Types              as S

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

runFIO =
  let runGLMEffects = if verbose then runEffectsVerboseIO else runEffectsIO
  in  runGLMEffects . GLM.asyncToIOFinal . GLM.runRandomIO

throwEither
  :: (P.Member (P.Error GLM.GLMError) r) => Either GLM.GLMError a -> P.Sem r a
throwEither x = case x of
  Left  e -> P.throw e
  Right x -> return x

throwMaybe
  :: (P.Member (P.Error GLM.GLMError) r) => T.Text -> Maybe a -> P.Sem r a
throwMaybe msg x = throwEither $ maybe (Left $ GLM.OtherGLMError $ msg) Right x

main :: IO ()
main = do
  hSetBuffering stdout NoBuffering
  let lmmControls = GLM.LMMControls GLM.LMM_BOBYQA 1e-6 --defaultLMMControls

  let glmmControls = GLM.GLMMControls
        GLM.UseCanonical
        10
        (GLM.PIRLSConvergenceCriterion GLM.PCT_Deviance 1e-9 30)

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
  let
    getObservation r =
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
    asGLMM lmms = GLM.GeneralizedLinearMixedModel
      (GLM.GeneralizedLinearMixedModelSpec lmms
                                           vW
                                           (GLM.Binomial vN)
                                           glmmControls
      )

  resultEither <- runFIO $ do
    let (vY, mX, rcM) = FL.fold
          (lmePrepFrame getObservation
                        fixedEffects
                        groups
                        getPredictor
                        groupLabels
          )
          frame
    rowClassifier  <- throwEither $ either (Left . GLM.OtherGLMError) Right rcM
    fitSpecByGroup <- GLM.fitSpecByGroup fixedEffects
                                         effectsByGroup
                                         rowClassifier
    let (n, p) = LA.size mX
        rcRows = VB.length $ GLM.rowInfos rowClassifier
    when verbose $ liftIO $ do
      putStrLn $ "classifiers=" ++ show rowClassifier
    when (rcRows /= n)
      $  P.throw
      $  GLM.OtherGLMError
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
    let
      mkLambda = makeLambda fitSpecByGroup
      (_, q)   = SLA.dim smZ
      lmms     = GLM.LinearMixedModelSpec
        (GLM.MixedModelSpec (GLM.RegressionModelSpec fixedEffects mX vY)
                            fitSpecByGroup
        )
        lmmControls
      randomEffectCalc = GLM.RandomEffectCalculated smZ mkLambda
      th0              = setCovarianceVector fitSpecByGroup 1 0 -- LA.fromList [2, 2]
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ SD.toDenseMatrix smZ
    let mdVerbosity = if verbose then MDVSimple else MDVNone
-- compare LMM and GLMM with ObservationDistribution set to Normal
{-
    liftIO $ putStrLn "LMM"
    ((th2_LMM, pd2_LMM, sigma2_LMM, vBetaU2_LMM, vb2_LMM, cs_LMM),_,_) <-
      minimizeDeviance mdVerbosity ML (asLMM mm) randomEffectCalc th0
    GLM.report (asLMM mm) smZ (bu_vBeta vBetaU2_LMM) (SD.toSparseVector vb2_LMM)
    let fes_LMM =
          GLM.fixedEffectStatistics (asLMM mm) sigma2_LMM cs_LMM vBetaU2_LMM
    liftIO $ putStrLn $ "FixedEffectStatistics: " ++ show fes_LMM
    epg_LMM <- GLM.effectParametersByGroup rowClassifier effectsByGroup vb2_LMM
    liftIO $ putStrLn $ "EffectParametersByGroup: " ++ show epg_LMM
    gec_LMM <- GLM.effectCovariancesByGroup effectsByGroup
                                            (asLMM mm)
                                            sigma2_LMM
                                            th2_LMM
    liftIO $ putStrLn $ "EffectCovariancesByGroup: " ++ show gec_LMM
    rebl_LMM <- GLM.randomEffectsByLabel epg_LMM rowClassifier
    liftIO
      $  putStrLn
      $  "Random Effects:\n"
      ++ (T.unpack $ GLM.printRandomEffectsByLabel rebl_LMM)
    let fLMM r = do
          let obs = getObservation r
          fitted <- GLM.fitted (asLMM mm)
                               getPredictor
                               groupLabels
                               fes_LMM
                               epg_LMM
                               rowClassifier
                               r
          return (obs, fitted)
    fitted_LMM <- traverse fLMM (FL.fold FL.list frame)
    liftIO $ putStrLn $ "Fitted:\n" ++ show fitted_LMM
-}
--
    liftIO $ putStrLn $ "GLMM"
    let mm = asGLMM lmms
    checkProblem mm randomEffectCalc
    ((th2_GLMM, pd2_GLMM, sigma2_GLMM, vBetaU2_GLMM, vb2_GLMM, cs_GLMM), vMuSol, cf) <-
      minimizeDeviance mdVerbosity ML mm randomEffectCalc th0
    liftIO $ do
      putStrLn $ "deviance=" ++ show pd2_GLMM
      putStrLn $ "beta=" ++ show (GLM.bu_vBeta vBetaU2_GLMM)
      putStrLn $ "u=" ++ show (GLM.bu_svU vBetaU2_GLMM)
      putStrLn $ "b=" ++ show vb2_GLMM
    GLM.report mm smZ (GLM.bu_vBeta vBetaU2_GLMM) (SD.toSparseVector vb2_GLMM)
    let fep_GLMM = GLM.fixedEffectParameters mm vBetaU2_GLMM
        fes_GLMM =
          GLM.fixedEffectStatistics mm sigma2_GLMM cs_GLMM vBetaU2_GLMM
    liftIO $ putStrLn $ "FixedEffectStatistics: " ++ show fes_GLMM
    epg <- GLM.effectParametersByGroup rowClassifier effectsByGroup vb2_GLMM
    liftIO $ putStrLn $ "EffectParametersByGroup: " ++ show epg
    gec <- GLM.effectCovariancesByGroup effectsByGroup mm sigma2_GLMM th2_GLMM
    liftIO $ putStrLn $ "EffectCovariancesByGroup: " ++ show gec
    rebl <- GLM.randomEffectsByLabel epg rowClassifier
    liftIO
      $  putStrLn
      $  "Random Effects:\n"
      ++ (T.unpack $ GLM.printRandomEffectsByLabel rebl)

    liftIO $ putStrLn $ "Boostrapping for confidence intervals"

    bootstraps <- GLM.parametricBootstrap mdVerbosity
                                          ML
                                          mm
                                          randomEffectCalc
                                          cf
                                          th2_GLMM
                                          vMuSol
                                          (sqrt sigma2_GLMM)
                                          10
                                          True
    let f r = do
          let obs = getObservation r
          fitted <- GLM.fitted mm
                               getPredictor
                               groupLabels
                               fep_GLMM
                               epg
                               rowClassifier
                               r
          fitted' <- GLM.fitted' mm
                                 getPredictor
                                 groupLabels
                                 fixedEffects
                                 effectsByGroup
                                 rowClassifier
                                 vBetaU2_GLMM
                                 vb2_GLMM
                                 r
{-                     
          bootWCI <- GLM.bootstrappedConfidence mm
                                                (Just . getPredictor r)
                                                (Just . groupLabels r)
                                                rowClassifier
                                                effectsByGroup
                                                (vBetaU2_GLMM, vb2_GLMM)
                                                bootstraps
                                                GLM.BCI_Accelerated
                                                (S.mkCL 0.95)
-}
          return (obs, fitted, fitted')

    fitted <- traverse f (FL.fold FL.list frame)
    liftIO $ putStrLn $ "Fitted:\n" ++ (L.intercalate "\n" $ fmap show fitted)
    liftIO $ putStrLn "Done"


  case resultEither of
    Left  e  -> putStrLn $ "Error: " ++ (show e)
    Right () -> putStrLn $ "Success!"




