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
  railFrame <- defaultLoadToFrame @'[Rail, Travel] railCSV (const True)
  let railFixedEffects :: GLM.FixedEffects ()
      railFixedEffects = GLM.InterceptOnly
      railGroups       = IS.fromList [RG_Rail]
      (vY, mX, rcM)    = FL.fold
        (lmePrepFrame (realToFrac . F.rgetField @Travel)
                      railFixedEffects
                      railGroups
                      railPredictor
                      railGroupLabels
        )
        railFrame
      groupEffectMap = M.fromList [(RG_Rail, IS.fromList [GLM.Intercept])]
      fixedEffects   = railFixedEffects
-}
{-
  sleepStudyFrame <- defaultLoadToFrame @'[Reaction, Days, Subject]
    sleepStudyCSV
    (const True)
  let
    sleepStudyFixedEffects :: GLM.FixedEffects SleepStudyPredictor
    sleepStudyFixedEffects = GLM.allFixedEffects True
    sleepStudyGroups = IS.fromList [SSG_Subject]
    (vY, mX, rcM)          = FL.fold
      (lmePrepFrame (realToFrac . F.rgetField @Reaction)
                    sleepStudyFixedEffects
                    sleepStudyGroups
                    sleepStudyPredictor
                    sleepStudyGroupLabels
      )
      sleepStudyFrame
    groupEffectMap = M.fromList
      [(SSG_Subject, IS.fromList [GLM.Intercept, GLM.Predictor SleepStudyDays])]
    fixedEffects = sleepStudyFixedEffects
-}

  oatsFrame <- defaultLoadToFrame @'[Block, Variety, Nitro, Yield]
    oatsCSV
    (const True)
  let oatsFixedEffects :: GLM.FixedEffects OatsPredictor
      oatsFixedEffects = GLM.allFixedEffects True -- model using OatsPredictors and with intercept
      oatsGroups       = IS.fromList [OG_Block, OG_VarietyBlock]
      (vY, mX, rcM)    = FL.fold
        (lmePrepFrame (realToFrac . F.rgetField @Yield)
                      oatsFixedEffects
                      oatsGroups
                      getOatsPredictor
                      oatsGroupLabels
        )
        oatsFrame
      groupEffectMap = M.fromList
        [ (OG_Block       , IS.fromList [GLM.Intercept])
        , (OG_VarietyBlock, IS.fromList [GLM.Intercept])
        ]
      fixedEffects = oatsFixedEffects

  resultEither <- runPIRLS_M $ do
    rowClassifier <- throwEither rcM
    groupInfoList <- throwMaybe "missing group in groupEffectMap" $ do
      traverse
          (\(grp, ie) ->
            M.lookup grp (GLM.groupSizes rowClassifier) >>= return . (grp, , ie)
          )
        $ M.toList groupEffectMap
    groupFSM <- throwEither $ fmap M.fromList $ traverse
      (\(grp, n, ige) ->
        makeGroupFitSpec n fixedEffects ige >>= return . (grp, )
      )
      groupInfoList
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
      putStrLn $ show $ fmap colsForGroup groupFSM
      putStrLn $ "y=" ++ show vY
      putStrLn "X="
      LA.disp 2 mX
      putStrLn $ "levels=" ++ show groupFSM
    smZ <- throwMaybe "Error making Z, the random effect model matrix"
      $ makeZ mX groupFSM rowClassifier
    let mkLambda         = makeLambda groupFSM
        (_, q)           = SLA.dim smZ
        mixedModel = MixedModel (RegressionModel fixedEffects mX vY) groupFSM
        randomEffectCalc = RandomEffectCalculated smZ mkLambda
        th0              = setCovarianceVector groupFSM 1 0 -- LA.fromList [2, 2]
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ SD.toDenseMatrix smZ
    checkProblem mixedModel randomEffectCalc
    let mdVerbosity = if verbose then MDVSimple else MDVNone
    (th2_ML, pd2_ML, vBeta2_ML, vu2_ML, vb2_ML) <- minimizeDeviance
      mdVerbosity
      ML
      mixedModel
      randomEffectCalc
      th0
    liftIO $ do
      putStrLn $ "ML Via method 2"
      putStrLn $ "deviance=" ++ show pd2_ML
      putStrLn $ "beta=" ++ show vBeta2_ML
      putStrLn $ "u=" ++ show vu2_ML
      putStrLn $ "b=" ++ show vb2_ML
    report p
           q
           groupFSM
           vY
           mX
           smZ
           (SD.toSparseVector vBeta2_ML)
           (SD.toSparseVector vb2_ML)
    (th2_REML, pd2_REML, vBeta2_REML, vu2_REML, vb2_REML) <- minimizeDeviance
      mdVerbosity
      REML
      mixedModel
      randomEffectCalc
      th0
    liftIO $ do
      putStrLn $ "REML Via method 2"
      putStrLn $ "deviance=" ++ show pd2_REML
      putStrLn $ "beta=" ++ show vBeta2_REML
      putStrLn $ "u=" ++ show vu2_REML
      putStrLn $ "b=" ++ show vb2_REML
    report p
           q
           groupFSM
           vY
           mX
           smZ
           (SD.toSparseVector vBeta2_REML)
           (SD.toSparseVector vb2_REML)
    when verbose $ do
      cholmodFactor                    <- cholmodAnalyzeProblem randomEffectCalc
      (pdTest, betaTest, uTest, bTest) <- liftIO $ profiledDeviance
        PDVAll
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




