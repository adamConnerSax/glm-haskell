{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import           Numeric.LinearMixedModel
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
import qualified Data.Sparse.Common            as SLA
import qualified Data.List                     as L
import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB


verbose = True

main :: IO ()
main = do
{-
  railFrame <- defaultLoadToFrame @'[Rail, Travel] railCSV (const True)
  let (vY, mX, (vRC, numInCat)) = FL.fold
        (lmePrepFrameOne (realToFrac . F.rgetField @Travel)
                         (const $ LA.fromList [1])
                         (F.rgetField @Rail)
        )
        railFrame
      levels = VB.fromList [LevelSpec numInCat True Nothing]
      th0    = setCovarianceVector levels 1 0 -- LA.fromList [2]
-}
{-
  sleepStudyFrame <- defaultLoadToFrame @'[Reaction, Days, Subject]
    sleepStudyCSV
    (const True)
  let (vY, mX, (vRC, numInCat)) = FL.fold
        (lmePrepFrameOne
          (realToFrac . F.rgetField @Reaction)
          ((\x -> LA.fromList [1, realToFrac x]) . F.rgetField @Days)
          (F.rgetField @Subject)
        )
        sleepStudyFrame
      levels = VB.fromList
        [LevelSpec numInCat True (Just $ VB.fromList [False, True])]
      th0 = setCovarianceVector levels 1 0 --LA.fromList [1, 1, 0]
-}

  oatsFrame <- defaultLoadToFrame @'[Block, Variety, Nitro, Yield]
    oatsCSV
    (const True)
  let (vY, mX, (vRC, numInCat1, numInCat2)) = FL.fold
        (lmePrepFrameTwo
          (realToFrac . F.rgetField @Yield)
          ((\x -> LA.fromList [1, realToFrac x]) . F.rgetField @Nitro)
          (\r -> (F.rgetField @Variety r, F.rgetField @Block r))
          (F.rgetField @Block)
        )
        oatsFrame
      levels = VB.fromList
        [LevelSpec numInCat1 True Nothing, LevelSpec numInCat2 True Nothing]
      th0 = setCovarianceVector levels 1 0 -- LA.fromList [2, 2]

  resultEither <- runPIRLS_M $ do
    let (n, p) = LA.size mX
        rcRows = VB.length vRC
    when verbose $ liftIO $ do
      putStrLn $ "classifiers=" ++ show vRC
    when (rcRows /= n)
      $  P.throw
      $  "Only "
      <> (T.pack $ show rcRows)
      <> " in vRC but there are "
      <> (T.pack $ show n)
      <> " rows in the data!"
    let rowClassifier n = vRC VB.! n
    when verbose $ liftIO $ do
      putStrLn $ show $ fmap colsForLevel levels
      putStrLn $ "y=" ++ show vY
      putStrLn "X="
      LA.disp 2 mX
      putStrLn $ "levels=" ++ show levels
    let smZ              = makeZ mX levels rowClassifier
        makeST           = makeSTF levels
        mkLambda         = makeLambda levels
        (_, q)           = SLA.dim smZ
        mixedModel       = MixedModel (RegressionModel mX vY) levels
        randomEffectCalc = RandomEffectCalculated smZ mkLambda
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ SD.toDenseMatrix smZ
    checkProblem mixedModel randomEffectCalc
    (th2_ML, pd2_ML, vBeta2_ML, vu2_ML, vb2_ML) <- minimizeDeviance2
      MDVSimple
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
           levels
           vY
           mX
           smZ
           (SD.toSparseVector vBeta2_ML)
           (SD.toSparseVector vb2_ML)
    (th2_REML, pd2_REML, vBeta2_REML, vu2_REML, vb2_REML) <- minimizeDeviance2
      MDVSimple
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
           levels
           vY
           mX
           smZ
           (SD.toSparseVector vBeta2_REML)
           (SD.toSparseVector vb2_REML)
{-           
    cholmodFactor                    <- cholmodAnalyzeProblem randomEffectCalc
    (pdTest, betaTest, uTest, bTest) <- liftIO $ profiledDeviance2
      PDVAll
      cholmodFactor
      REML
      mixedModel
      randomEffectCalc
      (LA.fromList [0.855, 1.127])
    liftIO $ do
      putStrLn $ "pdTest=" ++ show pdTest
      putStrLn $ "betaTest=" ++ show betaTest
      putStrLn $ "uTest=" ++ show uTest
      putStrLn $ "bTest=" ++ show bTest
-}
  case resultEither of
    Left  err -> putStrLn $ "Error: " ++ (T.unpack err)
    Right ()  -> putStrLn $ "Success!"




