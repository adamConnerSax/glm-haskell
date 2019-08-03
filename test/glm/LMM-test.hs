{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import           Numeric.LinearMixedModel
import qualified Numeric.SparseDenseConversions as SD
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
      levels =
        VB.fromList [LevelSpec numInCat True Nothing]
      th0 = setCovarianceVector levels 1 0 -- LA.fromList [2]
-}

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
      levels = VB.fromList [LevelSpec numInCat True (Just $ VB.fromList [False, True])]
      th0    = setCovarianceVector levels 1 0 --LA.fromList [1, 1, 0]

{-
  oatsFrame <- defaultLoadToFrame @'[Block, Variety, Nitro, Yield]
    oatsCSV
    (const True)
  let (vY, mX, (vRC, numInCat1, numInCat2)) = FL.fold
        (lmePrepFrameTwo
          (realToFrac . F.rgetField @Yield)
          ((\x -> LA.fromList [1, realToFrac x]) . F.rgetField @Nitro)
          (F.rgetField @Block)
          (\r -> (F.rgetField @Variety r, F.rgetField @Block r))
        )
        oatsFrame
      levels =
        VB.fromList [LevelSpec numInCat1 True Nothing, LevelSpec numInCat2 True Nothing)]
      th0 = setCovarianceVector levels 1 0 -- LA.fromList [2, 2]
-}
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
    let smZ    = makeZ mX levels rowClassifier
        makeST = makeSTF levels
        mkLambda = makeLambda levels
        (_, q) = SLA.dim smZ
        mixedModel = MixedModel (RegressionModel mX vY) levels
        randomEffectCalc = RandomEffectCalculated smZ mkLambda
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ SD.toDenseMatrix smZ
    checkProblem mixedModel randomEffectCalc
    let smA = makeA mixedModel randomEffectCalc
    when verbose $ liftIO $ do
      putStrLn "A"
      LA.disp 2 $ SD.toDenseMatrix smA
    let (s, t) = makeST th0
    when verbose $ liftIO $ do
      putStrLn "S"
      LA.disp 2 $ SD.toDenseMatrix s
      putStrLn "T"
      LA.disp 2 $ SD.toDenseMatrix t
    let aStar = makeAStar p q smA makeST th0
    when verbose $ liftIO $ do
      putStrLn "A*"
      LA.disp 2 $ SD.toDenseMatrix aStar
    (th_ML, perm_ML, (pd_ML, ldL2_ML, r_ML, ldL_ML)) <- minimizeDeviance
      ML
      p
      q
      n
      levels
      smA
      makeST
      th0
    liftIO $ putStrLn $ "ML Solution: profiled Deviance=" ++ (show pd_ML)
    liftIO $ putStrLn $ "ML Solution: r=" ++ show r_ML
    liftIO $ putStrLn $ "ML Solution: ldL2=" ++ show ldL2_ML
    when verbose $ liftIO $ do
      putStrLn "perm="
      LA.disp 0 $ SD.toDenseMatrix perm_ML
      putStrLn "L="
      LA.disp 2 $ SD.toDenseMatrix $ ldL_ML
    (beta_ML, b_ML, bS_ML) <- parametersFromSolution ML
                                                     p
                                                     q
                                                     makeST
                                                     th_ML
                                                     (perm_ML, ldL_ML)
                                                     r_ML
    liftIO $ do
      putStrLn $ "ML Fixed  (beta) =" ++ show (SD.toDenseVector beta_ML)
      putStrLn $ "ML Random (th) =" ++ show th_ML
      putStrLn $ "ML Random (u, AKA b*) =" ++ show (SD.toDenseVector bS_ML)
      putStrLn $ "ML Random (b) =" ++ show (SD.toDenseVector b_ML)
    report p q levels vY mX smZ beta_ML b_ML

    (th2_ML, pd2_ML, vBeta2_ML, vu2_ML, vb2_ML) <- minimizeDeviance2 ML mixedModel randomEffectCalc th0
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

    (th_REML, perm_REML, (pd_REML, ldL2_REML, r_REML, ldL_REML)) <-
      minimizeDeviance REML p q n levels smA makeST th0
    liftIO $ putStrLn $ "REML Solution: profiled Deviance=" ++ (show pd_REML)
    liftIO $ putStrLn $ "REML Solution: r=" ++ show r_REML
    liftIO $ putStrLn $ "REML Solution: ldL2=" ++ show ldL2_REML
    when verbose $ liftIO $ do
      putStrLn "perm="
      LA.disp 0 $ SD.toDenseMatrix perm_REML
      putStrLn "L="
      LA.disp 2 $ SD.toDenseMatrix $ ldL_REML
    (beta_REML, b_REML, bS_REML) <- parametersFromSolution
      REML
      p
      q
      makeST
      th_REML
      (perm_REML, ldL_REML)
      r_REML
    liftIO $ do
      putStrLn $ "REML Fixed  (beta) =" ++ show (SD.toDenseVector beta_REML)
      putStrLn $ "REML Random (th) =" ++ show th_REML
      putStrLn $ "ML Random (u, AKA b*) =" ++ show (SD.toDenseVector bS_REML)
      putStrLn $ "REML Random (b) =" ++ show (SD.toDenseVector b_REML)
    report p q levels vY mX smZ beta_REML b_REML

    (th2_REML, pd2_REML, vBeta2_REML, vu2_REML, vb2_REML) <- minimizeDeviance2 REML mixedModel randomEffectCalc th0
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

  case resultEither of
    Left  err -> putStrLn $ "Error: " ++ (T.unpack err)
    Right ()  -> putStrLn $ "Success!"




