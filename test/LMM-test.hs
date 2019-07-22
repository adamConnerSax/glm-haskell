{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import           Numeric.LinearMixedModel

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

{-
yV :: LA.Vector Double = LA.fromList $ L.take 20 $ L.iterate (+ 2) 0
--thV :: LA.Vector Double = LA.fromList [1, 1, 0.1, 2, 2, -0.1]

xM :: LA.Matrix Double = LA.matrix 2 $ L.take 40 $ L.iterate (+ 1) 0
levels :: VB.Vector (Int, Bool, Maybe (VB.Vector Bool)) =
  VB.fromList [(10, True, Just $ VB.fromList [False, True])
--  , ( 2
--    , True
--    , Nothing {-Just $ VB.fromList [True, False]-}
--    )
                                                           ]
rows =
  L.take 20 $ L.iterate (\[x, y] -> [(x + 1) `mod` 10, y + 1 `mod` 5]) [1, 0]
rowClassifier n = fmap VB.fromList rows !! n
-}



verbose = True

main :: IO ()
main = do

  railFrame <- defaultLoadToFrame @'[Rail, Travel] railCSV (const True)
  let (vY, mX, (vRC, numInCat)) = FL.fold
        (lmePrepFrameOne (realToFrac . F.rgetField @Travel)
                         (const $ LA.fromList [1])
                         (F.rgetField @Rail)
        )
        railFrame
      levels =
        VB.fromList [(numInCat, True, Nothing :: Maybe (VB.Vector Bool))]
      th0 = LA.fromList [2]

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
      levels = VB.fromList [(numInCat, True, Just $ VB.fromList [False, True])]
      th0    = LA.fromList [2, 2, 0.1]
-}
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
      levels :: VB.Vector (Int, Bool, Maybe (VB.Vector Bool)) =
        VB.fromList [(numInCat1, True, Nothing), (numInCat2, True, Nothing)]
      th0 = LA.fromList [2, 2]
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
    let (_, q) = SLA.dim smZ
    when verbose $ liftIO $ do
      putStrLn $ "Z="
      LA.disp 2 $ asDense smZ
    checkProblem mX vY smZ
    let smA = makeA mX vY smZ
    when verbose $ liftIO $ do
      putStrLn "A"
      LA.disp 2 $ asDense smA
    let makeST = makeSTF levels

    let (s, t) = makeST th0
    when verbose $ liftIO $ do
      putStrLn "S"
      LA.disp 2 $ asDense s
      putStrLn "T"
      LA.disp 2 $ asDense t
    let aStar = makeAStar p q smA makeST th0
    when verbose $ liftIO $ do
      putStrLn "A*"
      LA.disp 2 $ asDense aStar
    (th_ML, (pd_ML, ldL2_ML, r_ML, ldL_ML)) <- minimizeDeviance ML
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
    when verbose $ liftIO $ LA.disp 2 ldL_ML
    let (beta_ML, b_ML, bS_ML) =
          parametersFromSolution ML p q makeST th_ML ldL_ML r_ML
    liftIO $ do
      putStrLn $ "ML Fixed  (beta) =" ++ show beta_ML
      putStrLn $ "ML Random (th) =" ++ show th_ML
      putStrLn $ "ML Random (u, AKA b*) =" ++ show bS_ML
      putStrLn $ "ML Random (b) =" ++ show b_ML
    report p q levels vY mX smZ beta_ML b_ML
    (th_REML, (pd_REML, ldL2_REML, r_REML, ldL_REML)) <- minimizeDeviance
      REML
      p
      q
      n
      levels
      smA
      makeST
      th0
    liftIO $ putStrLn $ "REML Solution: profiled Deviance=" ++ (show pd_REML)
    liftIO $ putStrLn $ "REML Solution: r=" ++ show r_REML
    liftIO $ putStrLn $ "REML Solution: ldL2=" ++ show ldL2_REML
    when verbose $ liftIO $ LA.disp 2 ldL_REML
    let (beta_REML, b_REML, bS_REML) =
          parametersFromSolution REML p q makeST th_REML ldL_REML r_REML
    liftIO $ do
      putStrLn $ "REML Fixed  (beta) =" ++ show beta_REML
      putStrLn $ "REML Random (th) =" ++ show th_REML
      putStrLn $ "ML Random (u, AKA b*) =" ++ show bS_REML
      putStrLn $ "REML Random (b) =" ++ show b_REML
    report p q levels vY mX smZ beta_REML b_REML
  case resultEither of
    Left  err -> putStrLn $ "Error: " ++ (T.unpack err)
    Right ()  -> putStrLn $ "Success!"




