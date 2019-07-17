{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Main where

import           Numeric.PIRLS

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

thV :: LA.Vector Double = LA.fromList [2]



main :: IO ()
main = do
  railFrame <- defaultLoadToFrame @'[Rail, Travel] railCSV (const True)
  let (vY, mX, (vRC, numInCat)) = FL.fold
        (lmePrepFrameOne (realToFrac . F.rgetField @Travel)
                         (const $ LA.fromList [1])
                         (F.rgetField @Rail)
        )
        railFrame
  resultEither <- runPIRLS_M $ do
    let (n, p) = LA.size mX
        rcRows = VB.length vRC
    liftIO $ print vRC
    when (rcRows /= n)
      $  P.throw
      $  "Only "
      <> (T.pack $ show rcRows)
      <> " in vRC but there are "
      <> (T.pack $ show n)
      <> " rows in the data!"
    let rowClassifier n = vRC VB.! n
        levels =
          VB.fromList [(numInCat, True, Nothing :: Maybe (VB.Vector Bool))]
    liftIO $ putStrLn $ show $ fmap colsForLevel levels
    liftIO $ putStrLn $ "y=" ++ show vY
    liftIO $ putStrLn "X="
    liftIO $ LA.disp 2 mX
    liftIO $ putStrLn $ "levels=" ++ show levels
    let smZ    = makeZ mX levels rowClassifier
    let (_, q) = SLA.dim smZ
    liftIO $ putStrLn $ "Z="
    liftIO $ LA.disp 2 $ asDense smZ
--  liftIO $ putStrLn "making A"
    checkProblem mX vY smZ
    let smA = makeA mX vY smZ
    liftIO $ putStrLn "A"
    liftIO $ LA.disp 2 $ asDense smA
    let makeST = makeSTF levels
    let (s, t) = makeST thV
    liftIO $ putStrLn "S"
    liftIO $ LA.disp 2 $ asDense s
    liftIO $ putStrLn "T"
    liftIO $ LA.disp 2 $ asDense t
    let aStar = makeAStar p q smA makeST thV
    liftIO $ putStrLn "A*"
    liftIO $ LA.disp 2 $ asDense aStar
    let (pd, ldL2, r, _) = profiledDeviance p q n smA makeST thV
    liftIO $ putStrLn $ "profiled Deviance=" ++ (show pd)
    liftIO $ putStrLn $ "r=" ++ show r
    liftIO $ putStrLn $ "ldL2=" ++ show ldL2
    (thS, (pdS, ldL2S, rS, ldL)) <- minimizeDeviance p q n levels smA makeST thV
    liftIO $ putStrLn $ "Solution: profiled Deviance=" ++ (show pdS)
    liftIO $ putStrLn $ "Solution: r=" ++ show rS
    liftIO $ putStrLn $ "Solution: ldL2=" ++ show ldL2S
    liftIO $ LA.disp 2 ldL
    let (beta, b) = parametersFromSolution p q makeST thS ldL rS
    liftIO $ do
      putStrLn $ "Fixed  (beta) =" ++ show beta
      putStrLn $ "Random (b) =" ++ show b
  case resultEither of
    Left  err -> putStrLn $ "Error: " ++ (T.unpack err)
    Right ()  -> putStrLn $ "Success!"




