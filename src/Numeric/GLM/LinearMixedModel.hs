{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.GLM.LinearMixedModel
  ( module Numeric.GLM.LinearMixedModel
  , module Numeric.GLM.MixedModel
  )
where

import qualified Numeric.GLM.Types             as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.SparseDenseConversions
                                               as SD
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
--import qualified Data.Sparse.SpVector          as SLA
import qualified Numeric.LinearAlgebra.Class   as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import           Numeric.LinearAlgebra.Sparse   ( (##)
                                                , (#^#)
                                                , (<#)
                                                , (#>)
                                                , (-=-)
                                                , (-||-)
                                                )
import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL

import           System.IO.Unsafe               ( unsafePerformIO )

--import           Data.Either                    ( partitionEithers )
import qualified Data.List                     as L
--import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS


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


-- ugh.  But I dont know a way in NLOPT to have bounds on some not others.
thetaLowerBounds :: GroupFitSpecs -> NL.Bounds
thetaLowerBounds levels =
  let negInfinity :: Double = negate $ 1 / 0
  in  NL.LowerBounds $ setCovarianceVector levels 0 negInfinity --FL.fold fld levels

data MinimizeDevianceVerbosity = MDVNone | MDVSimple

minimizeDeviance
  :: SemC r
  => MinimizeDevianceVerbosity
  -> DevianceType
  -> MixedModel
  -> RandomEffectCalculated
  -> CovarianceVector -- ^ initial guess for theta
  -> P.Sem
       r
       ( LA.Vector Double
       , Double
       , LA.Vector Double
       , LA.Vector Double
       , LA.Vector Double
       ) -- ^ (theta, profiled_deviance, beta, u, b)
minimizeDeviance verbosity dt mixedModel@(MixedModel _ levels) reCalc@(RandomEffectCalculated smZ mkLambda) th0
  = do
    cholmodAnalysis <- cholmodAnalyzeProblem reCalc
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      pd x = profiledDeviance pdv cholmodAnalysis dt mixedModel reCalc x
      obj x = unsafePerformIO $ fmap (\(d, _, _, _) -> d) $ pd x
      stop           = NL.ObjectiveAbsoluteTolerance 1e-6 NL.:| []
      thetaLB        = thetaLowerBounds levels
      algorithm      = NL.BOBYQA obj [thetaLB] Nothing
      problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
      expThetaLength = FL.fold
        (FL.premap
          (\l -> let e = effectsForGroup l in e + (e * (e - 1) `div` 2))
          FL.sum
        )
        levels
    when (LA.size th0 /= expThetaLength)
      $  P.throw
      $  "guess for theta has "
      <> (T.pack $ show $ LA.size th0)
      <> " entries but should have "
      <> (T.pack $ show expThetaLength)
      <> "."
    let eSol = NL.minimizeLocal problem th0
    case eSol of
      Left  result                       -> P.throw (T.pack $ show result)
      Right (NL.Solution pdS thS result) -> liftIO $ do
        putStrLn
          $  "Solution ("
          ++ show result
          ++ ") reached! At th="
          ++ show thS
        (pd, vBeta, svu, svb) <- pd thS
        return (thS, pd, vBeta, SD.toDenseVector svu, SD.toDenseVector svb)


report
  :: (LA.Container LA.Vector Double, SemC r)
  => Int -- ^ p
  -> Int -- ^ q
  -> GroupFitSpecs
  -> LA.Vector Double -- ^ y
  -> LA.Matrix Double -- ^ X
  -> SLA.SpMatrix Double -- ^ Z
  -> SLA.SpVector Double -- ^ beta
  -> SLA.SpVector Double -- ^ b
  -> P.Sem r ()
report p q groupFSs vY mX smZ svBeta svb = do
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
  let numberedGroups = zip [0 ..] (VB.toList groupFSs)
  liftIO $ mapM_
    (\(lN, l) -> putStrLn ("Level " ++ show lN) >> groupReport l vb)
    numberedGroups

