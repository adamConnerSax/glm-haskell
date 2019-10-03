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

import qualified Data.IndexedSet               as IS
import qualified Numeric.GLM.Types             as GLM
import           Numeric.GLM.MixedModel
import qualified Numeric.SparseDenseConversions
                                               as SD

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Map                      as M
import           Data.Maybe                     ( catMaybes )
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA

import qualified Numeric.LinearAlgebra.Class   as SLA
{-
import qualified Data.Sparse.SpVector          as SLA
import qualified Numeric.LinearAlgebra.Sparse  as SLA
import           Numeric.LinearAlgebra.Sparse   ( (##)
                                                , (#^#)
                                                , (<#)
                                                , (#>)
                                                , (-=-)
                                                , (-||-)
                                                )
-}
import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL


import           System.IO.Unsafe               ( unsafePerformIO )

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

profiledDevianceLMM
  :: ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVector -- ^ theta
  -> IO
       ( Double
       , Double
       , LA.Vector Double
       , SLA.SpVector Double
       , SLA.SpVector Double
       , CholeskySolutions
       ) -- ^ (pd, sigma^2, beta, u, b) 
profiledDevianceLMM verbosity cf dt mm@(MixedModel (RegressionModel _ mX vY) _) reCalc@(RandomEffectCalculated smZ mkLambda) vTh
  = do
    let n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        lambda = mkLambda vTh
        smZS   = smZ SLA.## lambda --(smT SLA.## smS)
        smZSt  = SLA.transpose smZS
    cs@(CholeskySolutions smLth smRzx mRx svBeta svu) <-
      cholmodCholeskySolutionsLMM cf mm reCalc vTh
    when (verbosity == PDVAll) $ do
      putStrLn "Z"
      LA.disp 2 $ SD.toDenseMatrix smZ
      putStrLn "Lambda"
      LA.disp 2 $ SD.toDenseMatrix lambda
      putStrLn "Lth"
      LA.disp 2 $ SD.toDenseMatrix smLth
      putStrLn "Rzx"
      LA.disp 2 $ SD.toDenseMatrix smRzx
      putStrLn "Rx"
      LA.disp 2 mRx
    let svb     = lambda SLA.#> svu -- (smT SLA.## smS) SLA.#> svu
        vBeta   = SD.toDenseVector svBeta
        vDev    = vY - (mX LA.#> vBeta) - (SD.toDenseVector $ smZS SLA.#> svu)
        rTheta2 = (vDev LA.<.> vDev) + (svu SLA.<.> svu)
    let logLth        = logDetTriangularSM smLth
        (dof, logDet) = case dt of
          ML   -> (realToFrac n, logLth)
          REML -> (realToFrac (n - p), logLth + (logDetTriangularM mRx))
        pd     = (2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof)))
        sigma2 = rTheta2 / dof
    when (verbosity == PDVAll) $ do
      let (_, _, smP) = cf
      putStrLn $ "smP="
      LA.disp 1 $ SD.toDenseMatrix smP
      putStrLn $ "rTheta^2=" ++ show rTheta2
      putStrLn $ "2 * logLth=" ++ show (2 * logLth)
      putStrLn $ "2 * logDet=" ++ show (2 * logDet)
    when (verbosity == PDVAll || verbosity == PDVSimple) $ do
      putStrLn $ "pd(th=" ++ (show vTh) ++ ") = " ++ show pd
    return (pd, sigma2, vBeta, svu, svb, cs)

cholmodCholeskySolutionsLMM
  :: CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVector
  -> IO CholeskySolutions
cholmodCholeskySolutionsLMM cholmodFactor mixedModel randomEffCalcs vTh = do
  let (cholmodC, cholmodF, _) = cholmodFactor
      (MixedModel             (RegressionModel _ mX vY) levels  ) = mixedModel
      (RandomEffectCalculated smZ mkLambda) = randomEffCalcs
      lambda                  = mkLambda vTh
      smZS                    = smZ SLA.## lambda
  cholmodCholeskySolutions' cholmodFactor smZS mX NormalEquationsLMM mixedModel


-- ugh.  But I dont know a way in NLOPT to have bounds on some not others.
thetaLowerBounds :: FitSpecByGroup g -> NL.Bounds
thetaLowerBounds groupFSM =
  let negInfinity :: Double = negate $ 1 / 0
  in  NL.LowerBounds $ setCovarianceVector groupFSM 0 negInfinity --FL.fold fld levels

data MinimizeDevianceVerbosity = MDVNone | MDVSimple

data SolutionComponents = SolutionComponents { mRX :: LA.Matrix Double, vBeta :: LA.Vector Double, vTheta :: LA.Vector Double, vb :: LA.Vector Double }

-- b is an Enumeration of Effects/Predictors and g is an enumeration of groups
minimizeDevianceLMM
  :: SemC r
  => MinimizeDevianceVerbosity
  -> DevianceType
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVector -- ^ initial guess for theta
  -> P.Sem
       r
       ( LA.Vector Double
       , Double
       , Double
       , LA.Vector Double
       , LA.Vector Double
       , LA.Vector Double
       , CholeskySolutions
       ) -- ^ (theta, profiled_deviance, sigma2, beta, u, b, cholesky blocks)
minimizeDevianceLMM verbosity dt mixedModel@(MixedModel _ levels) reCalc@(RandomEffectCalculated smZ mkLambda) th0
  = do
    cholmodAnalysis <- cholmodAnalyzeProblem reCalc
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      pd x = profiledDevianceLMM pdv cholmodAnalysis dt mixedModel reCalc x
      obj x = unsafePerformIO $ fmap (\(d, _, _, _, _, _) -> d) $ pd x
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
        (pd, sigma2, vBeta, svu, svb, cs) <- pd thS
        return
          ( thS
          , pd
          , sigma2
          , vBeta
          , SD.toDenseVector svu
          , SD.toDenseVector svb
          , cs
          )

fixedEffectStatistics
  :: (Ord b, Enum b, Bounded b)
  => GLM.FixedEffects b
  -> Double
  -> CholeskySolutions
  -> GLM.FixedEffectStatistics b
fixedEffectStatistics fe sigma2 (CholeskySolutions _ _ mRX svBeta _) =
  let means = SD.toDenseVector svBeta
      covs  = LA.scale sigma2 $ (LA.inv mRX) LA.<> (LA.inv $ LA.tr mRX)
  in  GLM.FixedEffectStatistics fe means covs

effectParametersByGroup
  :: (Ord g, Show g, Show b)
  => GLM.RowClassifier g
  -> GLM.EffectsByGroup g b
  -> LA.Vector Double
  -> Either T.Text (GLM.EffectParametersByGroup g b)
effectParametersByGroup rc ebg vb = do
  let groups = IS.members $ GLM.groupIndices rc
      groupOffset group = do
        size     <- GLM.groupSize rc group
        nEffects <- IS.size <$> GLM.groupEffects ebg group
        return $ size * nEffects
  offsets <- FL.prescan FL.sum <$> traverse groupOffset groups
  let
    ep (group, offset) = do
      size    <- GLM.groupSize rc group
      effects <- GLM.groupEffects ebg group
      let nEffects = IS.size effects
          parameterMatrix =
            LA.fromRows
              $ VS.chunksOf nEffects
              $ VS.take (size * nEffects)
              $ VS.drop offset vb
      return (group, GLM.EffectParameters effects parameterMatrix)
  fmap M.fromList $ traverse ep $ zip groups offsets


effectCovariancesByGroup
  :: (Ord g, Show g, Enum b, Bounded b, Show b)
  => GLM.EffectsByGroup g b
  -> Double
  -> LA.Vector Double
  -> Either T.Text (GLM.EffectCovariancesByGroup g b)
effectCovariancesByGroup ebg sigma2 vTh = do
  let
    groups = M.keys ebg
    nElements n = n * (n + 1) `div` 2
    groupOffset group = nElements . IS.size <$> GLM.groupEffects ebg group
    ltIndices n = [ (r, c) | c <- [0 .. (n - 1)], r <- [c .. (n - 1)] ]
    ltAssocs n v = zip (ltIndices n) (VS.toList v)
    tF ((r, c), x) = if r > c then Just ((c, r), x) else Nothing
    tAssocs l = catMaybes $ fmap tF l
    makeLT n v = let l = ltAssocs n v in LA.assoc (n, n) 0 $ l ++ (tAssocs l)
  offsets <- FL.prescan FL.sum <$> traverse groupOffset groups
  let ecs (group, offset) = do
        effects <- GLM.groupEffects ebg group
        let nEffects = IS.size effects
            cm =
              LA.scale sigma2
                $ makeLT nEffects
                $ VS.take (nElements nEffects)
                $ VS.drop offset vTh
        return (group, GLM.GroupEffectCovariances effects cm)
  fmap M.fromList $ traverse ecs $ zip groups offsets


report
  :: (LA.Container LA.Vector Double, SemC r)
  => Int -- ^ p
  -> Int -- ^ q
  -> FitSpecByGroup g
  -> LA.Vector Double -- ^ y
  -> LA.Matrix Double -- ^ X
  -> SLA.SpMatrix Double -- ^ Z
  -> SLA.SpVector Double -- ^ beta
  -> SLA.SpVector Double -- ^ b
  -> P.Sem r ()
report p q groupFSM vY mX smZ svBeta svb = do
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
  let numberedGroups = zip [0 ..] (FL.fold FL.list groupFSM)
  liftIO $ mapM_
    (\(lN, l) -> putStrLn ("Level " ++ show lN) >> groupReport l vb)
    numberedGroups

