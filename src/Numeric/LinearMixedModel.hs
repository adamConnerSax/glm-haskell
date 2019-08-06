{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.LinearMixedModel
  ( module Numeric.LinearMixedModel
  , module Numeric.MixedModel
  )
where

import           Numeric.MixedModel
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
makeA :: MixedModel -> RandomEffectCalculated -> SLA.SpMatrix Double
makeA (MixedModel (RegressionModel mX vY) _) (RandomEffectCalculated smZ _) =
  let (n, p)     = LA.size mX
      yRows      = LA.size vY
      (zRows, q) = SLA.dim smZ
      smX        = SD.toSparseMatrix mX
      svY        = SD.toSparseVector vY
      zTz        = smZ SLA.#^# smZ
      zTx        = smZ SLA.#^# smX
      xTz        = smX SLA.#^# smZ
      xTx        = smX SLA.#^# smX
      m_zTy      = SLA.fromColsL $ [(-1) SLA..* (SLA.transpose smZ SLA.#> svY)]
      m_xTy      = SLA.fromColsL $ [(-1) SLA..* (SLA.transpose smX SLA.#> svY)]
      m_yTz = SLA.transpose $ SLA.fromColsL $ [(-1) SLA..* (svY SLA.<# smZ)]
      m_yTx = SLA.transpose $ SLA.fromColsL $ [(-1) SLA..* (svY SLA.<# smX)]
      yTy        = SLA.fromListSM (1, 1) [(0, 0, svY SLA.<.> svY)]
      c1         = (zTz SLA.-=- xTz) SLA.-=- m_yTz
      c2         = (zTx SLA.-=- xTx) SLA.-=- m_yTx
      c3         = (m_zTy SLA.-=- m_xTy) SLA.-=- yTy
  in  (c1 SLA.-||- c2) SLA.-||- c3


makeAStar
  :: Int -- ^ p (cols of X)
  -> Int -- ^ q (cols of Z) 
  -> SLA.SpMatrix Double -- ^ A
  -> (LA.Vector Double -> (SLA.SpMatrix Double, SLA.SpMatrix Double))
  -> CovarianceVector -- ^ theta
  -> SLA.SpMatrix Double
makeAStar p q smA mkST vTh =
  let (smS, smT) = mkST vTh
      tTs        = smT #^# smS
      st         = smS ## smT
      m1 =
        (tTs -||- SLA.zeroSM q p -||- SLA.zeroSM q 1)
          -=- (SLA.zeroSM p q -||- SLA.eye p -||- SLA.zeroSM p 1)
          -=- (SLA.zeroSM 1 q -||- SLA.zeroSM 1 p -||- SLA.eye 1)
      m2 =
        (st -||- SLA.zeroSM q p -||- SLA.zeroSM q 1)
          -=- (SLA.zeroSM p q -||- SLA.eye p -||- SLA.zeroSM p 1)
          -=- (SLA.zeroSM 1 q -||- SLA.zeroSM 1 p -||- SLA.eye 1)
      m3 =
        (SLA.eye q -||- SLA.zeroSM q p -||- SLA.zeroSM q 1)
          -=- (SLA.zeroSM p q -||- SLA.zeroSM p p -||- SLA.zeroSM p 1)
          -=- (SLA.zeroSM 1 q -||- SLA.zeroSM 1 p -||- SLA.zeroSM 1 1)
  in  (m1 ## smA ## m2) SLA.^+^ m3


makeAStar'
  :: SemC r
  => LA.Matrix Double  -- ^ X
  -> LA.Vector Double -- ^ y
  -> SLA.SpMatrix Double -- ^ Z
  -> (  CovarianceVector
     -> P.Sem r (SLA.SpMatrix Double, SLA.SpMatrix Double)
     )
  -> CovarianceVector -- ^ theta
  -> P.Sem r (SLA.SpMatrix Double)
makeAStar' mX vY smZ mkST vTh = do
  (smS, smT) <- mkST vTh
  let smX    = SD.toSparseMatrix mX
      svY    = SD.toSparseVector vY
      (_, q) = SLA.dim smZ
      smZS   = smZ ## smT ## smS
  liftIO $ do
    putStrLn "Z*="
    LA.disp 2 $ SD.toDenseMatrix smZS
  let zsTzs  = smZS #^# smZS
      zsTx   = smZS #^# smX
      xTzs   = smX #^# smZS
      xTx    = smX #^# smX
      m_zsTy = SLA.fromColsL $ [(-1) SLA..* (SLA.transpose smZS SLA.#> svY)]
      m_xTy  = SLA.fromColsL $ [(-1) SLA..* (SLA.transpose smX SLA.#> svY)]
      m_yTzs = SLA.transpose $ SLA.fromColsL $ [(-1) SLA..* (svY SLA.<# smZS)]
      m_yTx  = SLA.transpose $ SLA.fromColsL $ [(-1) SLA..* (svY SLA.<# smX)]
      yTy    = SLA.fromListSM (1, 1) [(0, 0, svY SLA.<.> svY)]
      c1     = ((zsTzs SLA.^+^ SLA.eye q) -=- xTzs) -=- m_yTzs
      c2     = (zsTx -=- xTx) -=- m_yTx
      c3     = (m_zsTy -=- m_xTy) -=- yTy
  return $ (c1 SLA.-||- c2) SLA.-||- c3



-- since A* is (p + q + 1) x (p + q + 1) so is L
profiledDeviance
  :: CholmodFactor
  -> DevianceType
  -> Int -- ^ p (cols of X)
  -> Int -- ^ q (cols of Z)
  -> Int -- ^ n (rows of X or Z)
  -> SLA.SpMatrix Double -- ^ A
  -> (  CovarianceVector
     -> (SLA.SpMatrix Double, SLA.SpMatrix Double)
     )
  -> CovarianceVector -- ^ theta
  -> (Double, Double, Double, SLA.SpMatrix Double)
profiledDeviance (cholmodC, cholmodF) dt p q n smA mkST th
  = let
      smAS = makeAStar p q smA mkST th
      cholL =
        CH.unsafeSpMatrixCholesky cholmodC cholmodF CH.SquareSymmetricLower smAS -- yikes!
      --diagL = SLA.extractDiag cholL
      getDiagX n = SLA.lookupWD_SM cholL (n, n)
      r = getDiagX $ q + p
      logDiagF n = FL.fold (FL.premap (log . getDiagX) FL.sum) [0 .. (n - 1)]
      (dof, logDiag) = case dt of
        ML   -> (n, logDiagF q)
        REML -> (n - p, logDiagF (q + p))
      ldL2 = realToFrac 2 * logDiag
      d =
        let rDOF = realToFrac dof
        in  ldL2 + rDOF * (1 + log (2 * pi * r * r / rDOF))
    in
      (d, ldL2, r, cholL)



-- ugh.  But I dont know a way in NLOPT to have bounds on some not others.
thetaLowerBounds :: Levels -> NL.Bounds
thetaLowerBounds levels =
  let negInfinity :: Double = negate $ 1 / 0
  in  NL.LowerBounds $ setCovarianceVector levels 0 negInfinity --FL.fold fld levels

minimizeDeviance
  :: SemC r
  => DevianceType
  -> Int -- ^ p (cols of X)
  -> Int -- ^ q (cols of Z)
  -> Int -- ^ n (rows of X or Z)
  -> Levels
  -> SLA.SpMatrix Double -- ^ A
  -> (  CovarianceVector
     -> (SLA.SpMatrix Double, SLA.SpMatrix Double)
     )
  -> CovarianceVector -- ^ initial guess for theta
  -> P.Sem
       r
       ( LA.Vector Double
       , ( Double
         , Double
         , Double
         , SLA.SpMatrix Double
         )
       )
minimizeDeviance dt p q n levels smA mkST th0 = do
  cholmodC <- liftIO CH.allocCommon
  liftIO $ CH.startC cholmodC
  cholmodF <- liftIO $ CH.spMatrixAnalyze cholmodC CH.SquareSymmetricLower smA
  let
    pd x = profiledDeviance (cholmodC, cholmodF) dt p q n smA mkST x
    obj x = (\(d, _, _, _) -> d) $ pd x
    stop           = NL.ObjectiveRelativeTolerance 1e-6 NL.:| []
    thetaLB        = thetaLowerBounds levels
    algorithm      = NL.BOBYQA obj [thetaLB] Nothing
    problem        = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
    expThetaLength = FL.fold
      (FL.premap
        (\l -> let e = effectsForLevel l in e + (e * (e - 1) `div` 2))
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
    Right (NL.Solution pdS thS result) -> do
      liftIO
        $  putStrLn
        $  show "Solution ("
        ++ show result
        ++ ") reached! At th="
      liftIO $ print thS
      return $ (thS, pd thS)

data MinimizeDevianceVerbosity = MDVNone | MDVSimple

minimizeDeviance2
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
minimizeDeviance2 verbosity dt mixedModel@(MixedModel _ levels) reCalc@(RandomEffectCalculated smZ mkLambda) th0
  = do
    (cholmodC, cholmodF) <- cholmodAnalyzeProblem reCalc
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      pd x = profiledDeviance2 pdv (cholmodC, cholmodF) dt mixedModel reCalc x
      obj x = unsafePerformIO $ fmap (\(d, _, _, _) -> d) $ pd x
      stop           = NL.ObjectiveAbsoluteTolerance 1e-6 NL.:| []
      thetaLB        = thetaLowerBounds levels
      algorithm      = NL.BOBYQA obj [thetaLB] Nothing
      problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
      expThetaLength = FL.fold
        (FL.premap
          (\l -> let e = effectsForLevel l in e + (e * (e - 1) `div` 2))
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
        (pd, vBeta, vu, vb) <- pd thS
        return (thS, pd, vBeta, vu, vb)


report
  :: (LA.Container LA.Vector Double, SemC r)
  => Int -- ^ p
  -> Int -- ^ q
  -> Levels
  -> LA.Vector Double -- ^ y
  -> LA.Matrix Double -- ^ X
  -> SLA.SpMatrix Double -- ^ Z
  -> SLA.SpVector Double -- ^ beta
  -> SLA.SpVector Double -- ^ b
  -> P.Sem r ()
report p q levels vY mX smZ svBeta svb = do
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
    levelReport l@(LevelSpec n b _) b' = do
      putStrLn $ show n ++ " groups"
      when b $ reportStats "Intercept" $ VS.take n b'
      let (numSlopes, bS) = case b of
            True  -> (effectsForLevel l - 1, VS.drop n b')
            False -> (effectsForLevel l, b')
      mapM_
        (\s -> reportStats ("Slope " <> (T.pack $ show s))
                           (VS.take n $ VS.drop (n * s) bS)
        )
        [0 .. (numSlopes - 1)]
  liftIO $ putStrLn $ "p=" ++ show p ++ "; q=" ++ show q
  liftIO $ reportStats "Residual" vEps
  let numberedLevels = zip [0 ..] (VB.toList levels)
  liftIO $ mapM_
    (\(lN, l) -> putStrLn ("Level " ++ show lN) >> levelReport l vb)
    numberedLevels

---

parametersFromSolution
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double, SemC r)
  => DevianceType -- ^ unused for now??
  -> Int
  -> Int
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a)) -- ^ makeST
  -> LA.Vector a -- ^ solution theta
  -> (SLA.SpMatrix a, SLA.SpMatrix a) -- ^ solution (P,L)
  -> a -- ^ solution r
  -> P.Sem
       r
       (SLA.SpVector a, SLA.SpVector a, SLA.SpVector a) -- ^ beta, b and b* (also called u)                       
parametersFromSolution dt p q makeST th (permSM, ldL) r = do
  let sv :: SLA.SpVector Double =
        SLA.fromListDenseSV (p + q + 1) $ L.replicate (p + q) 0 ++ [r]
  solV <- liftIO $ SLA.triUpperSolve (SLA.transposeSM ldL SLA.## permSM) sv
  let bStar  = SLA.takeSV q $ solV
      beta   = SLA.takeSV p $ SLA.dropSV q $ solV
      (s, t) = makeST th
      b      = (t SLA.## s) SLA.#> bStar
  return (beta, b, bStar)
