{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE DeriveDataTypeable  #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.MixedModel where

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH
import qualified Numeric.SparseDenseConversions
                                               as SD
import qualified Numeric.GLM.ProblemTypes      as GLM
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Data.IndexedSet               as IS

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Knit.Effect.Logger            as P

--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import qualified Control.Exception             as X
import           Control.Monad                  ( when
                                                , join
                                                )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Data.Array                    as A
import           Data.Function                  ( (&) )
import qualified Data.List                     as L
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA
import           Data.Typeable                  ( Typeable )

import           Numeric.LinearAlgebra.Sparse   ( (#>) )

import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL
import           System.IO.Unsafe               ( unsafePerformIO )



import qualified Data.Map                      as M
import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

import           System.IO                      ( hSetBuffering
                                                , stdout
                                                , BufferMode(..)
                                                )

runEffectsVerboseIO :: GLM.GLMEffects a -> IO (Either GLM.GLMError a)
runEffectsVerboseIO action =
  action
    & P.filteredLogEntriesToIO P.logAll
    & P.errorToIOFinal
    & P.embedToFinal
    & P.runFinal


runEffectsIO :: GLM.GLMEffects a -> IO (Either GLM.GLMError a)
runEffectsIO action =
  action
    & P.filteredLogEntriesToIO P.nonDiagnostic
    & P.errorToIOFinal
    & P.embedToFinal
    & P.runFinal


unsafePerformEffects :: GLM.GLMEffects a -> a -> a
unsafePerformEffects action def =
  action
    & runEffectsIO
    & fmap (either X.throwIO return)
    & join
    & (\x -> X.catch
        x
        (\(e :: X.SomeException) -> putStrLn ("Error: " ++ show e) >> return def
        )
      )
    & unsafePerformIO

lmmOptimizerToNLOPT
  :: GLM.LMMOptimizer
  -> (  NL.Objective
     -> [NL.Bounds]
     -> Maybe NL.InitialStep
     -> NL.LocalAlgorithm
     )
lmmOptimizerToNLOPT GLM.LMM_BOBYQA       = NL.BOBYQA
lmmOptimizerToNLOPT GLM.LMM_NELDERMEAD   = NL.NELDERMEAD
lmmOptimizerToNLOPT GLM.LMM_SBPLX        = NL.SBPLX
lmmOptimizerToNLOPT GLM.LMM_NEWUOA_BOUND = NL.NEWUOA_BOUND

data DevianceType = ML | REML deriving (Show, Eq)

effectsForGroup :: GLM.GroupFitSpec -> Int
effectsForGroup (GLM.GroupFitSpec _ b vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForGroup #-}

colsForGroup :: GLM.GroupFitSpec -> Int
colsForGroup l = GLM.nCategories l * effectsForGroup l
{-# INLINABLE colsForGroup #-}

data ParameterEstimates =
  ParameterEstimates
  { pEstimate :: LA.Vector Double
  , pCovariance :: LA.Matrix Double
  }

type FixedParameterEstimates = ParameterEstimates
type GroupParameterEstimates = VB.Vector FixedParameterEstimates

data MinimizeDevianceVerbosity = MDVNone | MDVSimple

minimizeDeviance
  :: GLM.EffectsIO r
  => MinimizeDevianceVerbosity
  -> DevianceType
  -> GLM.MixedModel b g
  -> GLM.RandomEffectCalculated
  -> GLM.CovarianceVec -- ^ initial guess for theta
  -> P.Sem
       r
       ( ( GLM.CovarianceVec
         , Double
         , Double
         , GLM.BetaU
         , LA.Vector Double
         , CholeskySolutions
         )
       , LA.Vector Double
       , CholmodFactor
       ) -- ^ ((theta, profiled_deviance, sigma2, (beta, u), b, Cholesky blocks), mu, Cholmod Analysis)
minimizeDeviance verbosity dt mm reCalc th0 =
  P.wrapPrefix "minimizeDeviance" $ do
    P.logLE P.Diagnostic "Starting..."
    when (GLM.generalized mm && (dt == REML)) $ P.throw $ GLM.OtherGLMError
      "Can't use REML for generalized models."
    P.logLE P.Diagnostic "Beginning Cholmod analysis of Z matrix structure."
    cholmodAnalysis <- cholmodAnalyzeProblem reCalc
    P.logLE P.Diagnostic "Finished Cholmod analysis of Z matrix structure."
    res@(vTh, _, dev2, betaU, _, _) <- minimizeDevianceInner verbosity
                                                             dt
                                                             mm
                                                             reCalc
                                                             cholmodAnalysis
                                                             th0
    let
      mX    = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      zStar = GLM.makeZS reCalc vTh
      vEta  = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom betaU)
      vMu = LA.cmap (GLM.invLink $ GLM.linkFunction (linkFunctionType mm)) vEta
    return (res, vMu, cholmodAnalysis)
{-    
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      (objectiveOs, finalOs) = if GLM.generalized mm
        then (Optim_GLMM_PIRLS, Optim_GLMM_Final)
        else (Optim_LMM, Optim_LMM)
      pd
        :: GLM.EffectsIO r
        => OptimizationStage
        -> GLM.CovarianceVec
        -> P.Sem
             r
             ( Double
             , Double
             , GLM.BetaU
             , SLA.SpVector Double
             , CholeskySolutions
             )
      pd os x = profiledDeviance pdv cholmodAnalysis dt mm reCalc os x
      obj x = unsafePerformEffects
        (fmap (\(d, _, _, _, _) -> d) $ pd objectiveOs x)
        0
      levels    = GLM.mmsFitSpecByGroup $ GLM.mixedModelSpec mm
      thetaLB   = thetaLowerBounds levels
      algorithm = (lmmOptimizerToNLOPT $ GLM.lmmOptimizer $ GLM.lmmControls mm)
        obj
        [thetaLB]
        Nothing
      stop =
        NL.ObjectiveAbsoluteTolerance
            (GLM.lmmOptimizerTolerance $ GLM.lmmControls mm)
          NL.:| []
      problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
      expThetaLength = FL.fold
        (FL.premap
          (\l -> let e = effectsForGroup l in e + (e * (e - 1) `div` 2))
          FL.sum
        )
        levels
    when (LA.size th0 /= expThetaLength)
      $  P.throw
      $  GLM.OtherGLMError
      $  "guess for theta has "
      <> (T.pack $ show $ LA.size th0)
      <> " entries but should have "
      <> (T.pack $ show expThetaLength)
      <> "."
    let eSol = NL.minimizeLocal problem th0
    case eSol of
      Left result -> P.throw (GLM.OtherGLMError $ T.pack $ show result)
      Right (NL.Solution pdS thS result) -> do
        P.logLE P.Info
          $  "Solution ("
          <> (T.pack $ show result)
          <> ") reached! At th="
          <> (T.pack $ show thS)
        (pdVal, sigma2, betaU, svb, cs) <- pd finalOs thS
        return (thS, pdVal, sigma2, betaU, SD.toDenseVector svb, cs)
-}

minimizeDevianceInner
  :: GLM.EffectsIO r
  => MinimizeDevianceVerbosity
  -> DevianceType
  -> GLM.MixedModel b g
  -> GLM.RandomEffectCalculated
  -> CholmodFactor
  -> GLM.CovarianceVec -- ^ initial guess for theta
  -> P.Sem
       r
       ( GLM.CovarianceVec
       , Double
       , Double
       , GLM.BetaU
       , LA.Vector Double
       , CholeskySolutions
       ) -- ^ (theta, profiled_deviance, sigma2, beta, u, b, cholesky blocks)
minimizeDevianceInner verbosity dt mm reCalc cf th0 = P.wrapPrefix "minimizeDevianceInner" $ do
  let
    pdv = case verbosity of
      MDVNone   -> PDVNone
      MDVSimple -> PDVSimple
    (objectiveOs, finalOs) = if GLM.generalized mm
      then (Optim_GLMM_PIRLS, Optim_GLMM_Final)
      else (Optim_LMM, Optim_LMM)
    pd
      :: GLM.EffectsIO r
      => OptimizationStage
      -> GLM.CovarianceVec
      -> P.Sem
           r
           ( Double
           , Double
           , GLM.BetaU
           , SLA.SpVector Double
           , CholeskySolutions
           )
    pd os x = profiledDeviance pdv cf dt mm reCalc os x
    obj x =
      unsafePerformEffects (fmap (\(d, _, _, _, _) -> d) $ pd objectiveOs x) 0
    levels    = GLM.mmsFitSpecByGroup $ GLM.mixedModelSpec mm
    thetaLB   = thetaLowerBounds levels
    algorithm = (lmmOptimizerToNLOPT $ GLM.lmmOptimizer $ GLM.lmmControls mm)
      obj
      [thetaLB]
      Nothing
    stop =
      NL.ObjectiveAbsoluteTolerance
          (GLM.lmmOptimizerTolerance $ GLM.lmmControls mm)
        NL.:| []
    problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
    eSol    = NL.minimizeLocal problem th0
  case eSol of
    Left result -> P.throw (GLM.OtherGLMError $ T.pack $ show result)
    Right (NL.Solution pdS thS result) -> do
      P.logLE P.Diagnostic
        $  "Solution ("
        <> (T.pack $ show result)
        <> ") reached! At th="
        <> (T.pack $ show thS)
      (pdVal, sigma2, betaU, svb, cs) <- pd finalOs thS
      return (thS, pdVal, sigma2, betaU, SD.toDenseVector svb, cs)


-- ugh.  But I dont know a way in NLOPT to have bounds on some not others.
thetaLowerBounds :: GLM.FitSpecByGroup g -> NL.Bounds
thetaLowerBounds groupFSM =
  let negInfinity :: Double = negate $ 1 / 0
  in  NL.LowerBounds $ setCovarianceVector groupFSM 0 negInfinity --FL.fold fld levels

setCovarianceVector
  :: GLM.FitSpecByGroup g -> Double -> Double -> GLM.CovarianceVec
setCovarianceVector groupFSM diag offDiag = FL.fold fld groupFSM
 where
  fld = FL.Fold
    (\bs l ->
      let e       = effectsForGroup l
          entries = e * (e + 1) `div` 2
          col n = n `div` e
          row n = col n + n `mod` e
          set n = if (row n == col n) then diag else offDiag
      in  bs ++ fmap set (take entries $ iterate (+ 1) 0)
    )
    []
    LA.fromList

{-
Z is the random effects model matrix
like X, the fixed effect model matrix,
it has a row for each observation.
Z has a column for each random effect:
A column of ones for a random intercept and
a column with the predictor from X for a random
slope.
-}
makeZ
  :: forall g r
   . (Enum g, Bounded g, A.Ix g, GLM.Effects r)
  => GLM.FixedPredictors
  -> GLM.FitSpecByGroup g
  -> GLM.RowClassifier g
  -> P.Sem r GLM.RandomEffectModelMatrix
makeZ mX groupFSM rc = do
  let
    maybeError :: T.Text -> Maybe a -> P.Sem r a
    maybeError err = maybe (P.throw $ GLM.OtherGLMError $ err) return
    (nO, nP) = LA.size mX
    k        = FL.fold FL.length groupFSM -- number of levels
--    groupSize g = fmap nCategories $ M.lookup g groupFSM 
    q        = FL.fold FL.sum $ fmap colsForGroup groupFSM -- total number of columns in Z
    predictor rowIndex fixedEffectIndex =
      mX `LA.atIndex` (rowIndex, fixedEffectIndex)
    -- construct Z for level as a fold over rows
    entries
      :: Int -> g -> GLM.GroupFitSpec -> Int -> P.Sem r [(Int, Int, Double)]
    entries colOffset group groupFS rowIndex = do
      categoryN <-
        maybeError "Error in categorNumberFromRowIndex"
          $ GLM.categoryNumberFromRowIndex rc rowIndex group
      let
        entryOffset = colOffset + (effectsForGroup groupFS * categoryN)
        (intercept, slopeOffset) = if GLM.groupIntercept groupFS
          then ([(rowIndex, entryOffset, 1)], entryOffset + 1)
          else ([], entryOffset)
        slopeF :: FL.Fold (Bool, Int) [(Int, Int, Double)]
        slopeF = FL.Fold
          (\(mes, ziCol) (b, fixedEffectIndex) -> if b
            then
              ( (rowIndex, ziCol, predictor rowIndex fixedEffectIndex) : mes
              , ziCol + 1
              )
            else (mes, ziCol)
          )
          ([], slopeOffset)
          fst
        slopes = case GLM.groupSlopes groupFS of
          Just slopeV -> FL.fold slopeF $ VB.zip slopeV (VB.generate nP id)
          Nothing     -> []
      return $ intercept ++ slopes
    ziFoldM
      :: g
      -> Int
      -> GLM.GroupFitSpec
      -> FL.FoldM (P.Sem r) Int [(Int, Int, Double)]
    ziFoldM group groupOffset groupFS =
      let q = colsForGroup groupFS
      in  FL.FoldM
            (\mes n -> fmap (mes ++) $ entries groupOffset group groupFS n)
            (return [])
            return
    groups       = IS.members $ GLM.groupIndices rc -- [minBound .. maxBound] --Numbers = VB.generate (VB.length groupFSs) id
    groupOffsets = FL.prescan FL.sum $ fmap colsForGroup $ M.elems groupFSM
  zisM <- FL.foldM
    ( sequenceA
    $ fmap (\(group, gOffset, gSpec) -> ziFoldM group gOffset gSpec)
    $ zip3 groups groupOffsets (M.elems groupFSM)
    )
    [0 .. (nO - 1)]
  return $ SLA.fromListSM (nO, q) $ concat zisM

-- TODO: Add levels checks
-- TODO: Add Lambda checks
-- TODO: parse, don't validate!  Which I think means typed sizes here.
checkProblem
  :: GLM.Effects r
  => GLM.MixedModel b g
  -> GLM.RandomEffectCalculated
  -> P.Sem r ()
checkProblem mm reCalc = do
  let mX         = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vY         = GLM.rmsObservations $ GLM.regressionModelSpec mm
      smZ        = GLM.recModelMatrix reCalc
      (n, _)     = LA.size mX
      yRows      = LA.size vY
      (zRows, _) = SLA.dim smZ
  when (zRows /= n)
    $  P.throw
    $  GLM.OtherGLMError
    $  (T.pack $ show n)
    <> " cols in X but "
    <> (T.pack $ show zRows)
    <> " cols in Z"
  when (yRows /= n)
    $  P.throw
    $  GLM.OtherGLMError
    $  (T.pack $ show n)
    <> " cols in X but "
    <> (T.pack $ show yRows)
    <> " entries in Y"
  return () --(MixedModel (RegressionModel _ mX vY) _) (RandomEffectCalculated smZ _)
-- For each level, i, with p_i random effects,
-- Lambda is built from lower-triangular blocks
-- which are p_i x p_i.  Theta contains the values
-- for each block, in col major order
makeLambda :: GLM.FitSpecByGroup g -> (GLM.CovarianceVec -> GLM.LambdaMatrix)
makeLambda groupFSM =
  let blockF vTH = FL.Fold
        (\(bs, vx) l' -> ((l', vx) : bs, VS.drop (effectsForGroup l') vx))
        ([], vTH)
        (reverse . fst)
      blockData vTh = FL.fold (blockF vTh) groupFSM
      templateBlock l vTh =
        let pl  = effectsForGroup l
            lTh = VS.toList vTh
            lts vx = fmap (\((r, c), v) -> (r, c, v)) $ zip
              ([ (r, c) | c <- [0 .. (pl - 1)], r <- [c .. (pl - 1)] ])
              vx
        in  SLA.fromListSM (pl, pl) $ lts lTh
      perGroup (l, vx) = replicate (GLM.nCategories l) $ templateBlock l vx
      allDiags vTh = concat $ fmap perGroup $ blockData vTh
  in  (\vTh -> SLA.fromBlocksDiag (allDiags vTh))

xTx :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTx smX = (smX SLA.#^# smX)

xTxPlusI :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTxPlusI smX = (smX SLA.#^# smX) SLA.^+^ (SLA.eye $ SLA.ncols smX)

logDetTriangularSM :: RealFloat a => SLA.SpMatrix a -> a
logDetTriangularSM smX =
  let f (_, _, x) = log x
      diag (r, c, _) = (r == c)
  in  FL.fold (FL.premap f FL.sum) $ L.filter diag $ SLA.toListSM smX

logDetTriangularM :: (RealFloat a, LA.Container VS.Vector a) => LA.Matrix a -> a
logDetTriangularM mX =
  let f n = log $ mX `LA.atIndex` (n, n)
      (rows, _) = LA.size mX
  in  FL.fold (FL.premap f FL.sum) [0 .. (rows - 1)]

type CholmodFactor = ( CH.ForeignPtr CH.Common -- ^ pre-allocated CHOMOD common space
                     , CH.ForeignPtr CH.Factor -- ^ precomputed pattern work on Z
                     , GLM.PMatrix -- ^ permutation matrix from above
                     )

cholmodMakeCommon :: GLM.EffectsIO r => P.Sem r (CH.ForeignPtr CH.Common)
cholmodMakeCommon = liftIO $ do
  cholmodC <- CH.allocCommon
  CH.startC cholmodC
  CH.setFinalLL 1 cholmodC
  return cholmodC

cholmodAnalyzeProblem
  :: GLM.EffectsIO r => GLM.RandomEffectCalculated -> P.Sem r CholmodFactor
cholmodAnalyzeProblem reCalc = do
  cholmodC        <- cholmodMakeCommon
  (cholmodF, smP) <-
    liftIO
    $ CH.spMatrixAnalyzeWP cholmodC CH.SquareSymmetricLower
    $ xTx
    $ GLM.recModelMatrix reCalc
--  liftIO $ CH.printFactor cholmodF "After analyze, before factorize" cholmodC
--  liftIO $ putStrLn "smP=" >> (LA.disp 1 $ SD.toDenseMatrix smP)
  return (cholmodC, cholmodF, smP)

data ProfiledDevianceVerbosity = PDVNone | PDVSimple | PDVAll deriving (Show, Eq)
data OptimizationStage = Optim_LMM | Optim_GLMM_PIRLS | Optim_GLMM_Final deriving (Show, Eq)

-- this one is IO since we'll need to unsafePerformIO it
profiledDeviance
  :: GLM.EffectsIO r
  => ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> GLM.MixedModel b g
  -> GLM.RandomEffectCalculated
  -> OptimizationStage
  -> GLM.CovarianceVec -- ^ theta
  -> P.Sem
       r
       ( Double
       , Double
       , GLM.BetaU
       , SLA.SpVector Double
       , CholeskySolutions
       ) -- ^ (pd, sigma^2, beta, u, b) 
profiledDeviance verbosity cf dt mm reCalc os vTh =
  P.wrapPrefix "profiledDeviance" $ do
    P.logLE P.Diagnostic $ "here"
    liftIO $ hSetBuffering stdout NoBuffering
    let mX     = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
        smZ    = GLM.recModelMatrix reCalc
        lambda = GLM.recMkLambda reCalc vTh
        vY     = GLM.rmsObservations $ GLM.regressionModelSpec mm
        n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        zStar  = GLM.makeZS reCalc vTh

    P.logLE P.Diagnostic "Starting. Calling getCholeskySolutions."
    (cs@(CholeskySolutions smLth smRzx mRxx), betaU) <- getCholeskySolutions
      cf
      mm
      zStar
      os
      vTh
    when (verbosity == PDVAll) $ liftIO $ do
      putStrLn "Z"
      LA.disp 2 $ SD.toDenseMatrix smZ
      putStrLn "Lambda"
      LA.disp 2 $ SD.toDenseMatrix lambda
      putStrLn "Lth"
      LA.disp 2 $ SD.toDenseMatrix smLth
      putStrLn "Rzx"
      LA.disp 2 $ SD.toDenseMatrix smRzx
      putStrLn "Rxx"
      LA.disp 2 mRxx
    let osFinal = if (GLM.generalized mm) then Optim_GLMM_Final else Optim_LMM
        vEta    = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom betaU)
    (pd, rTheta2) <- profiledDeviance' verbosity
                                       dt
                                       mm
                                       zStar
                                       vTh
                                       osFinal
                                       cs
                                       vEta
                                       (GLM.bu_svU betaU)
    let dof = case dt of
          ML   -> realToFrac n
          REML -> realToFrac (n - p)
        sigma2 = rTheta2 / dof
    when (verbosity == PDVAll) $ liftIO $ do
      let (_, _, smP) = cf
      putStrLn $ "smP="
      LA.disp 1 $ SD.toDenseMatrix smP
      putStrLn $ "rTheta^2=" ++ show rTheta2
      let logLth = logDetTriangularSM smLth
          logDet = case dt of
            ML   -> logLth
            REML -> logLth + logDetTriangularM mRxx
      putStrLn $ "2 * logLth=" ++ show (2 * logLth)
      putStrLn $ "2 * logDet=" ++ show (2 * logDet)
    when (verbosity == PDVAll || verbosity == PDVSimple) $ liftIO $ do
      putStrLn $ "pd(th=" ++ (show vTh) ++ ") = " ++ show pd
    let svb = lambda SLA.#> (GLM.bu_svU betaU)
    return (pd, sigma2, betaU, svb, cs)



profiledDeviance'
  :: GLM.EffectsIO r
  => ProfiledDevianceVerbosity
  -> DevianceType
  -> GLM.MixedModel b g
  -> GLM.ZStarMatrix
  -> GLM.CovarianceVec
  -> OptimizationStage
  -> CholeskySolutions
  -> GLM.EtaVec
  -> GLM.UVec
  -> P.Sem r (Double, Double) -- profiled deviance, deviance + u'u 
profiledDeviance' pdv dt mm zStar vTh os chol vEta svU =
  P.wrapPrefix "profiledDeviance'" $ do
    P.logLE P.Diagnostic
      $  "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vTh)
      <> ")..."
--    P.logLE P.Diagnostic $ "vEta=" <> (T.pack $ show vEta)
--    P.logLE P.Diagnostic $ "svU=" <> (T.pack $ show svU)
    let
      mX       = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vY       = GLM.rmsObservations $ GLM.regressionModelSpec mm
      (od, vW) = case mm of
        GLM.LinearMixedModel _ ->
          (GLM.Normal, (LA.fromList $ L.replicate (LA.size vY) 1.0))
        GLM.GeneralizedLinearMixedModel glmmSpec ->
          ( GLM.glmmsObservationsDistribution glmmSpec
          , GLM.glmmsWeights glmmSpec
          )
      logLth = logDetTriangularSM (lTheta chol)
--      smZS   = GLM.makeZS reCalc vTh
      vMu = LA.cmap (GLM.invLink $ GLM.linkFunction (linkFunctionType mm)) vEta
      uDotu = svU SLA.<.> svU
    P.logLE P.Diagnostic
      $  "logLth="
      <> (T.pack $ show logLth)
      <> "; u.u="
      <> (T.pack $ show uDotu)
    (pd, rTheta2) <- case mm of
      GLM.LinearMixedModel _ -> do
        when (os /= Optim_LMM)
          $ P.throw
          $ GLM.OtherGLMError
              "In profiledDeviance': OptimzationStage must be Optim_LMM for all LMM calls"
        let
          n             = LA.size vY
          vW            = VS.replicate n 1.0
          devResidual   = GLM.deviance (GLM.Normal) vW vY vMu
          rTheta2       = devResidual + uDotu
          (_  , p     ) = LA.size mX
          (dof, logDet) = case dt of
            ML -> (realToFrac n, logLth)
            REML ->
              (realToFrac (n - p), logLth + (logDetTriangularM $ rXX chol))
        P.logLE P.Diagnostic $ "LMM devResid=: " <> (T.pack $ show devResidual)
        return
          $ ((2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof))), rTheta2)
      GLM.GeneralizedLinearMixedModel _ -> do
        let devResidual = GLM.deviance od vW vY vMu --GLM.devianceCondAbs od vW vY vMu
            aic         = GLM.aicR od vW vY vMu devResidual
            rTheta2     = devResidual + uDotu
        P.logLE P.Diagnostic
          $  "GLMM: devResid= "
          <> (T.pack $ show devResidual)
          <> "; aicR="
          <> (T.pack $ show aic)
        case os of
          Optim_LMM ->
            P.throw
              $ GLM.OtherGLMError
              $ "In profiledDeviance': OptimizationStage was Optim_LMM in a GLMM call"
          Optim_GLMM_PIRLS -> return $ (devResidual + uDotu, rTheta2)
          Optim_GLMM_Final -> return $ (aic + uDotu + 2 * logLth, rTheta2)
    P.logLE P.Diagnostic $ "; pd=" <> (T.pack $ show pd)
    return (pd, rTheta2)

upperTriangular :: Int -> Int -> a -> Bool
upperTriangular r c _ = (r <= c)

lowerTriangular :: Int -> Int -> a -> Bool
lowerTriangular r c _ = (r >= c)

--type LTheta = SLA.SpMatrix Double
--type Rzx = SLA.SpMatrix Double
--type Rxx = LA.Matrix Double
--type Beta = LA.Vector Double

data CholeskySolutions = CholeskySolutions { lTheta :: GLM.LMatrix
                                           , rZX :: GLM.RzxMatrix
                                           , rXX :: GLM.RxxMatrix
                                           } deriving (Show)

data NormalEquations = NormalEquationsLMM | NormalEquationsGLMM GLM.WMatrix GLM.LinkFunction GLM.EtaVec GLM.UVec

normalEquationsRHS
  :: GLM.EffectsIO r
  => NormalEquations
  -> GLM.UMatrix
  -> GLM.VMatrix
  -> GLM.Observations
  -> P.Sem r (SLA.SpVector Double, LA.Vector Double)
normalEquationsRHS NormalEquationsLMM smUt mVt vY =
  P.wrapPrefix "normalEquationsRHS" $ do
    let svRhsZ = smUt SLA.#> SD.toSparseVector vY
        vRhsX  = mVt LA.#> vY
--    P.logLE P.Diagnostic $ "RHS=" <> (T.pack $ show (svRhsZ, vRhsX))
    return (svRhsZ, vRhsX) --(SLA.transpose smU SLA.#> SD.toSparseVector vY, LA.tr mV LA.#> vY)

normalEquationsRHS (NormalEquationsGLMM vVarW lf vEta svU) smUt mVt vY =
  P.wrapPrefix "normalEquationsRHS" $ do
{-  
    P.logLE P.Diagnostic "Starting..."

    P.logLE P.Diagnostic $ "vVarW=" <> (T.pack $ show vVarW)

    P.logLE P.Diagnostic $ "svU=" <> (T.pack $ show svU)
    P.logLE P.Diagnostic $ "smU=" <> (T.pack $ show smU)
    P.logLE P.Diagnostic $ "mV=" <> (T.pack $ show mV)
    P.logLE P.Diagnostic $ "vMu=" <> (T.pack $ show vMu)
    P.logLE P.Diagnostic $ "vY=" <> (T.pack $ show vY)
-}
    let vMu         = VS.map (GLM.invLink lf) vEta
--        vMuEta      = VS.map (GLM.derivInv lf) vEta
--        vYMinusMu   = vY - vMu
        vWtYMinusMu = VS.zipWith3 (\wt y mu -> sqrt wt * (y - mu)) vVarW vY vMu
  {-
      vWrkResids  = VS.zipWith3 (\y mu muEta -> (y - mu) / muEta) vY vMu vMuEta
      vWrkResp    = VS.zipWith (+) vEta vWrkResids
      vWtWrkResp  = VS.zipWith3
        (\wt muEta wrkResp -> muEta * (sqrt wt) * wrkResp)
        vVarW
        vMuEta
        vWrkResp
-}
--        vZ          = VS.zipWith (*) vMuEta vWtYMinusMu
        vX          = vWtYMinusMu
{-      
    P.logLE P.Diagnostic
      $  "vWtYminusMu="
      <> (T.pack $ show vWtYMinusMu)
      <> "\nvWtWrkResp="
      <> (T.pack $ show vWtWrkResp)
-}
    let svUX   = smUt SLA.#> (SD.toSparseVector vX)
        svRhsZ = svUX SLA.^-^ svU
        vRhsX  = mVt LA.#> vX
--    P.logLE P.Diagnostic $ "RHS=" <> (T.pack $ show (svRhsZ, vRhsX))
    return (svRhsZ, vRhsX)

cholmodCholeskySolutionsLMM
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModelSpec b g
  -> GLM.ZStarMatrix
  -> P.Sem r (CholeskySolutions, GLM.BetaU)
cholmodCholeskySolutionsLMM cholmodFactor mixedModelSpec zStar = do
  let mX = GLM.rmsFixedPredictors $ GLM.mmsRegressionModelSpec mixedModelSpec
      vY = GLM.rmsObservations $ GLM.mmsRegressionModelSpec mixedModelSpec
      levels = GLM.mmsFitSpecByGroup mixedModelSpec
--      smZS   = GLM.makeZS reCalc vTh
  cholmodCholeskySolutions' cholmodFactor
                            (SLA.transpose $ GLM.smZS zStar)
                            mX
                            NormalEquationsLMM
                            mixedModelSpec

observationsDistribution :: GLM.MixedModel b g -> GLM.ObservationsDistribution
observationsDistribution (GLM.LinearMixedModel _) = GLM.Normal
observationsDistribution (GLM.GeneralizedLinearMixedModel glmmSpec) =
  GLM.glmmsObservationsDistribution glmmSpec

useLink :: GLM.MixedModel b g -> GLM.UseLink
useLink (GLM.LinearMixedModel _) = GLM.UseCanonical
useLink (GLM.GeneralizedLinearMixedModel glmmSpec) =
  GLM.glmLink $ GLM.glmmsControls glmmSpec

linkFunctionType :: GLM.MixedModel b g -> GLM.LinkFunctionType
linkFunctionType mm = case useLink mm of
  GLM.UseOther x   -> x
  GLM.UseCanonical -> GLM.canonicalLink $ observationsDistribution mm

getCholeskySolutions
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModel b g
  -> GLM.ZStarMatrix
  -> OptimizationStage
  -> GLM.CovarianceVec
  -> P.Sem r (CholeskySolutions, GLM.BetaU)
getCholeskySolutions cf mm@(GLM.LinearMixedModel _) zStar _ vTh =
  P.wrapPrefix "getCholeskySolutions (LMM)" $ do
    P.logLE P.Diagnostic "Starting..."
    P.logLE P.Diagnostic $ "theta=" <> (T.pack $ show vTh)
    P.logLE P.Diagnostic "Calling cholmodCholeskySolutionsLMM"
    res <- cholmodCholeskySolutionsLMM cf (GLM.mixedModelSpec mm) zStar
    P.logLE P.Diagnostic "Finished"
    return res
-- TODO: There is much low hanging fruit here in terms of not re-calcing things
--
getCholeskySolutions cf mm@(GLM.GeneralizedLinearMixedModel glmmSpec) zStar os vTh
  = P.wrapPrefix "getCholeskySolutions (GLMM)" $ do
    P.logLE P.Diagnostic
      $  "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vTh)
      <> ")..."
    let
      mX          = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vY          = GLM.rmsObservations $ GLM.regressionModelSpec mm
      vW          = GLM.glmmsWeights glmmSpec
      od          = GLM.glmmsObservationsDistribution glmmSpec
      GLM.GLMMControls ul maxHalvings pirlsCC = GLM.glmmsControls glmmSpec
      n           = LA.size vY
      (_, q)      = SLA.dim (GLM.smZS zStar)
      (_, _, smP) = cf
      vEtaExact   = computeEta0 (linkFunctionType mm) vY
      svU0        = SD.toSparseVector $ LA.fromList $ L.replicate q 0 -- FIXME (maybe 0 is right here??)
      vBeta0      = betaFrom mX zStar svU0 vEtaExact -- approximate fit with no random effects
      vEta0       = GLM.denseLinearPredictor
        mX
        zStar
        (GLM.LP_ComputeFrom $ GLM.BetaU vBeta0 svU0)
{-        
    P.logLE P.Diagnostic
      $  "vEtaExact ="
      <> (T.pack $ show vEtaExact)
      <> "\nsvU0="
      <> (T.pack $ show svU0)
      <> "\nvBeta0="
      <> (T.pack $ show vBeta0)
      <> "\nvEta0="
      <> (T.pack $ show vEta0)
-}
    let
      lf = GLM.linkFunction $ linkFunctionType mm
      pdFunction os chol x y =
        fst <$> profiledDeviance' PDVNone ML mm zStar vTh os chol x y
      --incS betaU dBetaU x = addBetaU betaU (scaleBetaU x dBetaU)
      iterateUntilConverged n pdCurrent vEta svU = do
        P.logLE P.Diagnostic
          $  "iterateUntilConverged (n="
          <> (T.pack $ show n)
          <> "; pd="
          <> (T.pack $ show pdCurrent)
          <> ")"
        (pd', vEta', svU', dBetaU, chol, smU) <- updateEtaBetaU
          cf
          mm
          maxHalvings
          zStar
          (pdFunction os)
          vEta
          svU
        cc <- case GLM.pirlsConvergenceType pirlsCC of
          GLM.PCT_Eta -> do
            let dEta = (vEta - vEta') -- GLM.diffLinearPredictor mX zStar lp lp'
            return $ (LA.norm_2 dEta) / (LA.norm_2 vEta')
          GLM.PCT_Deviance   -> return $ abs $ (pdCurrent - pd') / pd'
          GLM.PCT_Orthogonal -> P.throw
            $ GLM.OtherGLMError "ConvergeOrthognal not implemented yet."
{-          
        occ <- pirlsConvergence glmm
                                smZS
                                smP
                                (lTheta chol)
                                smU
                                vEta
                                vEta'
                                svU'
                                (bu_svU dBetaU)
-}
        case (cc < (GLM.pirlsTolerance pirlsCC), n) of
          (True, _) ->
            P.logLE P.Diagnostic "Finished." >> return (vEta', svU', chol)
          (False, 1) -> P.throw
            (GLM.OtherGLMError "Too many iterations in getCholeskySolutions.")
          (False, m) -> do
            P.logLE P.Diagnostic $ "Not converged.  cc=" <> (T.pack $ show cc)
            iterateUntilConverged (m - 1) pd' vEta' svU'
    (vEta', svU', chol) <- iterateUntilConverged (GLM.pirlsMaxSteps pirlsCC)
                                                 (1.0 / 0.0) -- Infinity :: Double
                                                 vEta0
                                                 svU0
    let vBeta' = betaFrom mX zStar svU' vEta'
    P.logLE P.Diagnostic "Finished."
    return (chol, GLM.BetaU vBeta' svU')

updateEtaBetaU
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModel b g
  -> Int
  -> GLM.ZStarMatrix
  -> (  CholeskySolutions
     -> GLM.EtaVec
     -> GLM.UVec
     -> P.Sem r Double
     )
  -> GLM.EtaVec
  -> GLM.UVec
  -> P.Sem
       r
       ( Double
       , GLM.EtaVec
       , GLM.UVec
       , GLM.BetaU
       , CholeskySolutions
       , SLA.SpMatrix Double
       )
updateEtaBetaU cf mm maxHalvings zStar pdFunction vEta svU =
  P.wrapPrefix "updateEtaBetaU" $ do
    P.logLE P.Diagnostic $ "Starting..." ---(Eta=" <> (T.pack $ show vEta) <> ")"
    let mX = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
    (dBetaU, chol, smU)             <- compute_dBetaU cf mm zStar vEta svU
--    P.logLE P.Diagnostic $ "d" <> (T.pack $ show dBetaU)
    (pdFinal, vEta', svU', dBetaU') <- refine_dBetaU mm
                                                     maxHalvings
                                                     zStar
                                                     (pdFunction chol)
                                                     vEta
                                                     svU
                                                     dBetaU
    P.logLE P.Diagnostic
      $  "Finished (||dU||^2 = "
      <> (T.pack $ show $ (GLM.bu_svU dBetaU' SLA.<.> GLM.bu_svU dBetaU'))
      <> "; ||dBeta||^2="
      <> (T.pack $ show $ (GLM.bu_vBeta dBetaU' LA.<.> GLM.bu_vBeta dBetaU'))
      <> ")." --- (Eta=" <> (T.pack $ show vEta') <> ")"
    return (pdFinal, vEta', svU', dBetaU', chol, smU)

data ShrinkPD = Shrunk Double GLM.EtaVec GLM.UVec GLM.BetaU | NotShrunk deriving (Show)

refine_dBetaU
  :: GLM.EffectsIO r
  => GLM.MixedModel b g
  -> Int -- max step halvings
  -> GLM.ZStarMatrix
  -> (GLM.EtaVec -> GLM.UVec -> P.Sem r Double) -- profiled deviance
  -> GLM.EtaVec
  -> GLM.UVec
  -> GLM.BetaU
  -> P.Sem r (Double, GLM.EtaVec, GLM.UVec, GLM.BetaU)
refine_dBetaU mm maxHalvings zStar pdF vEta svU dBetaU =
  P.wrapPrefix "refineDBetaU" $ do
    P.logLE P.Diagnostic $ "Starting..."
    {-        
          <>  "\nvEta="
          <> (T.pack $ show vEta)
          <> "\nsvU'="
          <> (T.pack $ show svU')
          <> "\ndBetaU="
          <> (T.pack $ show dBetaU)
-}
    let --pdF x y = fst <$> profiledDeviance' PDVNone ML glmm reCalc vTh chol x y
        mX  = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
        vY  = GLM.rmsObservations $ GLM.regressionModelSpec mm
        tol = GLM.lmmOptimizerTolerance $ GLM.lmmControls mm
--        smZS = GLM.smZS zStar
    pd0 <- pdF vEta svU --vBeta svU
    let
      vdEta = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom dBetaU)
      checkOne x = do
        let deltaBetaU = GLM.scaleBetaU x dBetaU
            svU'       = svU SLA.^+^ (GLM.bu_svU deltaBetaU) --svBetaU' = incS svBetaU svdBetaU x
            --vBeta'     = vBeta + (GLM.bu_vBeta deltaBetaU)
            vEta'      = vEta + (LA.scale x vdEta)
        pdNew <- pdF vEta' svU'
        return $ if pdNew <= (pd0 + tol) -- sometimes dBetaU is basically 0 because we are at a minimum.  So we need a little breathing room.
          then Shrunk pdNew vEta' svU' deltaBetaU
          else NotShrunk
      check triesLeft x = do
        P.logLE P.Diagnostic
          $  "check (triesLeft="
          <> (T.pack $ show triesLeft)
          <> "; x="
          <> (T.pack $ show x)
          <> ")..."
        co <- checkOne x
        case co of
          NotShrunk -> if triesLeft > 1
            then
              (do
                check (triesLeft - 1) (x / 2)
              )
            else
              P.throw
              $  GLM.OtherGLMError
              $  "Too many step-halvings in getCholeskySolutions. (dBetaU = "
              <> (T.pack $ show dBetaU)
          sh@(Shrunk pdFinal vEtaFinal svUFinal deltaBetaU) -> do
--            P.logLE P.Diagnostic $ "refined: " <> (T.pack $ show sh)
            return (pdFinal, vEtaFinal, svUFinal, deltaBetaU)
    check maxHalvings 1.0


-- glmer: updateXwts but also with lambda
spUV
  :: GLM.MixedModel b g
  -> GLM.ZStarMatrix
  -> GLM.EtaVec
  -> (GLM.UMatrix, GLM.VMatrix)
spUV mm@(GLM.LinearMixedModel _) zStar _ =
  let mX = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
  in  (GLM.smZS zStar, mX)

spUV mm@(GLM.GeneralizedLinearMixedModel glmmSpec) zStar vEta =
  let mX    = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vW    = GLM.glmmsWeights glmmSpec
      smZS  = GLM.smZS zStar --GLM.makeZS reCalc vTh
      vWUV  = weightsForUV mm vW vEta
      smWUV = SLA.mkDiagonal (VS.length vWUV) $ VS.toList vWUV
      smU   = smWUV SLA.## smZS
      mV    = (LA.diag vWUV) LA.<> mX
  in  (smU, mV)

-- glmer: sqrtWrkWt
weightsForUV :: GLM.MixedModel b g -> GLM.WMatrix -> GLM.EtaVec -> GLM.WMatrix
weightsForUV mm vW vEta =
  let lf        = GLM.linkFunction $ linkFunctionType mm
      vMu       = VS.map (GLM.invLink lf) vEta
      vdMudEta  = VS.map (GLM.derivInv lf) vEta
      vVariance = GLM.scaledVariance (observationsDistribution mm) vMu
      vVSW      = GLM.varianceScaledWeights (observationsDistribution mm) vW vMu --    
  in  VS.zipWith (\muEta vsw -> muEta * sqrt vsw) vdMudEta vVSW
--    VS.zipWith3 (\w muEta var -> muEta * sqrt (w / var)) vW vdMudEta vVariance

compute_dBetaU
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModel b g
  -> GLM.ZStarMatrix
  -> GLM.EtaVec
  -> GLM.UVec
  -> P.Sem
       r
       (GLM.BetaU, CholeskySolutions, SLA.SpMatrix Double)
compute_dBetaU cf mm@(GLM.GeneralizedLinearMixedModel glmmSpec) zStar vEta svU
  = P.wrapPrefix "compute_dBetaU" $ do
    P.logLE P.Diagnostic "Starting..."
    let mX        = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
        smZS      = GLM.smZS zStar --GLM.makeZS reCalc vTh
        vW        = GLM.glmmsWeights glmmSpec
        lf        = GLM.linkFunction $ linkFunctionType mm
        (smU, mV) = spUV mm zStar vEta
        vMu       = VS.map (GLM.invLink lf) vEta
        vVarW = GLM.varianceScaledWeights (observationsDistribution mm) vW vMu
        neqs      = NormalEquationsGLMM vVarW lf vEta svU
    (chol, dBetaU) <- cholmodCholeskySolutions' cf
                                                (SLA.transpose smU)
                                                mV
                                                neqs
                                                (GLM.mixedModelSpec mm)
    P.logLE P.Diagnostic "Finished."
    return (dBetaU, chol, smU)

compute_dBetaU _ (GLM.LinearMixedModel _) _ _ _ = P.throw
  $ GLM.OtherGLMError "Called compute_dBetaU given an LMM instead of a GLMM."
{-
pirlsConvergence
  :: EffectsIO r
  => GeneralizedLinearMixedModel b g
  -> ZStarMatrix -- Z*Lambda(theta)
  -> PMatrix
  -> LMatrix
  -> UMatrix
  -> EtaVec
  -> EtaVec
  -> UVec
  -> UVec
  -> P.Sem r Double
pirlsConvergence (LMM _ _) _ _ _ _ _ _ _ _ =
  P.throw $ OtherGLMError "LMM given to pirlsConvergence"

pirlsConvergence glmm@(GLMM mm _ _ glmmC) smZS smP smL smU vEta vEta' svU svdU
  = P.wrapPrefix "pirlsConvergence" $ do
    let (GLMMControls _ _ _ pirlsC) = glmmC
--    P.logLE P.Diagnostic $ "U=" <> (T.pack $ show svU)
--    P.logLE P.Diagnostic $ "dU=" <> (T.pack $ show svdU)
    case pirlsC of
      ConvergeEta _ _ -> do
        let dEta = vEta' - vEta
--        P.logLE P.Diagnostic $ "Eta=" <> (T.pack $ show vEta')
--        P.logLE P.Diagnostic $ "dEta=" <> (T.pack $ show dEta)
        return $ (LA.norm_2 dEta) / (LA.norm_2 vEta')
      ConvergeOrthogonal _ _ -> do
        let MixedModel (RegressionModel _ mX vY) _ = mm
            lft        = linkFunctionType glmm
            vMu        = VS.map (GLM.invLink $ GLM.linkFunction lft) vEta
            svYMinusMu = SD.toSparseVector $ (vY - vMu)
            x1         = SLA.transpose smL SLA.#> svdU
            x2         = smL SLA.#> x1
            x3 = smP SLA.#> (SLA.transpose smU SLA.#> svYMinusMu) SLA.^-^ svU
            q'         = realToFrac $ snd $ SLA.dim smZS
            n'         = realToFrac $ LA.size vY
--        P.logLE P.Diagnostic $ "L'dU=" <> (T.pack $ show x1)
--        P.logLE P.Diagnostic $ "LL'dU=" <> (T.pack $ show x2)
--        P.logLE P.Diagnostic $ "(y - mu)" <> (T.pack $ show svYMinusMu)
--        P.logLE P.Diagnostic $ "P(smU')*(y - mu) - u=" <> (T.pack $ show x3)
        return $ (SLA.norm2 x1 / sqrt q') / (SLA.norm2 x2 / sqrt (n' - q'))
-}
-- compute a starting point for the linear predictor.  Just assume the given observations but fix any out of bounds
-- That is, assume mu = y* and then eta = g (mu) = g (y*) where
-- y* is y, but with some care at the boundaries.  That is, y but with all elements in the domain of g.
computeEta0 :: GLM.LinkFunctionType -> GLM.Observations -> GLM.EtaVec
computeEta0 lft vY =
  let lF = GLM.link $ GLM.linkFunction lft
      lowerBound lb x = if x < lb then lb else x
      upperBound ub x = if x > ub then ub else x
      bounded lb ub = lowerBound lb . upperBound ub
      f x = case lft of
        GLM.IdentityLink    -> x
        GLM.LogisticLink    -> bounded 0.00001 0.99999 x --minimum (maximum (0.0001, x), 0.99999)
        GLM.ExponentialLink -> lowerBound 0.0001 x
  in  VS.map (lF . f) vY

-- NB: eta = X <*> Beta + Z* <*> u
-- X <*> Beta = eta - Z* <*> u
-- X'X <*> Beta = X' <*> (eta - Z* <*> u)
-- Beta = inverse(X'X) <*> X' <*> (eta - Z* <*> u)
-- and X'X is symmetric positive definite so we solve the above via Cholesky
betaFrom
  :: GLM.FixedPredictors
  -> GLM.ZStarMatrix
  -> GLM.UVec
  -> GLM.EtaVec
  -> GLM.BetaVec
betaFrom mX zStar svU vEta =
  let vZSu = SD.toDenseVector $ (GLM.smZS zStar) SLA.#> svU
      vRhs = LA.tr mX LA.#> (vEta - vZSu)
      cXtX = LA.chol $ LA.trustSym $ LA.tr mX LA.<> mX
  in  head $ LA.toColumns $ LA.cholSolve cXtX (LA.asColumn vRhs)

cholmodCholeskySolutions'
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.UMatrix -- Ut (ZSt in linear case)
  -> GLM.VMatrix -- V (X in linear case)
  -> NormalEquations
  -> GLM.MixedModelSpec b g
  -> P.Sem r (CholeskySolutions, GLM.BetaU)
cholmodCholeskySolutions' cholmodFactor smUt mV nes mixedModelSpec =
  P.wrapPrefix "cholmodCholeskySolutions'" $ do
    let (cholmodC, cholmodF, smP) = cholmodFactor
        vY = GLM.rmsObservations $ GLM.mmsRegressionModelSpec mixedModelSpec
        n                         = LA.size vY
        (_, p)                    = LA.size mV
        (q, _)                    = SLA.dim smUt
    -- Cholesky factorize to get L_theta *and* update factor for solving with it
    P.logLE P.Diagnostic "Starting..."
    -- TODO: Cholmod has built in support for factorizing XtX + a * I.  Use it.
    let cfs = CH.FactorizeAtAPlusBetaI 1
    liftIO $ CH.spMatrixFactorizeP cholmodC cholmodF cfs CH.UnSymmetric smUt
    (svRhsZ, vRhsX) <- normalEquationsRHS nes smUt (LA.tr mV) vY
  {-  P.logLE P.Diagnostic
      $  "Normal Equation RHS:\nsmPRhsZ="
      <> (T.pack $ show smPRhsZ)
      <> "\nvRhsX="
      <> (T.pack $ show vRhsX)
    P.logLE P.Diagnostic "Calling solveSparse for svC."
-}
    smPRhsZ         <- liftIO
      $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P (SD.svColumnToSM svRhsZ)
--    P.logLE P.Diagnostic "After smPRhsZ..."
    svC <-
      SLA.toSV
        <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPRhsZ)
    let smUtV = smUt SLA.#~# (SD.toSparseMatrix mV)
--    P.logLE P.Diagnostic "Calling solveSparse for smPUtV."
    smPUtV <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smUtV -- PU'V
--    P.logLE P.Diagnostic "Calling solveSparse for smRzx."
    smRzx  <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPUtV -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
    -- compute Rxx
    -- NB: c is defined as Lu + Rzx(Beta) which allows solving the normal equations in pieces
    let vTvMinusRzxTRzx =
          (LA.tr mV) LA.<> mV - (SD.toDenseMatrix $ smRzx SLA.#^# smRzx)
        mRxx            = LA.chol $ LA.trustSym $ vTvMinusRzxTRzx
    -- now we have the entire Cholesky decomposition.  Solve the normal equations
        vRzxtC          = SD.toDenseVector $ (SLA.transposeSM smRzx) SLA.#> svC
        betaSols        = LA.cholSolve mRxx (LA.asColumn $ vRhsX - vRzxtC)
        vBeta           = head $ LA.toColumns $ betaSols
        svBeta          = SD.toSparseVector vBeta
        svRzxBeta       = smRzx SLA.#> svBeta
        smCMinusRzxBeta = SD.svColumnToSM (svC SLA.^-^ svRzxBeta)
--    P.logLE P.Diagnostic "Calling solveSparse for smPu."
    smPu <- liftIO
      $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Lt smCMinusRzxBeta
--    P.logLE P.Diagnostic "Calling solveSparse for svu."
    svu <-
      SLA.toSV
        <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Pt $ smPu)
    -- NB: This has to happen after the solves because it unfactors the factor...
    -- NB: It does *not* undo the analysis
--    P.logLE P.Diagnostic "Calling choleskyFactorSM for smLth."
    smLth <- liftIO $ CH.choleskyFactorSM cholmodF cholmodC
    P.logLE P.Diagnostic
      $  "min(u)="
      <> (T.pack $ show $ FL.fold FL.minimum svu)
      <> "; max(u)="
      <> (T.pack $ show $ FL.fold FL.maximum svu)
      <> "; min(beta)="
      <> (T.pack $ show $ VS.minimum vBeta)
      <> "; max(beta)="
      <> (T.pack $ show $ VS.maximum vBeta)
    P.logLE P.Diagnostic "Finished."
    return $ (CholeskySolutions smLth smRzx mRxx, GLM.BetaU vBeta svu)



---- Unused
{-
cholmodCholeskySolutions_Old
  :: CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> IO (CholeskySolutions, BetaU)
cholmodCholeskySolutions_Old cholmodFactor mixedModel randomEffCalcs vTh = do
  let (cholmodC, cholmodF, _) = cholmodFactor
      (MixedModel             (RegressionModel _ mX vY) levels  ) = mixedModel
      (RandomEffectCalculated smZ mkLambda) = randomEffCalcs
      lambda                  = mkLambda vTh
      smZS                    = smZ SLA.## lambda
      smZSt                   = SLA.transpose smZS
      n                       = LA.size vY
      (_, p)                  = LA.size mX
      (_, q)                  = SLA.dim smZ
  -- Cholesky factorize to get L_theta *and* update factor for solving with it
  CH.spMatrixFactorize cholmodC cholmodF CH.SquareSymmetricLower
    $ SLA.filterSM lowerTriangular
    $ xTxPlusI
    $ smZS
    -- compute Rzx

  let svZSty = smZSt SLA.#> SD.toSparseVector vY
  smPZSty <- CH.solveSparse cholmodC
                            cholmodF
                            CH.CHOLMOD_P
                            (SD.svColumnToSM svZSty)
  svCu <- SLA.toSV <$> CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPZSty
  let smZStX = smZSt SLA.#~# (SD.toSparseMatrix mX)
  smPZStX <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smZStX -- P(Z*)'X
  smRzx   <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPZStX -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
  --  let smRxTRx = (SD.toSparseMatrix $ LA.tr mX LA.<> mX) SLA.^-^ (SLA.transposeSM smRzx SLA.## smRzx)
  --      beta = SLA.luSolve (SD.toSparseVector (mX LA.#> vY) SLA.^-^ (smRzx SLA.#> svCu)
  -- compute Rx
-- NB: cu is defined as Lu + Rzx(Beta) which allows solving the equations in pieces
  let xTxMinusRzxTRzx =
        (LA.tr mX) LA.<> mX - (SD.toDenseMatrix $ smRzx SLA.#^# smRzx)
      mRx              = LA.chol $ LA.trustSym $ xTxMinusRzxTRzx
      vXty             = (LA.tr mX) LA.#> vY
      vRzxtCu          = SD.toDenseVector $ (SLA.transposeSM smRzx) SLA.#> svCu
      vXtyMinusRzxtCu  = vXty - vRzxtCu
      betaSols         = LA.cholSolve mRx (LA.asColumn vXtyMinusRzxtCu)
      vBeta            = head $ LA.toColumns $ betaSols
      svBeta           = SD.toSparseVector vBeta
      svRzxBeta        = smRzx SLA.#> svBeta
--      smRzxBeta        = SD.svColumnToSM svRzxBeta
      smCuMinusRzxBeta = SD.svColumnToSM (svCu SLA.^-^ svRzxBeta)
  smPu  <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Lt smCuMinusRzxBeta
  svu   <- SLA.toSV <$> (CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Pt $ smPu)
  -- NB: This has to happen after the solves because it unfactors the factor...
  -- NB: It does *not* undo the analysis
  smLth <- CH.choleskyFactorSM cholmodF cholmodC
  return (CholeskySolutions smLth smRzx mRx, BetaU svBeta svu)
-}
