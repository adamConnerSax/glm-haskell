{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE DeriveDataTypeable  #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators #-}
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
import qualified Polysemy.State                as P
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
import           Data.Maybe                     (fromMaybe)
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA
--import           Data.Typeable                  ( Typeable )

--import           Numeric.LinearAlgebra.Sparse   ( (#>) )

import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL
import           System.IO.Unsafe               ( unsafePerformIO )



import qualified Data.Map                      as M
--import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

import qualified Say 

import           System.IO                      ( hSetBuffering
                                                , stdout
                                                , hFlush
                                                , BufferMode(..)
                                                )

--import qualified Debug.Trace as Trace

runEffectsVerboseIO :: GLM.GLMEffects a -> IO (Either GLM.GLMError a)
runEffectsVerboseIO action =
  (P.catch action glmExceptionLogger)
  & runLogOnGLMException
  & P.filteredLogEntriesToIO P.logAll
  & P.errorToIOFinal
  & P.embedToFinal
  & P.runFinal


runEffectsIO :: GLM.GLMEffects a -> IO (Either GLM.GLMError a)
runEffectsIO action =
  (P.catch action glmExceptionLogger)
  & runLogOnGLMException
  & P.filteredLogEntriesToIO P.nonDiagnostic
  & P.errorToIOFinal
  & P.embedToFinal
  & P.runFinal


unsafePerformEffects :: Bool -> GLM.GLMEffects a -> a -> a
unsafePerformEffects verbose action def = 
  let run = if verbose then runEffectsVerboseIO else runEffectsIO 
  in (P.wrapPrefix "unsafePerformEffects" action)
     & run
     & fmap (either X.throwIO return)
     & join
     & (\x -> X.catch
         x
         (\(e :: X.SomeException) -> Say.say ("Error (during unsafePerformEffects in objective): " <> (T.pack $ show e)) >> return def
         )
       )
     & unsafePerformIO


logOnGLMException :: P.Member (P.State [T.Text]) r => T.Text -> P.Sem r ()
logOnGLMException t = P.modify (\msgs -> t : msgs)

runLogOnGLMException :: P.Sem (P.State [T.Text] ': r) a -> P.Sem r a
runLogOnGLMException = fmap snd . P.runState []

glmExceptionLogToText :: [T.Text] -> T.Text
glmExceptionLogToText =  T.intercalate "\n" . reverse

getGLMExceptionLog :: P.Member (P.State [T.Text]) r => P.Sem r T.Text
getGLMExceptionLog = glmExceptionLogToText <$> P.get

glmExceptionLogger :: GLM.Effects r => GLM.GLMError -> P.Sem r a
glmExceptionLogger e = do
   l <- getGLMExceptionLog
   P.logLE P.Info $ "GLM Error: " <> (T.pack $ show e)
   P.logLE P.Diagnostic $ "Exception Log:\n" <> l
   P.throw e

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
lmmOptimizerToNLOPT GLM.LMM_COBYLA = \obj bds mIn -> NL.COBYLA obj bds [] [] mIn

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

data MinimizeDevianceVerbosity = MDVNone | MDVSimple deriving (Eq, Ord)

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
         , GLM.BetaVec
         , GLM.MaybeZeroUVec
         , GLM.MaybeZeroVec (LA.Vector Double)
         , CholeskySolutions
         )
       , LA.Vector Double
       , CholmodFactor
       ) -- ^ ((theta, profiled_deviance, sigma2, beta, Maybe u, b, Cholesky blocks), mu, Cholmod Analysis)
minimizeDeviance verbosity dt mm reCalc th0 = 
  flip P.catch (\(e :: GLM.GLMError) -> do
                     l <- getGLMExceptionLog
                     P.logLE P.Diagnostic $ "Exception Log:\n" <> l
                     P.throw e
                 )
  $ P.wrapPrefix "minimizeDeviance" $ do
      P.logLE P.Diagnostic "Starting..."
      when (GLM.generalized mm && (dt == REML)) $ P.throw $ GLM.OtherGLMError
        "Can't use REML for generalized models."
      P.logLE P.Diagnostic "Beginning Cholmod analysis of Z matrix structure."
      cholmodAnalysis <- cholmodAnalyzeProblem reCalc
      P.logLE P.Diagnostic "Finished Cholmod analysis of Z matrix structure."
      res@(vTh, _, dev2, vBeta, mzvu, _, _) <- minimizeDevianceInner verbosity
                                               dt
                                               mm
                                               reCalc
                                               cholmodAnalysis
                                               th0
      let
        mX    = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
        zStar = GLM.makeZS reCalc vTh
        vEta  = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom vBeta mzvu)
        vMu = LA.cmap (GLM.invLink $ GLM.linkFunction (linkFunctionType mm)) vEta
      return (res, vMu, cholmodAnalysis)



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
       , GLM.BetaVec
       , GLM.MaybeZeroUVec
       , GLM.MaybeZeroVec (LA.Vector Double)
       , CholeskySolutions
       ) -- ^ (theta, profiled_deviance, sigma2, beta, u, b, cholesky blocks)
minimizeDevianceInner verbosity dt mm reCalc cf th0 =
  P.wrapPrefix "minimizeDevianceInner" $ do
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
             , GLM.BetaVec
             , GLM.MaybeZeroUVec
             , GLM.MaybeZeroBVec
             , CholeskySolutions
             )
      pd os x = profiledDeviance pdv cf dt mm reCalc os (Just x)
      obj x = unsafePerformEffects
        (verbosity >= MDVSimple)
        (fmap (\(d, _, _, _, _, _) -> d) $ pd objectiveOs x)
        0
      levels    = GLM.mmsFitSpecByGroup $ GLM.mixedModelSpec mm
      thetaLB   = thetaLowerBounds levels
      algorithm = (lmmOptimizerToNLOPT $ GLM.lmmOptimizer $ GLM.lmmControls mm)
        obj
        [thetaLB]
        (NL.InitialStep <$> (GLM.lmmOptimizerInitialStep $ GLM.lmmControls mm))
      stop =
        NL.ObjectiveRelativeTolerance
        (GLM.lmmOptimizerTolerance $ GLM.lmmControls mm)
        NL.:|  [NL.ObjectiveAbsoluteTolerance (GLM.lmmOptimizerTolerance $ GLM.lmmControls mm)]
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
        (pdVal, sigma2, vBeta, mzvu, mzvb, cs) <- pd finalOs thS
        return (thS, pdVal, sigma2, vBeta, mzvu, fmap SD.toDenseVector mzvb, cs)


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

zCols :: GLM.FitSpecByGroup g -> Int
zCols = FL.fold FL.sum . fmap colsForGroup  
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
    q        = zCols groupFSM --FL.fold FL.sum $ fmap colsForGroup groupFSM -- total number of columns in Z
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
--  liftIO $ CH.printCommon "before analysis" cholmodC 
  (cholmodF, smP) <-
    liftIO
    $ CH.spMatrixAnalyzeWP cholmodC CH.UnSymmetric
--    $ xTx
    $ SLA.transpose
    $ GLM.recModelMatrix reCalc
--  liftIO $ CH.printFactor cholmodF "After analyze, before factorize" cholmodC
--  liftIO $ putStrLn "smP=" >> (LA.disp 1 $ SD.toDenseMatrix smP)
--  liftIO $ CH.printCommon "after analysis" cholmodC
--  liftIO $ CH.printFactor cholmodF  "after analysis" cholmodC
  return (cholmodC, cholmodF, smP)

data ProfiledDevianceVerbosity = PDVNone | PDVSimple | PDVAll deriving (Show, Eq)
data OptimizationStage = Optim_LMM | Optim_GLMM_PIRLS | Optim_GLMM_Final deriving (Show, Eq)

profiledDeviance
  :: GLM.EffectsIO r
  => ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> GLM.MixedModel b g
  -> GLM.RandomEffectCalculated
  -> OptimizationStage
  -> Maybe GLM.CovarianceVec -- ^ theta
  -> P.Sem
       r
       ( Double
       , Double
       , GLM.BetaVec
       , GLM.MaybeZeroUVec
       , GLM.MaybeZeroBVec
       , CholeskySolutions
       ) -- ^ (pd, sigma^2, beta, u, b) 
profiledDeviance verbosity cf dt mm reCalc os vThM =
  P.wrapPrefix "profiledDeviance" $ do
--    P.logLE P.Diagnostic $ "here"
    let mX     = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
        smZ    = GLM.recModelMatrix reCalc
        vY     = GLM.rmsObservations $ GLM.regressionModelSpec mm
        n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        vTh = fromMaybe (setCovarianceVector (GLM.mmsFitSpecByGroup $ GLM.mixedModelSpec mm) 0 0) vThM 
        lambda = GLM.recMkLambda reCalc vTh
        zStar = GLM.makeZS reCalc vTh

    P.logLE P.Diagnostic "Starting. Calling getCholeskySolutions."
    (cs@(CholeskySolutions smLth smRzx mRxx), vBeta, mzvu) <- getCholeskySolutions
      cf
      mm
      zStar
      os
      vThM
    when (verbosity == PDVAll) $ liftIO $ do -- needs to be in IO for some of the show functions :(
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
        vEta    = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom vBeta mzvu)
    (pd, rTheta2) <- profiledDeviance' verbosity
                                       dt
                                       mm
                                       zStar
                                       vThM
                                       osFinal
                                       cs
                                       vEta                                       
                                       mzvu
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
    let mzvb = fmap (\svU -> lambda SLA.#> svU) mzvu
    return (pd, sigma2, vBeta, mzvb, mzvb, cs)


profiledDeviance'
  :: GLM.EffectsIO r
  => ProfiledDevianceVerbosity
  -> DevianceType
  -> GLM.MixedModel b g
  -> GLM.ZStarMatrix
  -> Maybe GLM.CovarianceVec
  -> OptimizationStage
  -> CholeskySolutions
  -> GLM.EtaVec
  -> GLM.MaybeZeroUVec
  -> P.Sem r (Double, Double) -- profiled deviance, deviance + u'u 
profiledDeviance' pdv dt mm zStar vThM os chol vEta mzvu =
  P.wrapPrefix "profiledDeviance'" $ do
    P.logLE P.Diagnostic
      $  "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vThM)
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
      uDotu = GLM.useMaybeZeroVec 0 (\svU -> svU SLA.<.> svU) mzvu
    P.logLE P.Diagnostic
      $  "logLth="
      <> (T.pack $ show logLth)
      <> "; u.u="
      <> (T.pack $ show uDotu)
      <> "; min mu="
      <> (T.pack $ show $ VS.minimum vMu)
      <> "; max mu="
      <> (T.pack $ show $ VS.maximum vMu)
    (pd, rTheta2) <- case mm of
      GLM.LinearMixedModel _ -> do
        when (os /= Optim_LMM)
          $ P.throw
          $ GLM.OtherGLMError
              "In profiledDeviance': OptimizationStage must be Optim_LMM for all LMM calls"
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

data NormalEquations = NormalEquationsLMM | NormalEquationsGLMM GLM.WMatrix (VS.Vector Double) GLM.LinkFunction GLM.EtaVec GLM.MuVec GLM.MaybeZeroUVec

normalEquationsRHS
  :: GLM.EffectsIO r
  => NormalEquations
  -> Maybe GLM.UMatrix
  -> GLM.VMatrix
  -> GLM.Observations
  -> P.Sem r (GLM.MaybeZeroVec (SLA.SpVector Double), LA.Vector Double)
normalEquationsRHS NormalEquationsLMM smUtM mVt vY =
  P.wrapPrefix "normalEquationsRHS" $ do
    let mzvRhsZ = GLM.MaybeZeroVec $ fmap (\smUt -> smUt SLA.#> SD.toSparseVector vY) smUtM 
        vRhsX  = mVt LA.#> vY
--    P.logLE P.Diagnostic $ "RHS=" <> (T.pack $ show (svRhsZ, vRhsX))
    return (mzvRhsZ, vRhsX) --(SLA.transpose smU SLA.#> SD.toSparseVector vY, LA.tr mV LA.#> vY)

normalEquationsRHS (NormalEquationsGLMM vVSW vM lf vEta vMu mzvu) smUtM mVt vY =
  P.wrapPrefix "normalEquationsRHS" $ do
    P.logLE P.Diagnostic "Starting..."

    let vYMinusMu   = vY - vMu        
        vWtYMinusMu = VS.zipWith3 (\wt y mu -> sqrt wt * (y - mu)) vVSW vY vMu
--        v1 = VS.zipWith (\m yMinusMu -> yMinusMu / m) vM vYMinusMu
--        v2 = v1 + vEta
--        vR0 = VS.zipWith3 (\w m x -> sqrt w * m * x) vVarW vM v2
--        vR0 = VS.zipWith4 (\w m ym eta -> sqrt w * (ym + (m * eta))) vVarW vM vYMinusMu vEta -- for equations for u, beta
        vR0 = VS.zipWith (\w ym -> sqrt w * ym) vVSW vYMinusMu -- for equations for du, dBeta
        msvUt_r0 = GLM.MaybeZeroVec $ fmap (\smUt -> smUt SLA.#> (SD.toSparseVector vR0)) smUtM
        vVt_r0 = mVt LA.#> vR0
    let msvUX   = GLM.MaybeZeroVec $ fmap (\smUt -> smUt SLA.#> (SD.toSparseVector vWtYMinusMu)) smUtM
        g v = GLM.useMaybeZeroVec v (\x -> v SLA.^-^ x)
        msvRhsZ = fmap (\svUX -> g svUX mzvu) msvUX
        vRhsX  = mVt LA.#> vWtYMinusMu
--    P.logLE P.Diagnostic $ "r0=" <> (T.pack $ show vR0)
--    P.logLE P.Diagnostic $ "RHS=" <> (T.pack $ show (svUt_r0, vVt_r0))
--    return (msvUt_r0, vVt_r0)
    return (msvRhsZ, vRhsX)

cholmodCholeskySolutionsLMM
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModelSpec b g
  -> Maybe GLM.ZStarMatrix
  -> P.Sem r (CholeskySolutions, GLM.BetaVec, GLM.MaybeZeroUVec)
cholmodCholeskySolutionsLMM cholmodFactor mixedModelSpec zStarM = do
  let mX = GLM.rmsFixedPredictors $ GLM.mmsRegressionModelSpec mixedModelSpec
      vY = GLM.rmsObservations $ GLM.mmsRegressionModelSpec mixedModelSpec
      levels = GLM.mmsFitSpecByGroup mixedModelSpec
--      smZS   = GLM.makeZS reCalc vTh
  cholmodCholeskySolutions' cholmodFactor
                            (fmap (SLA.transpose . GLM.smZS) zStarM)
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
  -> Maybe GLM.CovarianceVec
  -> P.Sem r (CholeskySolutions, GLM.BetaVec, GLM.MaybeZeroUVec)
getCholeskySolutions cf mm@(GLM.LinearMixedModel _) zStar _ vThM =
  P.wrapPrefix "getCholeskySolutions (LMM)" $ do
    P.logLE P.Diagnostic "Starting..."
    P.logLE P.Diagnostic $ "theta=" <> (T.pack $ show vThM)
    P.logLE P.Diagnostic "Calling cholmodCholeskySolutionsLMM"
    res <- cholmodCholeskySolutionsLMM cf (GLM.mixedModelSpec mm) (fmap (const zStar) vThM)
    P.logLE P.Diagnostic "Finished"
    return res
-- TODO: There is much low hanging fruit here in terms of not re-calcing things
--
getCholeskySolutions cf mm@(GLM.GeneralizedLinearMixedModel glmmSpec) zStar os vThM
  = P.wrapPrefix "getCholeskySolutions (GLMM)" $ do
    P.logLE P.Diagnostic
      $  "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vThM)
      <> ")..."
    logOnGLMException $
      "Starting ("
      <> (T.pack $ show os)
      <> " @ theta="
      <> (T.pack $ show vThM)
      <> ")..."
    let
      mX           = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vY           = GLM.rmsObservations $ GLM.regressionModelSpec mm
      vW           = GLM.weights mm
      od           = GLM.glmmsObservationsDistribution glmmSpec
      GLM.GLMMControls ul maxHalvings pirlsCC = GLM.glmmsControls glmmSpec
      n            = LA.size vY
      (_, q)       = SLA.dim (GLM.smZS zStar)
      (_, _, smP)  = cf
      vEtaExact    = computeEta0 (linkFunctionType mm) vY
    vBetaAtZeroU <- betaFrom' mm zStar vEtaExact GLM.zeroVec -- approximate fit with no random effects
--      svU0        = SD.toSparseVector $ LA.fromList $ L.replicate q 0 -- FIXME (maybe 0 is right here??)
    let lf = GLM.linkFunction $ linkFunctionType mm
        pdFunction os chol x y =
          fst <$> profiledDeviance' PDVNone ML mm zStar vThM os chol x y
    (vEta0, vBeta0, mzvu0, _) <- pirls cf mm BetaAtZeroU zStar (GLM.pirlsMaxSteps pirlsCC) (1.0 / 0.0) (pdFunction os) vEtaExact vBetaAtZeroU GLM.zeroVec
{-
      vEta0       = GLM.denseLinearPredictor
        mX
        zStar
        (GLM.LP_ComputeFrom $ GLM.BetaU vBeta0 svU0)
-}        
    P.logLE P.Diagnostic
      $  "min(vEtaExact)="
      <> (T.pack $ show $ VS.minimum vEtaExact)
      <>  "; max(vEtaExact)="
      <> (T.pack $ show $ VS.maximum vEtaExact)
      <> "\nmin(svU0)="
      <> (T.pack $ show $ fmap (VB.minimum . SLA.toVector) mzvu0)
      <> "; max(svU0)="
      <> (T.pack $ show $ fmap (VB.maximum . SLA.toVector) mzvu0)
      <> "\nmin(vBeta0)="
      <> (T.pack $ show $ VS.minimum vBeta0)
      <> "; max(vBeta0)="
      <> (T.pack $ show $ VS.maximum vBeta0)
      <> "\nmin(vEta0)="
      <> (T.pack $ show $ VS.minimum vEta0)
      <> "; max(vEta0)="
      <> (T.pack $ show $ VS.maximum vEta0)



    
      --incS betaU dBetaU x = addBetaU betaU (scaleBetaU x dBetaU)
    
    (vEta', vBeta', mzvu', chol) <- pirls
                                   cf
                                   mm
                                   SolveBetaU
                                   zStar
                                   (GLM.pirlsMaxSteps pirlsCC) 
                                   (1.0 / 0.0) -- Infinity :: Double
                                   (pdFunction os)
                                   vEta0
                                   vBeta0
                                   mzvu0
--    let vBeta' = betaFrom mX zStar svU' vEta'
    P.logLE P.Diagnostic "Finished."
    return (chol, vBeta', mzvu')


data PIRLS_Type = SolveBetaU | BetaAtZeroU

printPIRLSType :: PIRLS_Type -> T.Text
printPIRLSType SolveBetaU = "Solving For Beta and u"
printPIRLSType BetaAtZeroU = "Solving For Beta at fixed u=0"


pirls ::  GLM.EffectsIO r
      => CholmodFactor
      -> GLM.MixedModel b g
      -> PIRLS_Type --GLM.ZStarMatrix -- ^ if this is nothing, we are only solving for beta
      -> GLM.ZStarMatrix
      -> Int
      -> Double 
      -> (CholeskySolutions
          -> GLM.EtaVec
          -> GLM.MaybeZeroUVec
          -> P.Sem r Double
         )
      -> GLM.EtaVec
      -> GLM.BetaVec
      -> GLM.MaybeZeroUVec
      -> P.Sem r (GLM.EtaVec, GLM.BetaVec, GLM.MaybeZeroUVec,  CholeskySolutions)

pirls _ mm@(GLM.LinearMixedModel _) _ _ _ _ _ _ _ _ = P.throw $ GLM.OtherGLMError "Called pirls with a non-generalized mixed model."

pirls cf mm@(GLM.GeneralizedLinearMixedModel glmmSpec) pirlsType zStar n pdCurrent pdF vEta vBeta mzvu =
  P.wrapPrefix ("PIRLS (" <> printPIRLSType pirlsType <> ")") $ do
    P.logLE P.Diagnostic
      $  "iterateUntilConverged (n="
        <> (T.pack $ show n)
        <> "; pd="
        <> (T.pack $ show pdCurrent)
    P.logLE P.Diagnostic
      $ "min(vEta)="
        <> (T.pack $ show $ VS.minimum vEta) 
        <> "; max(vEta)="
        <> (T.pack $ show $ VS.maximum vEta)
    P.logLE P.Diagnostic
      $ "min(svU)="
      <> (T.pack $ show $ fmap (VB.minimum . SLA.toVector) mzvu) 
      <> "; max(svU)="
      <> (T.pack $ show $ fmap (VB.maximum . SLA.toVector) mzvu)
      <> ")"
    let GLM.GLMMControls ul maxHalvings pirlsCC = GLM.glmmsControls glmmSpec
    (pd', vEta', vBeta', mzvu', chol, smU) <- updateEtaBetaU
                                                      cf
                                                      mm
                                                      maxHalvings
                                                      pirlsType
                                                      zStar
                                                      pdF
                                                      vEta
                                                      vBeta
                                                      mzvu

    cc <- case GLM.pirlsConvergenceType pirlsCC of
      GLM.PCT_Eta -> do
        let dEta = (vEta - vEta') -- GLM.diffLinearPredictor mX zStar lp lp'
        return $ (LA.norm_2 dEta) / (LA.norm_2 vEta')
      GLM.PCT_Deviance   -> return $ abs $ (pdCurrent - pd') / pd'
      GLM.PCT_Orthogonal -> P.throw $ GLM.OtherGLMError "ConvergeOrthognal not implemented yet."

    case (cc < (GLM.pirlsTolerance pirlsCC), n) of
      (True, _) -> do
        P.logLE P.Diagnostic "Converged."
        return (vEta', vBeta', mzvu', chol)
      (False, 1) -> P.throw $ GLM.OtherGLMError "Too many iterations in getCholeskySolutions."
      (False, m) -> do
        P.logLE P.Diagnostic $ "Not converged.  cc=" <> (T.pack $ show cc)
        pirls cf mm pirlsType zStar (m - 1) pd' pdF vEta' vBeta' mzvu'

updateEtaBetaU
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModel b g
  -> Int
  -> PIRLS_Type
  -> GLM.ZStarMatrix
  -> (  CholeskySolutions
     -> GLM.EtaVec
     -> GLM.MaybeZeroUVec
     -> P.Sem r Double
     )
  -> GLM.EtaVec
  -> GLM.BetaVec
  -> GLM.MaybeZeroUVec
  -> P.Sem
       r
       ( Double
       , GLM.EtaVec
       , GLM.BetaVec
       , GLM.MaybeZeroUVec
       , CholeskySolutions
       , Maybe GLM.UMatrix
       )
updateEtaBetaU cf mm maxHalvings pirlsType zStar pdFunction vEta vBeta mzvu =
  P.wrapPrefix "updateEtaBetaU" $ do
    P.logLE P.Diagnostic $ "Starting..." ---(Eta=" <> (T.pack $ show vEta) <> ")"
    let mX = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
--        vBeta = betaFrom mX zStar svU vEta        
--        betaU0 = (GLM.BetaU vBeta svU)
    P.logLE P.Diagnostic
      $ "min(vBeta)="
      <> (T.pack $ show $ VS.minimum vBeta) 
      <> "; max(vBeta)="
      <> (T.pack $ show $ VS.maximum vBeta)
    (vdBeta, mzvdu, chol, msmU, mV)             <- compute_dBetaU cf mm pirlsType zStar vEta mzvu
--    let dBetaU = GLM.diffBetaU betaU betaU0
    logOnGLMException $ "updateEtaBetaU: "
      <> "beta0="
      <> (T.pack $ show vBeta)
      <> "\nd(beta)="
      <> (T.pack $ show vdBeta )
      <> "\nEta="
      <> (T.pack $ show vEta)
      <> "\nDist="
      <> (T.pack $ show  $ observationsDistribution mm)
      <> "\nRxx="
      <> (T.pack $ show $ rXX chol)
      <> "\nVtV="
      <> (T.pack $ show (LA.tr mV LA.<> mV))
    (pdFinal, vEta', vdBeta', mzvdu') <- refine_dBetaU mm
                                         maxHalvings
                                         pirlsType
                                         zStar
                                         (pdFunction chol)
                                         vEta
                                         mzvu
                                         vdBeta
                                         mzvdu
    let vBeta' = vBeta + vdBeta'
        mzvu' = mzvu `GLM.addMZU` mzvdu' 
    P.logLE P.Diagnostic
      $  "Finished (||dU||^2 = "
      <> (T.pack $ show $ fmap (\svu -> svu  SLA.<.> svu) mzvdu')
      <> "; ||dBeta||^2="
      <> (T.pack $ show $ (vdBeta' LA.<.> vdBeta'))
      <> ")." --- (Eta=" <> (T.pack $ show vEta') <> ")"
    return (pdFinal, vEta', vBeta', mzvu', chol, msmU)

data ShrinkPD = Shrunk Double GLM.EtaVec GLM.BetaVec GLM.MaybeZeroUVec  | NotShrunk Double deriving (Show)

refine_dBetaU
  :: GLM.EffectsIO r
  => GLM.MixedModel b g
  -> Int -- max step halvings
  -> PIRLS_Type
  -> GLM.ZStarMatrix
  -> (GLM.EtaVec -> GLM.MaybeZeroUVec -> P.Sem r Double) -- profiled deviance
  -> GLM.EtaVec
  -> GLM.MaybeZeroUVec
  -> GLM.BetaVec
  -> GLM.MaybeZeroUVec
  -> P.Sem r (Double, GLM.EtaVec, GLM.BetaVec, GLM.MaybeZeroUVec)
refine_dBetaU mm maxHalvings pirlsType zStar pdF vEta mzvu vdBeta mzvdu =
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
    pd0 <- pdF vEta mzvu --vBeta svU
    let
      
      vdEta = GLM.denseLinearPredictor mX zStar (GLM.LP_ComputeFrom vdBeta mzvdu)
      checkOne x = do
        let mzvdu' = fmap (SLA.scale x) mzvdu
            mzvu'  = mzvu `GLM.addMZU` mzvdu'
            vdBeta'     = LA.scale x vdBeta
            vEta'      = vEta + (LA.scale x vdEta)
        pdNew <- pdF vEta' mzvu'
        return $ if (pdNew - pd0) / pd0 <= tol -- sometimes dBetaU is basically 0 because we are at a minimum.  So we need a little breathing room.
          then Shrunk pdNew vEta' vdBeta' mzvdu' 
          else NotShrunk pdNew
      check triesLeft x xs = do
        P.logLE P.Diagnostic
          $  "check (triesLeft="
          <> (T.pack $ show triesLeft)
          <> "; x="
          <> (T.pack $ show x)
          <> ")..."
        co <- checkOne x
        case co of
          NotShrunk pd -> if triesLeft > 1
                          then check (triesLeft - 1) (x / 2) (pd : xs)
                          else (do
                                   let dBetaNorm2 = LA.norm_2 $ vdBeta
                                       dUNorm2 = fmap SLA.norm2 mzvdu
                                   P.throw
                                     $  GLM.OtherGLMError
                                     $  ("Too many step-halvings in getCholeskySolutions: Norm2(dBeta) = "
                                         <> (T.pack $ show dBetaNorm2)
                                         <> "; Norm2(dU) = "
                                         <> (T.pack $ show dUNorm2)
                                         <> "; pd0="
                                         <> (T.pack $ show pd0)
                                         <> "; pdPath="
                                         <> (T.pack $ show $ reverse (pd : xs))                                         
                                        )
                               )
          sh@(Shrunk pdFinal vEtaFinal svUFinal deltaBetaU) -> do
--            P.logLE P.Diagnostic $ "refined: " <> (T.pack $ show sh)
            return (pdFinal, vEtaFinal, svUFinal, deltaBetaU)
    check maxHalvings 1.0 []


-- glmer: updateXwts but also with lambda, with mu tacked on to be handed up the chain
spUV
  :: GLM.EffectsIO r
  => GLM.MixedModel b g
  -> PIRLS_Type
  -> GLM.ZStarMatrix
  -> GLM.WMatrix
  -> VS.Vector Double
  -> P.Sem r (Maybe GLM.UMatrix, GLM.VMatrix)
spUV mm@(GLM.LinearMixedModel _) pirlsType zStar _ _ = P.wrapPrefix "spUV" $ do
  let mX = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
  return (Just $ GLM.smZS zStar, mX)

spUV mm@(GLM.GeneralizedLinearMixedModel _) pirlsType zStar vVSW vM = P.wrapPrefix "spUV" $ do
  let mX    = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vWUV  = VS.zipWith (\m vsw -> sqrt vsw * m) vM vVSW -- glmer: sqrtWrkWt ?  
      smWUV = SLA.mkDiagonal (VS.length vWUV) $ VS.toList vWUV
      mV    = (LA.diag vWUV) LA.<> mX
      msmU  = case pirlsType of
        SolveBetaU -> Just $ smWUV SLA.## (GLM.smZS zStar)
        BetaAtZeroU -> Nothing
  P.logLE P.Diagnostic $ "min(vWUV)=" <> (T.pack $ show $ VS.minimum vWUV ) <> "; max(vWUV)=" <> (T.pack $ show $ VS.maximum vWUV )
  return (msmU, mV)

{-
-- glmer: sqrtWrkWt
weightsForUV :: GLM.MixedModel b g -> GLM.WMatrix -> VS.Vector Double -> GLM.EtaVec -> GLM.MuVec -> GLM.WMatrix
weightsForUV mm vVSW vM vEta vMu = case useLink mm of
  GLM.UseCanonical -> GLM.canonicaldEtadMuOverVar (observationsDistribution mm) vEta
  GLM.UseOther _ ->     
    let lf       = GLM.linkFunction $ linkFunctionType mm
--
--        vdMudEta = VS.map (GLM.derivInv lf) vEta
  --      vVariance = GLM.scaledVariance (observationsDistribution mm) vMu
--        vVSW     = GLM.varianceScaledWeights (observationsDistribution mm) vW vMu --    
    in  VS.zipWith (\muEta vsw -> muEta * sqrt vsw) vdMudEta vVSW
--    VS.zipWith3 (\w muEta var -> muEta * sqrt (w / var)) vW vdMudEta vVariance
-}

compute_dBetaU
  :: GLM.EffectsIO r
  => CholmodFactor
  -> GLM.MixedModel b g
  -> PIRLS_Type
  -> GLM.ZStarMatrix
  -> GLM.EtaVec
  -> GLM.MaybeZeroUVec
  -> P.Sem
       r
       (GLM.BetaVec, GLM.MaybeZeroUVec, CholeskySolutions, Maybe GLM.UMatrix, GLM.VMatrix)
compute_dBetaU cf mm@(GLM.GeneralizedLinearMixedModel _) pirlsType zStar vEta mzvu
  = P.wrapPrefix "compute_dBetaU" $ do
    P.logLE P.Diagnostic "Starting..."
    let mX   = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
        vW   = GLM.weights mm
        lf   = GLM.linkFunction $ linkFunctionType mm
        vMu  = VS.map (GLM.invLink lf) vEta
        vM   = VS.map (GLM.derivInv lf) vEta
        vVSW = GLM.varianceScaledWeights (observationsDistribution mm) vW vMu
    P.logLE P.Diagnostic $ "min(vMu)=" <> (T.pack $ show $ VS.minimum vMu) <> "; max(vMu)=" <> (T.pack $ show $ VS.maximum vMu)
    P.logLE P.Diagnostic $ "min(vM)=" <> (T.pack $ show $ VS.minimum vM) <> "; max(vM)=" <> (T.pack $ show $ VS.maximum vM)    
    P.logLE P.Diagnostic $ "min(vVSW)=" <> (T.pack $ show $ VS.minimum vVSW) <> "; max(vVSW)=" <> (T.pack $ show $ VS.maximum vVSW)
    P.logLE P.Diagnostic "Computing U and V"
    (msmU, mV) <- spUV mm pirlsType zStar vVSW vM
--    P.logLE P.Diagnostic $ "mVtV=" <> (T.pack $ show $ LA.tr mV LA.<> mV)
    let neqs = NormalEquationsGLMM vVSW vM lf vEta vMu mzvu
    (chol, vdBeta, mzvdu) <- cholmodCholeskySolutions' cf
                             (fmap SLA.transpose msmU)
                             mV
                             neqs
                             (GLM.mixedModelSpec mm)
    P.logLE P.Diagnostic $ "min(dBeta)=" <> (T.pack $ show $ VS.minimum vdBeta) <> "; max(dBeta)=" <> (T.pack $ show $ VS.maximum vdBeta)
    P.logLE P.Diagnostic $ "min(du)="
      <> (T.pack $ show $ fmap (VB.minimum . SLA.toVector) mzvdu)
      <> "; max(dU)="
      <> (T.pack $ show $ fmap (VB.maximum . SLA.toVector) mzvdu)    
    P.logLE P.Diagnostic "Finished."
    return (vdBeta, mzvdu, chol, msmU, mV)

compute_dBetaU _ (GLM.LinearMixedModel _) _ _ _ _ = P.throw
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
      f x = case lft of
        GLM.IdentityLink    -> x
        GLM.LogisticLink    -> GLM.bounded 0.00001 0.99999 x --minimum (maximum (0.0001, x), 0.99999)
        GLM.ExponentialLink -> max 0.0001 x
  in  VS.map (lF . f) vY

-- NB: eta = X <*> Beta + Z* <*> u
-- X <*> Beta = eta - Z* <*> u
-- X'X <*> Beta = X' <*> (eta - Z* <*> u)
-- Beta = inverse(X'X) <*> X' <*> (eta - Z* <*> u)
-- and X'X is symmetric positive definite so we solve the above via Cholesky

betaFrom
  :: GLM.MixedModel b g
  -> GLM.ZStarMatrix
  -> GLM.MaybeZeroUVec
  -> GLM.EtaVec
  -> GLM.BetaVec
betaFrom mm zStar mzvu vEta =
  let smZS = GLM.smZS zStar
      (_, q) = SLA.dim smZS
      vZSu = GLM.useMaybeZeroVec (VS.replicate q 0) (\svU -> SD.toDenseVector $ smZS SLA.#> svU) mzvu
      mX = GLM.rmsFixedPredictors $ GLM.regressionModelSpec mm
      vRhs = LA.tr mX LA.#> (vEta - vZSu)
      xTx = LA.tr mX LA.<> mX      
      cXtX = LA.chol $ LA.trustSym $ xTx --LA.tr mX LA.<> mX
  in  head $ LA.toColumns $ LA.cholSolve cXtX (LA.asColumn vRhs)


betaFrom'
  :: GLM.EffectsIO r
  => GLM.MixedModel b g
  -> GLM.ZStarMatrix -- this is unused and thus a sign that a redesign is in order.
  -> GLM.EtaVec
  -> GLM.MaybeZeroUVec
  -> P.Sem r GLM.BetaVec
betaFrom' mm zStar vEta mzvu = do
  let  vM   = VS.map (GLM.derivInv lf) vEta
       vMu = VS.map (GLM.invLink lf) vEta
       vY = GLM.observations mm
       vW   = GLM.weights mm
       vVSW = GLM.varianceScaledWeights (observationsDistribution mm) vW vMu
       lf =  GLM.linkFunction $ linkFunctionType mm
  (_, mV) <- spUV mm BetaAtZeroU zStar vVSW vM
  (_, vRhsX) <- normalEquationsRHS (NormalEquationsGLMM vW vM lf vEta vMu mzvu) Nothing mV vY
  let vTv = (LA.tr mV) LA.<> mV
      mRxx =  LA.chol $ LA.trustSym vTv
      betaSols = LA.cholSolve mRxx (LA.asColumn vRhsX)
      vBeta = head $ LA.toColumns betaSols
  return vBeta

{-
glmBeta
  :: GLM.FixedPredictors
  -> GLM.UVec
  -> GLM.EtaVec
  -> GLM.BetaVec
glmBeta mX zStar svU vEta =
  let vZSu = SD.toDenseVector $ (GLM.smZS zStar) SLA.#> svU
      vRhs = LA.tr mX LA.#> (vEta - vZSu)
      xTx = LA.tr mX LA.<> mX      
      cXtX = LA.chol $ LA.trustSym $ xTx --LA.tr mX LA.<> mX
  in  head $ LA.toColumns $ LA.cholSolve cXtX (LA.asColumn vRhs)  
-}

cholmodCholeskySolutions'
  :: GLM.EffectsIO r
  => CholmodFactor
  -> Maybe GLM.UMatrix -- Ut (ZSt in linear case). Nothing indicates we are solving at fixed u = 0.
  -> GLM.VMatrix -- V (X in linear case)
  -> NormalEquations
  -> GLM.MixedModelSpec b g
  -> P.Sem r (CholeskySolutions, GLM.BetaVec, GLM.MaybeZeroUVec)
cholmodCholeskySolutions' cholmodFactor smUtM mV nes mixedModelSpec =
  P.wrapPrefix "cholmodCholeskySolutions'" $ do
    let (cholmodC, cholmodF, smP) = cholmodFactor
        vY = GLM.rmsObservations $ GLM.mmsRegressionModelSpec mixedModelSpec
        n                         = LA.size vY
        (_, p)                    = LA.size mV
--        (q, _)                    = SLA.dim smUt
    -- Cholesky factorize to get L_theta *and* update factor for solving with it
    P.logLE P.Diagnostic "factorizing..."
    -- TODO: Cholmod has built in support for factorizing XtX + a * I.  Use it.
    let cfs = CH.FactorizeAtAPlusBetaI 1
    
    (mzvRhsZ, vRhsX) <- normalEquationsRHS nes smUtM (LA.tr mV) vY
{-    
    P.logLE P.Diagnostic
      $  "Normal Equation RHS:\nsvRhsZ="
      <> (T.pack $ show svRhsZ)
      <> "\nvRhsX="
      <> (T.pack $ show vRhsX)
-}
--    P.logLE P.Diagnostic "Calling solveSparse for svC."
    let vTv = (LA.tr mV) LA.<> mV
    case smUtM of
      Nothing -> do
        P.logLE P.Diagnostic "CholmodCholeskySolutions' for beta only (at fixed u=0)"
        let mRxx = LA.chol $ LA.trustSym vTv
            betaSols = LA.cholSolve mRxx (LA.asColumn vRhsX)
            vBeta = head $ LA.toColumns betaSols
            q = zCols $ GLM.mmsFitSpecByGroup mixedModelSpec
        return (CholeskySolutions (SLA.eye q) (SLA.zeroSM q p) mRxx, vBeta, GLM.zeroVec)
      Just smUt -> do
        let (q, _)                    = SLA.dim smUt
            svRhsZ = GLM.useMaybeZeroVec (SLA.zeroSV q) id mzvRhsZ
        liftIO $ CH.spMatrixFactorizeP cholmodC cholmodF cfs CH.UnSymmetric smUt        
        smPRhsZ         <- liftIO
          $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P (SD.svColumnToSM svRhsZ)
  --    P.logLE P.Diagnostic "After smPRhsZ..."
        svC <-
          SLA.toSV
          <$> (liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_LD smPRhsZ)
        let smUtV = smUt SLA.#~# (SD.toSparseMatrix mV)
            --    P.logLE P.Diagnostic "Calling solveSparse for smPUtV."
        smPUtV <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smUtV -- PU'V
--    P.logLE P.Diagnostic "Calling solveSparse for smRzx."
        smRzx  <- liftIO $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_LD smPUtV -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
    -- compute Rxx
    -- NB: c is defined as Lu + Rzx(Beta) which allows solving the normal equations in pieces
        let vTvMinusRzxTRzx =
              ((LA.tr mV) LA.<> mV)  - (SD.toDenseMatrix $ smRzx SLA.#^# smRzx)
--    liftIO $ putStrLn $ "mV=" ++ show mV
--    P.logLE P.Diagnostic $ "Rzx=" <> (T.pack $ show smRzx)
--    P.logLE P.Diagnostic $ "vtVMinusRzxTRzx=" <> (T.pack $ show vTvMinusRzxTRzx)
        let mRxx            = LA.chol $ LA.trustSym $ vTvMinusRzxTRzx -- mRxx' mRxx = v'v - Rzx'Rzx
--    P.logLE P.Diagnostic $ "mRxx=" <> (T.pack $ show mRxx)
    
    -- now we have the entire Cholesky decomposition.  Solve the normal equations
        let vRzxtC          = SD.toDenseVector $ (SLA.transposeSM smRzx) SLA.#> svC
            vBetaRhs         = vRhsX - vRzxtC
--    P.logLE P.Diagnostic $ "vRhsX=" <> (T.pack $ show vRhsX)
--    P.logLE P.Diagnostic $ "vRzxtC=" <> (T.pack $ show vRzxtC)
--    P.logLE P.Diagnostic $ "vBetaRHS=" <> (T.pack $ show vBetaRhs)    
        let betaSols        = LA.cholSolve mRxx (LA.asColumn $ vBetaRhs)
            vBeta           = head $ LA.toColumns $ betaSols
            svBeta          = SD.toSparseVector vBeta
            svRzxBeta       = smRzx SLA.#> svBeta
            smCMinusRzxBeta = SD.svColumnToSM (svC SLA.^-^ svRzxBeta)
--    P.logLE P.Diagnostic "Calling solveSparse for smPu."
        smPu <- liftIO
          $ CH.solveSparse cholmodC cholmodF CH.CHOLMOD_DLt smCMinusRzxBeta
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
        return $ (CholeskySolutions smLth smRzx mRxx, vBeta, GLM.MaybeZeroVec (Just svu))



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
