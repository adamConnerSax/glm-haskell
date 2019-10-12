{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
module Numeric.GLM.MixedModel where

import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH
import qualified Numeric.SparseDenseConversions
                                               as SD
import qualified Numeric.GLM.Types             as GLM
import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Data.IndexedSet               as IS

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P

--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Data.Array                    as A
import qualified Data.List                     as L
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA

import           Numeric.LinearAlgebra.Sparse   ( (#>) )

import qualified Numeric.LinearAlgebra         as LA
import qualified Numeric.NLOPT                 as NL
import           System.IO.Unsafe               ( unsafePerformIO )



import qualified Data.Map                      as M
import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

type FixedPredictors = LA.Matrix Double
type Observations = LA.Vector Double



-- for group k, what group effects are we modeling?
-- first Bool is for intercept, then (optional) vector,
-- of same length as X has columns, to indicate which
-- predictors in X get random slopes
-- NB:  If X has a constant column, there is redundancy
-- here between the Bool in the tuple and the first bool
-- in the vector.
data GroupFitSpec = GroupFitSpec { nCategories :: Int
                                 , groupIntercept :: Bool
                                 , groupSlopes :: Maybe (VB.Vector Bool)
                                 } deriving (Show, Eq)

makeGroupFitSpec
  :: (Ord b, Enum b, Bounded b)
  => Int -- ^ number of items in group
  -> GLM.FixedEffects b
  -> GLM.IndexedEffectSet b
  -> Either T.Text GroupFitSpec
makeGroupFitSpec n fixedEffects indexedGroupEffects = do
  let indexedFixedEffects = GLM.indexedFixedEffectSet fixedEffects
  when (not $ IS.subset indexedGroupEffects indexedFixedEffects) $ Left
    "group contains effects not in fixed effects in \"makeGroupFitSpec\""
  let modeled        = isJust . IS.index indexedGroupEffects
      groupIntercept = modeled GLM.Intercept
      modeledNI x = if x == GLM.Intercept then False else modeled x
      groupSlopesL = fmap modeledNI $ IS.members indexedFixedEffects
      groupSlopes  = if True `L.elem` groupSlopesL
        then Just (VB.fromList groupSlopesL)
        else Nothing
  return $ GroupFitSpec n groupIntercept groupSlopes

type FitSpecByGroup g = M.Map g GroupFitSpec --VB.Vector GroupFitSpec

fitSpecByGroup
  :: (Ord g, Show g, Ord b, Enum b, Bounded b, Show b)
  => GLM.FixedEffects b
  -> GLM.EffectsByGroup g b
  -> GLM.RowClassifier g
  -> Either T.Text (FitSpecByGroup g)
fitSpecByGroup fixedEffects ebg rowClassifier = do
  let lookupGroupSize (grp, ie) =
        M.lookup grp (GLM.groupSizes rowClassifier) >>= return . (grp, , ie)
      makeSpec (grp, n, ige) =
        makeGroupFitSpec n fixedEffects ige >>= return . (grp, )
      lookupError =
        "Failed to find a group in the group effect set ("
          <> (T.pack $ show ebg)
          <> ") in  the rowClassifier ("
          <> (T.pack $ show rowClassifier)
  groupInfoList <-
    maybe (Left lookupError) Right $ traverse lookupGroupSize $ M.toList ebg
  fmap M.fromList $ traverse makeSpec groupInfoList

data RegressionModel b = RegressionModel (GLM.FixedEffects b) FixedPredictors Observations deriving (Show, Eq)
data MixedModel b g = MixedModel (RegressionModel b) (FitSpecByGroup g) deriving (Show, Eq)

data GeneralizedLinearMixedModel b g =
  LMM (MixedModel b g)
  | GLMM (MixedModel b g) WMatrix GLM.ObservationDistribution deriving (Show, Eq)

mixedModel :: GeneralizedLinearMixedModel b g -> MixedModel b g
mixedModel (LMM x     ) = x
mixedModel (GLMM x _ _) = x

type RandomEffectModelMatrix = SLA.SpMatrix Double
type CovarianceVec = LA.Vector Double
type LambdaMatrix = SLA.SpMatrix Double

type WMatrix = LA.Vector Double -- constant weights
type UMatrix = SLA.SpMatrix Double
type VMatrix = LA.Matrix Double
type LMatrix = SLA.SpMatrix Double -- Cholesky Z block
type RzxMatrix = SLA.SpMatrix Double
type RxxMatrix = LA.Matrix Double

type BetaVec = LA.Vector Double
type UVec = LA.Vector Double
type MuVec = LA.Vector Double

data BetaU vt = BetaU { betaVec :: vt Double, uVec :: vt Double }
type DenseBetaU = BetaU LA.Vector
type SparseBetaU = BetaU SLA.SpVector

toDenseBetaU :: SparseBetaU -> DenseBetaU
toDenseBetaU (BetaU svBeta svU) =
  BetaU (SD.toDenseVector svBeta) (SD.toDenseVector svU)

toSparseBetaU :: DenseBetaU -> SparseBetaU
toSparseBetaU (BetaU vBeta vU) =
  BetaU (SD.toSparseVector vBeta) (SD.toSparseVector vU)

unaryOpBetaU :: (vt Double -> vt Double) -> BetaU vt -> BetaU vt
unaryOpBetaU uop (BetaU b u) = BetaU (uop b) (uop u)

binaryOpBetaU
  :: (vt Double -> vt Double -> vt Double) -> BetaU vt -> BetaU vt -> BetaU vt
binaryOpBetaU bOp (BetaU bA uA) (BetaU bB uB) = BetaU (bOp bA bB) (bOp uA uB)

data RandomEffectCalculated = RandomEffectCalculated RandomEffectModelMatrix (CovarianceVec -> LambdaMatrix)

type PMatrix = SLA.SpMatrix Double -- should this be (SLA.SpMatrix Int) ??

data DevianceType = ML | REML deriving (Show, Eq)

effectsForGroup :: GroupFitSpec -> Int
effectsForGroup (GroupFitSpec _ b vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForGroup #-}

colsForGroup :: GroupFitSpec -> Int
colsForGroup l = nCategories l * effectsForGroup l
{-# INLINABLE colsForGroup #-}

data ParameterEstimates =
  ParameterEstimates
  { pEstimate :: LA.Vector Double
  , pCovariance :: LA.Matrix Double
  }

type FixedParameterEstimates = ParameterEstimates
type GroupParameterEstimates = VB.Vector FixedParameterEstimates

type SemC r = (MonadIO (P.Sem r), P.Member (P.Error T.Text) r)
runPIRLS_M :: P.Sem '[P.Error T.Text, P.Lift IO] a -> IO (Either T.Text a)
runPIRLS_M = P.runM . P.runError {-. P.runErrorAsAnother (T.pack . show)-}

data MinimizeDevianceVerbosity = MDVNone | MDVSimple

minimizeDeviance
  :: SemC r
  => MinimizeDevianceVerbosity
  -> DevianceType
  -> GeneralizedLinearMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec -- ^ initial guess for theta
  -> P.Sem
       r
       ( CovarianceVec
       , Double
       , Double
       , DenseBetaU
       , LA.Vector Double
       , CholeskySolutions
       ) -- ^ (theta, profiled_deviance, sigma2, beta, u, b, cholesky blocks)
minimizeDeviance verbosity dt gmm reCalc@(RandomEffectCalculated smZ mkLambda) th0
  = do
    cholmodAnalysis <- cholmodAnalyzeProblem reCalc
    let
      pdv = case verbosity of
        MDVNone   -> PDVNone
        MDVSimple -> PDVSimple
      pd x = profiledDeviance pdv cholmodAnalysis dt gmm reCalc x
      obj x = unsafePerformIO $ fmap (\(d, _, _, _, _) -> d) $ pd x
      stop                = NL.ObjectiveAbsoluteTolerance 1e-6 NL.:| []
      MixedModel _ levels = mixedModel gmm
      thetaLB             = thetaLowerBounds levels
      algorithm           = NL.BOBYQA obj [thetaLB] Nothing
      problem = NL.LocalProblem (fromIntegral $ LA.size th0) stop algorithm
      expThetaLength      = FL.fold
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
        (pd, sigma2, svBetaU, svb, cs) <- pd thS
        return (thS, pd, sigma2, toDenseBetaU svBetaU, SD.toDenseVector svb, cs)


-- ugh.  But I dont know a way in NLOPT to have bounds on some not others.
thetaLowerBounds :: FitSpecByGroup g -> NL.Bounds
thetaLowerBounds groupFSM =
  let negInfinity :: Double = negate $ 1 / 0
  in  NL.LowerBounds $ setCovarianceVector groupFSM 0 negInfinity --FL.fold fld levels

setCovarianceVector :: FitSpecByGroup g -> Double -> Double -> CovarianceVec
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
  :: forall g
   . (Enum g, Bounded g, A.Ix g)
  => FixedPredictors
  -> FitSpecByGroup g
  -> GLM.RowClassifier g
  -> Maybe RandomEffectModelMatrix
makeZ mX groupFSM rc = do
  let
    (nO, nP) = LA.size mX
    k        = FL.fold FL.length groupFSM -- number of levels
--    groupSize g = fmap nCategories $ M.lookup g groupFSM 
    q        = FL.fold FL.sum $ fmap colsForGroup groupFSM -- total number of columns in Z
    predictor rowIndex fixedEffectIndex =
      mX `LA.atIndex` (rowIndex, fixedEffectIndex)
    -- construct Z for level as a fold over rows
    entries :: Int -> g -> GroupFitSpec -> Int -> Maybe [(Int, Int, Double)]
    entries colOffset group groupFS rowIndex = do
      categoryN <- GLM.categoryNumberFromRowIndex rc rowIndex group
      let
        entryOffset = colOffset + (effectsForGroup groupFS * categoryN)
        (intercept, slopeOffset) = if groupIntercept groupFS
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
        slopes = case groupSlopes groupFS of
          Just slopeV -> FL.fold slopeF $ VB.zip slopeV (VB.generate nP id)
          Nothing     -> []
      return $ intercept ++ slopes
    ziFoldM
      :: g -> Int -> GroupFitSpec -> FL.FoldM Maybe Int [(Int, Int, Double)]
    ziFoldM group groupOffset groupFS =
      let q = colsForGroup groupFS
      in  FL.FoldM
            (\mes n -> fmap (mes ++) $ entries groupOffset group groupFS n)
            (Just [])
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
checkProblem :: SemC r => MixedModel b g -> RandomEffectCalculated -> P.Sem r ()
checkProblem (MixedModel (RegressionModel _ mX vY) _) (RandomEffectCalculated smZ _)
  = do
    let (n, _)     = LA.size mX
        yRows      = LA.size vY
        (zRows, _) = SLA.dim smZ
    when (zRows /= n)
      $  P.throw @T.Text
      $  (T.pack $ show n)
      <> " cols in X but "
      <> (T.pack $ show zRows)
      <> " cols in Z"
    when (yRows /= n)
      $  P.throw @T.Text
      $  (T.pack $ show n)
      <> " cols in X but "
      <> (T.pack $ show yRows)
      <> " entries in Y"
    return ()


-- For each level, i, with p_i random effects,
-- Lambda is built from lower-triangular blocks
-- which are p_i x p_i.  Theta contains the values
-- for each block, in col major order
makeLambda :: FitSpecByGroup g -> (CovarianceVec -> LambdaMatrix)
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
      perGroup (l, vx) = replicate (nCategories l) $ templateBlock l vx
      allDiags vTh = concat $ fmap perGroup $ blockData vTh
  in  (\vTh -> SLA.fromBlocksDiag (allDiags vTh))

xTxPlusI :: SLA.SpMatrix Double -> SLA.SpMatrix Double
xTxPlusI smX = (smX SLA.#~^# smX) SLA.^+^ (SLA.eye $ SLA.ncols smX)

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

type CholmodFactor = (CH.ForeignPtr CH.Common -- ^ pre-allocated CHOMOD common space
                     , CH.ForeignPtr CH.Factor -- ^ precomputed pattern work on Z
                     , PMatrix -- ^ permutation matrix from above
                     )

cholmodAnalyzeProblem
  :: SemC r => RandomEffectCalculated -> P.Sem r CholmodFactor
cholmodAnalyzeProblem (RandomEffectCalculated smZ _) = do
  cholmodC <- liftIO CH.allocCommon
  liftIO $ CH.startC cholmodC
  liftIO $ CH.setFinalLL 1 cholmodC -- I don't quite get this.  We should be able to solve LDx = b either way??
  (cholmodF, smP) <-
    liftIO $ CH.spMatrixAnalyzeWP cholmodC CH.SquareSymmetricLower $ xTxPlusI
      smZ
--  liftIO $ CH.printFactor cholmodF "After analyze, before factorize" cholmodC
--  liftIO $ putStrLn "smP=" >> (LA.disp 1 $ SD.toDenseMatrix smP)
  return (cholmodC, cholmodF, smP)

data ProfiledDevianceVerbosity = PDVNone | PDVSimple | PDVAll deriving (Show, Eq)

profiledDeviance
  :: ProfiledDevianceVerbosity
  -> CholmodFactor
  -> DevianceType
  -> GeneralizedLinearMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec -- ^ theta
  -> IO
       ( Double
       , Double
       , SparseBetaU
       , SLA.SpVector Double
       , CholeskySolutions
       ) -- ^ (pd, sigma^2, beta, u, b) 
profiledDeviance verbosity cf dt gmm reCalc@(RandomEffectCalculated smZ mkLambda) vTh
  = do
    let MixedModel (RegressionModel _ mX vY) _ = mixedModel gmm
        n      = LA.size vY
        (_, p) = LA.size mX
        (_, q) = SLA.dim smZ
        lambda = mkLambda vTh
        smZS   = smZ SLA.## lambda --(smT SLA.## smS)
--        smZSt  = SLA.transpose smZS
    (cs@(CholeskySolutions smLth smRzx mRxx), svBetaU) <- getCholeskySolutions
      cf
      gmm
      GLM.UseCanonical -- FIX
      reCalc
      vTh
    when (verbosity == PDVAll) $ do
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
    let
      svb    = lambda SLA.#> (uVec svBetaU) -- (smT SLA.## smS) SLA.#> svu
      vBetaU = toDenseBetaU svBetaU
--        vBeta = SD.toDenseVector svBeta
--        vU    = SD.toDenseVector svu
      (pd, rTheta2) =
        profiledDeviance' verbosity dt gmm reCalc vTh smLth mRxx vBetaU
      dof = case dt of
        ML   -> realToFrac n
        REML -> realToFrac (n - p)
      sigma2 = rTheta2 / dof
    when (verbosity == PDVAll) $ do
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
    when (verbosity == PDVAll || verbosity == PDVSimple) $ do
      putStrLn $ "pd(th=" ++ (show vTh) ++ ") = " ++ show pd
    return (pd, sigma2, svBetaU, svb, cs)

profiledDeviance'
  :: ProfiledDevianceVerbosity
  -> DevianceType
  -> GeneralizedLinearMixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> LMatrix
  -> RxxMatrix
  -> DenseBetaU
  -> (Double, Double) -- profiled deviance, deviance + u'u 
profiledDeviance' pdv dt gmm re vTh smLth mRxx vBetaU =
  let MixedModel (RegressionModel _ mX vY) _ = mixedModel gmm
      (          RandomEffectCalculated smZ                    mkLambda) = re
      observationDistribution = case gmm of
        LMM _       -> GLM.Normal
        GLMM _ _ od -> od
      vBeta         = betaVec vBetaU
      vU            = uVec vBetaU
      lambda        = mkLambda vTh
      smZS          = smZ SLA.#~# lambda
      linPred       = linearPredictor mX smZS vBetaU
      logLth        = logDetTriangularSM smLth
      dev2 = GLM.deviance observationDistribution GLM.UseCanonical vY linPred
      rTheta2       = dev2 + (vU LA.<.> vU)
      n             = LA.size vY
      (_  , p     ) = LA.size mX
      (dof, logDet) = case dt of
        ML   -> (realToFrac n, logLth)
        REML -> (realToFrac (n - p), logLth + (logDetTriangularM mRxx))
      pd = (2 * logDet) + (dof * (1 + log (2 * pi * rTheta2 / dof)))
  in  (pd, rTheta2)


type LinearPredictor = LA.Vector Double

linearPredictor
  :: FixedPredictors -> SLA.SpMatrix Double -> DenseBetaU -> LinearPredictor
linearPredictor mX smZS vBetaU =
  mX
    LA.#> (betaVec vBetaU)
    +     (SD.toDenseVector $ smZS SLA.#> SD.toSparseVector (uVec vBetaU))

upperTriangular :: Int -> Int -> a -> Bool
upperTriangular r c _ = (r <= c)

lowerTriangular :: Int -> Int -> a -> Bool
lowerTriangular r c _ = (r >= c)

type LTheta = SLA.SpMatrix Double
type Rzx = SLA.SpMatrix Double
type Rxx = LA.Matrix Double
type Beta = LA.Vector Double

data CholeskySolutions = CholeskySolutions { lTheta :: LTheta
                                           , rZX :: Rzx
                                           , rXX :: Rxx
                                           } deriving (Show)
{-
                                           , svBeta :: SLA.SpVector Double
                                           , svu :: SLA.SpVector Double
                                           } deriving (Show)
-}

data NormalEquations = NormalEquationsLMM | NormalEquationsGLMM WMatrix MuVec UVec

normalEquationsRHS
  :: NormalEquations
  -> UMatrix
  -> VMatrix
  -> LA.Vector Double
  -> (SLA.SpVector Double, LA.Vector Double)
normalEquationsRHS NormalEquationsLMM smU mV vY =
  (SLA.transpose smU SLA.#> SD.toSparseVector vY, LA.tr mV LA.#> vY)

normalEquationsRHS (NormalEquationsGLMM vW vMu vU) smU mV vY =
  let v1   = VS.zipWith (*) (VS.map sqrt vW) (vY - vMu)
      rhsZ = (SLA.transpose smU SLA.#> SD.toSparseVector v1)
        SLA.^-^ SD.toSparseVector vU
      rhsX = LA.tr mV LA.#> v1
  in  (rhsZ, rhsX)


cholmodCholeskySolutionsLMM
  :: CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> IO (CholeskySolutions, SparseBetaU)
cholmodCholeskySolutionsLMM cholmodFactor mixedModel randomEffCalcs vTh = do
  let (MixedModel             (RegressionModel _ mX vY) levels  ) = mixedModel
      (RandomEffectCalculated smZ mkLambda) = randomEffCalcs
      lambda = mkLambda vTh
      smZS   = smZ SLA.## lambda
  cholmodCholeskySolutions' cholmodFactor smZS mX NormalEquationsLMM mixedModel

linkFunction
  :: GeneralizedLinearMixedModel b g -> GLM.UseLink -> GLM.LinkFunction
linkFunction glmm ul = GLM.linkFunction $ case ul of
  GLM.UseOther x   -> x
  GLM.UseCanonical -> GLM.canonicalLink $ case glmm of
    LMM _       -> GLM.Normal
    GLMM _ _ od -> od

spUV
  :: GeneralizedLinearMixedModel b g
  -> GLM.UseLink
  -> RandomEffectCalculated
  -> CovarianceVec
  -> SparseBetaU
  -> (UMatrix, VMatrix)

spUV (LMM (MixedModel (RegressionModel _ mX _) _)) _ (RandomEffectCalculated smZ mkLambda) vTh _
  = (smZ SLA.#~# mkLambda vTh, mX)

spUV glmm@(GLMM (MixedModel (RegressionModel _ mX vY) _) vW _) ul (RandomEffectCalculated smZ mkLambda) vTh svBetaU
  = let (GLM.LinkFunction l inv dinv) = linkFunction glmm ul --GLM.linkFunction (GLM.canonicalLink od)
        spLambda                      = mkLambda vTh
        vEta =
          SD.toDenseVector (smZ #> uVec svBetaU)
            +     mX
            LA.#> (SD.toDenseVector $ betaVec svBetaU)
        vdMudEta  = VS.map dinv vEta
        mWdMudEta = LA.diag $ VS.zipWith (*) vW vdMudEta
        smU       = SD.toSparseMatrix mWdMudEta SLA.#~# smZ SLA.#~# spLambda
        mV        = mWdMudEta LA.<> mX
    in  (smU, mV)

getCholeskySolutions
  :: CholmodFactor
  -> GeneralizedLinearMixedModel b g
  -> GLM.UseLink
  -> RandomEffectCalculated
  -> CovarianceVec
  -> IO (CholeskySolutions, SparseBetaU)
getCholeskySolutions cf (LMM mm) _ reCalc vTh =
  cholmodCholeskySolutionsLMM cf mm reCalc vTh

{-
-- TODO: There is much low hanging fruit here in terms of not re-calcing things
-- 
getCholeskySolutions cf glmm@(GLMM mm vW od) ul reCalc vTH =
  let
    MixedModel (RegressionModel _ mX vY) levels   = mm
    (RandomEffectCalculated smZ mkLambda) = reCalc
    smZS = smZ SLA.#~# mkLamda vTh
    n = LA.size vY
    (_,q) = SLA.dim smZS  
    vBeta0 = LA.fromList $ iterate n 0 -- FIXME (fit fixed effects only and use that for beta0)
    vU0  = LA.fromList $ iterate n 0 -- FIXME (maybe 0 is right here??)
    lf = linkFunction glmm ul
    deltaCCS svBetaU = do
      let (smU, mV) = spUV glmm ul reCalc vTh svBetaU
          vEta = linearPredictor mX smZS (toDenseBetaU svBetaU)
          vMu = VS.map (invLink lf) vEta
          neqs = NormalEquationsGLMM vW vMu (SD.toDenseVector $ vecU svBetaU)          
      (chol, svBetaU') <- cholmodCholeskySolutions' cf smU mV neqs mm
      return (svBetaU', lTheta solns, rXX solns)
    incS svBetaU svdBetaU x =
      binaryOpBetaU (SLA.^+^) svBetaU (unaryOpBetaU (x .*) svdBetaU)
    newBetaU (svBeta, svU) (svdBeta, svdU) smLth mRxx = 
      let pdF = profiledDeviance' PDVNone ML glmm reCalc vTh smLth mRxx
          pd0 = pdF svBeta svU
          checkOne x = 
            let (svBeta', svU') = incS (svBeta, svU) (svdBeta, svdU) x
            in if pdF svBeta' svU' < pd0 then (Just (svBeta', svU')) else Nothing
          check triesLeft x = case checkOne x of
            Nothing -> if n > 1 then check (n-1) x/2 else Nothing
            Just y -> Just y
      in check 10 1.0

orthogonalConvergence
  :: GeneralizedMixedModel b g
  -> GLM.UseLink
  -> SLA.SpMatrix Double -- Z*Lambda(theta)
  -> PMatrix
  -> LMatrix
  -> VMatrix
  -> UVec
  -> BetaVec
  -> SLA.SpVector Double
orthogonalConvergence glmm ul smZS smP smL vU vBeta svdU =
  let MixedModel (RegressionModel _ mX vY) _ = mixedModel glmm
      lft = case ul of
        UseOther x -> x
        UseCanonical -> canonicalLink $ case glmm of
          LMM -> Normal
          GLMM _ _ od -> od
      vEta = linearPredictor mX smZS vBeta vU
      vMu = VS.map (invLink $ linkFunction lft) vEta
      x1 = SLA.transpose smL SLA.#> svdU
      x2 = smP SLA.#~# (SLA.


glmmCholeskyStep
  :: CholmodFactor
  -> NormalEquations
  -> MixedModel
  -> WMatrix
  -> LinkFunctionType
  -> SLA.SpVector Double
  -> SLA.SpVector Double
  -> IO CholeskySolutions
glmmCholeskyStep mm vW lft svBeta svU =
-}


cholmodCholeskySolutions'
  :: CholmodFactor
  -> UMatrix -- U (ZS in linear case)
  -> VMatrix -- V (X in linear case)
  -> NormalEquations
  -> MixedModel b g
  -> IO (CholeskySolutions, SparseBetaU)
cholmodCholeskySolutions' cholmodFactor smU mV nes mixedModel = do
  let (cholmodC, cholmodF, _) = cholmodFactor
      (MixedModel (RegressionModel _ _ vY) _) = mixedModel
      n                       = LA.size vY
      (_, p)                  = LA.size mV
      (_, q)                  = SLA.dim smU
  -- Cholesky factorize to get L_theta *and* update factor for solving with it
  CH.spMatrixFactorize cholmodC cholmodF CH.SquareSymmetricLower
    $ SLA.filterSM lowerTriangular
    $ xTxPlusI
    $ smU
  -- compute Rzx
  let (svRhsZ, vRhsX) = normalEquationsRHS nes smU mV vY
  smPRhsZ <- CH.solveSparse cholmodC
                            cholmodF
                            CH.CHOLMOD_P
                            (SD.svColumnToSM svRhsZ)
  svC <- SLA.toSV <$> CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPRhsZ
  let smUtV = (SLA.transpose smU) SLA.#~# (SD.toSparseMatrix mV)
  smPUtV <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_P smUtV -- PU'V
  smRzx  <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_L smPUtV -- NB: If decomp was LL' then D is I but if it was LDL', we need D...
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
  smPu  <- CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Lt smCMinusRzxBeta
  svu   <- SLA.toSV <$> (CH.solveSparse cholmodC cholmodF CH.CHOLMOD_Pt $ smPu)
  -- NB: This has to happen after the solves because it unfactors the factor...
  -- NB: It does *not* undo the analysis
  smLth <- CH.choleskyFactorSM cholmodF cholmodC
  return $ (CholeskySolutions smLth smRzx mRxx, BetaU svBeta svu)



---- Unused
cholmodCholeskySolutions_Old
  :: CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVec
  -> IO (CholeskySolutions, SparseBetaU)
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
