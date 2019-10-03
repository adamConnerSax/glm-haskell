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
import qualified Data.IndexedSet               as IS

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P

--import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Colonnade                     as C

import qualified Data.Array                    as A
import qualified Data.List                     as L
import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
import qualified Data.Sparse.Common            as SLA

{-
import qualified Numeric.LinearAlgebra.Class   as SLA
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
--import qualified Numeric.NLOPT                 as NL



import qualified Data.Map                      as M
import           Data.Maybe                     ( isJust )
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

type FixedPredictors = LA.Matrix Double
type Observations = LA.Vector Double

data RegressionModel b = RegressionModel (GLM.FixedEffects b) FixedPredictors Observations

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

data MixedModel b g = MixedModel (RegressionModel b) (FitSpecByGroup g)

type RandomEffectModelMatrix = SLA.SpMatrix Double
type CovarianceVector = LA.Vector Double
type Lambda = SLA.SpMatrix Double

type WMatrix = LA.Vector Double -- constant weights
type UMatrix = SLA.SpMatrix Double
type VMatrix = LA.Matrix Double

type UVec = LA.Vector Double
type MuVec = LA.Vector Double

data RandomEffectCalculated = RandomEffectCalculated RandomEffectModelMatrix (CovarianceVector -> Lambda)

type PermutationMatrix = SLA.SpMatrix Double -- should this be (SLA.SpMatrix Int) ??

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

setCovarianceVector :: FitSpecByGroup g -> Double -> Double -> LA.Vector Double
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
makeLambda :: FitSpecByGroup g -> (CovarianceVector -> Lambda)
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
                     , PermutationMatrix -- ^ permutation matrix from above
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
                                           , rX :: Rxx
                                           , svBeta :: SLA.SpVector Double
                                           , svu :: SLA.SpVector Double
                                           } deriving (Show)

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


cholmodCholeskySolutions'
  :: CholmodFactor
  -> UMatrix -- U (ZS in linear case)
  -> VMatrix -- V (X in linear case)
  -> NormalEquations
  -> MixedModel b g
  -> IO CholeskySolutions
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
  return $ CholeskySolutions smLth smRzx mRxx svBeta svu

-- fitted values of input data
-- like "fitted" in lme4
fitted
  :: (Ord g, Show g, Show b, Ord b, Enum b, Bounded b)
  => (r -> b -> Double)
  -> (r -> g -> T.Text)
  -> GLM.FixedEffectStatistics b
  -> GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> r
  -> Either T.Text Double -- the map in effectParametersByGroup and the map in RowClassifier might not match
fitted getPred getLabel (GLM.FixedEffectStatistics fe vFE _) ebg rowClassifier row
  = do
    let applyEffect b e r = case b of
          GLM.Intercept   -> e
          GLM.Predictor b -> e * getPred r b
        applyEffects ies v r =
          FL.fold FL.sum
            $ fmap
                (\(index, effect) -> applyEffect effect (v `LA.atIndex` index) r
                )
            $ IS.toIndexedList ies
        groups = M.keys ebg
        applyGroupEffects r group = do
          labelIndex <- GLM.labelIndex rowClassifier group (getLabel r group)
          (GLM.EffectParameters ies mEP) <- GLM.effectParameters group ebg
          return $ applyEffects ies (LA.toRows mEP L.!! labelIndex) r
        fixedEffects = applyEffects (GLM.indexedFixedEffectSet fe) vFE row
    groupEffects <- traverse (applyGroupEffects row) $ M.keys ebg
    return $ fixedEffects + FL.fold FL.sum groupEffects


-- Like "ranef" in lme4
randomEffectsByLabel
  :: (Ord g, Show g, Show b, Enum b, Bounded b)
  => GLM.EffectParametersByGroup g b
  -> GLM.RowClassifier g
  -> Either
       T.Text
       (M.Map g (GLM.IndexedEffectSet b, [(T.Text, LA.Vector Double)]))
randomEffectsByLabel ebg rowClassifier = do
  let byLabel group (GLM.EffectParameters ies mEP) = do
        indexedLabels <- M.toList <$> GLM.labelMap rowClassifier group
        return
          $ (ies, fmap (\(l, r) -> (l, LA.toRows mEP L.!! r)) indexedLabels)
  M.traverseWithKey byLabel ebg


colRandomEffectsByLabel
  :: Show b
  => GLM.IndexedEffectSet b
  -> C.Colonnade C.Headed (T.Text, LA.Vector Double) T.Text
colRandomEffectsByLabel ies =
  let
    indexedEffects = IS.toIndexedList ies
    colLabel       = C.headed "Category" fst
    colEffect (index, effect) = C.headed
      (T.pack $ show effect)
      (T.pack . show . flip LA.atIndex index . snd)
  in
    mconcat $ (colLabel : fmap colEffect indexedEffects)


printRandomEffectsByLabel
  :: (Show g, Show b)
  => M.Map g (GLM.IndexedEffectSet b, [(T.Text, LA.Vector Double)])
  -> T.Text
printRandomEffectsByLabel rebg =
  let printOne (g, (ies, x)) =
        (T.pack $ show g)
          <> ":\n"
          <> (T.pack $ C.ascii (fmap T.unpack $ colRandomEffectsByLabel ies) x)
          <> "\n"
  in  mconcat $ fmap printOne $ M.toList rebg


---- Unused
cholmodCholeskySolutions_Old
  :: CholmodFactor
  -> MixedModel b g
  -> RandomEffectCalculated
  -> CovarianceVector
  -> IO CholeskySolutions
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
  return $ CholeskySolutions smLth smRzx mRx svBeta svu
