{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.PIRLS where

import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import           Polysemy.ConstraintAbsorber.MonadCatch
                                                ( absorbMonadThrow
                                                , SomeException
                                                , Exception(..)
                                                )
import qualified GLM.Internal.Log              as Log
import qualified Control.Foldl                 as FL
import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )


import qualified Data.Sparse.SpMatrix          as SLA
import qualified Data.Sparse.SpVector          as SLA
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



import           Data.Either                    ( partitionEithers )
import qualified Data.List                     as L
import qualified Data.Sequence                 as Seq
import qualified Data.Text                     as T
import qualified Data.Vector                   as VB
import qualified Data.Vector.Storable          as VS

type LevelSpec = (Int, Bool, Maybe (VB.Vector Bool)) -- level sizes and effects

effectsForLevel :: LevelSpec -> Int
effectsForLevel (_, b, vbM) =
  (if b then 1 else 0) + (maybe 0 (length . VB.filter id) vbM)
{-# INLINABLE effectsForLevel #-}

colsForLevel :: LevelSpec -> Int
colsForLevel l@(qL, _, _) = qL * effectsForLevel l
{-# INLINABLE colsForLevel #-}

type Levels = VB.Vector LevelSpec

-- classify row into its levels
-- the vector has a level number for each level
type RowClassifier = Int -> VB.Vector Int

-- for group k, what group effects are we modeling?
-- first Bool is for intercept, then (optional) vector,
-- of same length as X has columns, to indicate which
-- predictors in X get random slopes
-- NB:  If X has a constant column, there is redundancy
-- here between the Bool in the tuple and the first bool
-- in the vector.
--type LevelEffectsSpec = [(Bool, Maybe (LA.Vector Bool))]

type SemC r = (MonadIO (P.Sem r), {- P.Member (P.Error SomeException) r,-} P.Member (P.Error T.Text) r)
runPIRLS_M
  :: P.Sem '[{-P.Error SomeException,-}
             P.Error T.Text, P.Lift IO] a
  -> IO (Either T.Text a)
runPIRLS_M = P.runM . P.runError {-. P.runErrorAsAnother (T.pack . show)-}

makeZ
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double)
  => LA.Matrix a
  -> Levels
  -> RowClassifier
  -> SLA.SpMatrix a
makeZ mX levels rc =
  let
    (nO, nP) = LA.size mX
    k        = VB.length levels -- number of levels
    levelSize n = let (s, _, _) = levels VB.! n in s

    q          = FL.fold FL.sum $ fmap colsForLevel levels -- total number of columns in Z
    -- build list of items to put in Z, level by level
    obsIndices = Seq.fromFunction nO id
    intercept startingCol level =
      let spIndex obsIndex = (obsIndex, rc obsIndex VB.! level)
          spEntry obsIndex =
            let spi = spIndex obsIndex in (fst spi, startingCol + snd spi, 1)
      in  fmap spEntry obsIndices
    slope startingCol level predIndex =
      let spIndex obsIndex = (obsIndex, rc obsIndex VB.! level)
          spValue obsIndex = mX `LA.atIndex` (obsIndex, predIndex)
          spEntry obsIndex =
            let spi = spIndex obsIndex
            in  (fst spi, startingCol + snd spi, spValue obsIndex)
      in  fmap spEntry obsIndices
    zEntriesForSlope startingCol level vb =
      let predIndices = fmap fst $ VB.filter snd $ VB.zip
            (VB.generate (VB.length vb) id)
            vb
          startingCols = VB.generate
            (VB.length predIndices)
            (\n -> startingCol + n * (levelSize level))
          newStart = startingCol + (VB.length predIndices * (levelSize level))
          entriesF (pi, sc) = slope sc level pi
      in  ( FL.fold FL.mconcat $ fmap entriesF $ VB.zip predIndices startingCols
          , newStart
          )
    zEntriesForLevel startingCol level (qL, b, vbM) =
      let (interceptZs, newStartI) = if b
            then (intercept startingCol level, startingCol + qL)
            else (Seq.empty, startingCol)
          (slopeZs, newStartS) = case vbM of
            Nothing -> (Seq.empty, newStartI)
            Just vb -> zEntriesForSlope newStartI level vb
      in  (interceptZs <> slopeZs, newStartS)
    zFold = FL.Fold
      (\(zs, sc, l) x ->
        let (newZs, newStart) = zEntriesForLevel sc l x
        in  (zs <> newZs, newStart, l + 1)
      )
      (Seq.empty, 0, 0)
      id
    (zEntries, numCols, _) = FL.fold zFold levels
  in
    SLA.fromListSM (nO, q) zEntries

asDense sp =
  LA.matrix (snd $ SLA.dimSM sp) $ fmap (\(_, _, x) -> x) $ SLA.toDenseListSM sp

toSparseMatrix
  :: (LA.Container LA.Vector a, RealFrac a) => LA.Matrix a -> SLA.SpMatrix a
toSparseMatrix x
  = let
      (r, c)     = LA.size x
      colIndices = take c $ L.iterate (+ 1) 0
      rowIndices = take r $ L.iterate (+ 1) 0
      indexRow rI =
        fmap (\(cI, val) -> (rI, cI, val)) . zip colIndices . LA.toList
      indexedRows =
        concat . fmap (\(rI, rV) -> indexRow rI rV) . zip rowIndices . LA.toRows
    in
      SLA.fromListSM (r, c) $ indexedRows x

toSparseVector
  :: (LA.Container LA.Vector a, RealFrac a) => LA.Vector a -> SLA.SpVector a
toSparseVector v = SLA.fromListDenseSV (LA.size v) $ LA.toList v

checkProblem
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double, SemC r)
  => LA.Matrix a -- ^ X
  -> LA.Vector a -- ^ y
  -> SLA.SpMatrix a -- ^ Z
  -> P.Sem r ()
checkProblem mX vY smZ = do
  let (n, p)     = LA.size mX
      yRows      = LA.size vY
      (zRows, q) = SLA.dim smZ
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
makeA
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double)
  => LA.Matrix a -- ^ X
  -> LA.Vector a -- ^ y
  -> SLA.SpMatrix a -- ^ Z 
  -> SLA.SpMatrix a
makeA mX vY smZ =
  let (n, p)     = LA.size mX
      yRows      = LA.size vY
      (zRows, q) = SLA.dim smZ
      smX        = toSparseMatrix mX
      svY        = toSparseVector vY
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



-- S is the diagonal matrix of covariances in theta
-- T is unit-lower-triangular of off-diagonal covariances
makeSTF
  :: (Num a, VS.Storable a)
  => Levels
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a))
makeSTF levels
  = let
      f
        :: Num a
        => LevelSpec
        -> ([a] -> [(Int, Int, a)], [a] -> [(Int, Int, a)], Int, Int)
      f l@(n, _, _) =
        let
          e = effectsForLevel l
          s' th =
            fmap (\(x, v) -> (x, x, v)) $ zip (take e $ L.iterate (+ 1) 0) th
          tlt' th = fmap (\((r, c), v) -> (r, c, v)) $ zip
            ([ (r, c) | c <- [0 .. (e - 1)], r <- [(c + 1) .. (e - 1)] ])
            th
          t' th =
            tlt' th ++ fmap (\x -> (x, x, 1)) (L.take e $ L.iterate (+ 1) 0)
        in
          (s', t', e, n)
      makers = fmap f levels
      diagOffset o = fmap (\(r, c, v) -> (r + o, c + o, v))
      diagCopiesFrom i e n ivs =
        mconcat $ fmap (flip diagOffset ivs) $ L.take n $ L.iterate (+ e) i
    in
      \thV ->
        let thL = LA.toList thV
            fld = FL.Fold
              (\(sL, tL, thL', offset) (s', t', e, n) ->
                ( sL ++ diagCopiesFrom offset e n (s' thL')
                , tL ++ diagCopiesFrom offset e n (t' $ L.drop e thL')
                , L.drop (e + (e * (e - 1) `div` 2)) thL'
                , offset + (e * n)
                )
              )
              ([], [], thL, 0)
              id
            (s, t, _, qT) = FL.fold fld makers
        in  (SLA.fromListSM (qT, qT) s, SLA.fromListSM (qT, qT) t)


makeAStar
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double)
  => Int -- ^ p (cols of X)
  -> Int -- ^ q (cols of Z) 
  -> SLA.SpMatrix a -- ^ A
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a))
  -> LA.Vector a -- ^ theta
  -> SLA.SpMatrix a
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
  :: (LA.Container LA.Vector a, RealFrac a, SemC r, a ~ Double)
  => LA.Matrix a  -- ^ X
  -> LA.Vector a -- ^ y
  -> SLA.SpMatrix a -- ^ Z
  -> (LA.Vector a -> P.Sem r (SLA.SpMatrix a, SLA.SpMatrix a))
  -> LA.Vector a -- ^ theta
  -> P.Sem r (SLA.SpMatrix a)
makeAStar' mX vY smZ mkST vTh = do
  (smS, smT) <- mkST vTh
  let smX    = toSparseMatrix mX
      svY    = toSparseVector vY
      (_, q) = SLA.dim smZ
      smZS   = smZ ## smT ## smS
  liftIO $ do
    putStrLn "Z*="
    LA.disp 2 $ asDense smZS
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

data DevianceType = ML | REML deriving (Show, Eq)

-- since A* is (p + q + 1) x (p + q + 1) so is L
profiledDeviance
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double)
  => DevianceType
  -> Int -- ^ p (cols of X)
  -> Int -- ^ q (cols of Z)
  -> Int -- ^ n (rows of X or Z)
  -> SLA.SpMatrix a -- ^ A
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a))
  -> LA.Vector a -- ^ theta
  -> (a, a, a, LA.Matrix a)
profiledDeviance dt p q n smA mkST th
  = let
      smAS  = makeAStar p q smA mkST th
      cholL = LA.tr $ LA.chol $ LA.trustSym $ asDense smAS
      getDiag :: Int -> Double
      getDiag n = cholL `LA.atIndex` (n, n)
      r = getDiag $ p + q
      logDiagF n = FL.fold (FL.premap (log . getDiag) FL.sum) [0 .. (n - 1)]
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
thetaLowerBounds levels = NL.LowerBounds $ FL.fold fld levels
 where
  fld = FL.Fold
    (\bs l ->
      let e = effectsForLevel l
      in  bs ++ L.replicate e 0 ++ replicate (e * (e - 1) `div` 2) (negate 1e5)
    )
    []
    LA.fromList

minimizeDeviance
  :: (LA.Container LA.Vector a, RealFrac a, SemC r, a ~ Double)
  => DevianceType
  -> Int -- ^ p (cols of X)
  -> Int -- ^ q (cols of Z)
  -> Int -- ^ n (rows of X or Z)
  -> Levels
  -> SLA.SpMatrix a -- ^ A
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a))
  -> LA.Vector a -- ^ initial guess for theta
  -> P.Sem r (LA.Vector a, (a, a, a, LA.Matrix a))
minimizeDeviance dt p q n levels smA mkST th0 = do
  let
    pd x = profiledDeviance dt p q n smA mkST x
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

parametersFromSolution
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double)
  => DevianceType -- ^ unused for now??
  -> Int
  -> Int
  -> (LA.Vector a -> (SLA.SpMatrix a, SLA.SpMatrix a)) -- ^ makeST
  -> LA.Vector a -- ^ solution theta
  -> LA.Matrix a -- ^ solution L
  -> a -- ^ solution r
  -> (LA.Vector a, LA.Vector a, LA.Vector a) -- ^ beta, b and b* (also called u)                       
parametersFromSolution dt p q makeST th ldL r =
  let v        = LA.fromList $ L.replicate (p + q) 0 ++ [r]
      solution = LA.takesV [q, p] $ LA.flatten $ LA.triSolve LA.Upper
                                                             (LA.tr ldL)
                                                             (LA.asColumn v)
      beta   = solution L.!! 1
      bStar  = solution L.!! 0
      (s, t) = makeST th
      b      = ((asDense t) LA.<> (asDense s)) LA.#> bStar
  in  (beta, b, bStar)

report
  :: (LA.Container LA.Vector a, RealFrac a, a ~ Double, SemC r)
  => Int -- ^ p
  -> Int -- ^ q
  -> Levels
  -> LA.Vector a -- ^ y
  -> LA.Matrix a -- ^ X
  -> SLA.SpMatrix a -- ^ Z
  -> LA.Vector a -- ^ beta
  -> LA.Vector a -- ^ b
  -> P.Sem r ()
report p q levels vY mX smZ vBeta vb = do
  let
    mZ = asDense smZ
    reportMean0 prefix v = do
      let mean = meanV v
          var  = (v LA.<.> v) / (realToFrac $ LA.size v)
      putStrLn $ (T.unpack prefix) ++ ": mean (should be 0)=" ++ show mean
      putStrLn
        $  (T.unpack prefix)
        ++ ": variance (assuming mean 0)="
        ++ show var
      putStrLn $ (T.unpack prefix) ++ ": std. dev (assuming mean 0)=" ++ show
        (sqrt var)

    vEps = vY - ((mX LA.#> vBeta) + (mZ LA.#> vb))
    meanV v = LA.sumElements v / realToFrac (LA.size v)
    mEps   = meanV vEps --LA.sumElements vEps / realToFrac (LA.size vEps)
    varEps = vEps LA.<.> vEps -- assumes mEps is 0.  Which it should be!!                
    levelReport l@(n, b, _) b' = do
      putStrLn $ show n ++ " groups"
      when b $ reportMean0 "Intercept" $ VS.take n b'
      let (numSlopes, bS) = case b of
            True  -> (effectsForLevel l - 1, VS.drop n b')
            False -> (effectsForLevel l, b')
      mapM_
        (\s -> reportMean0 ("Slope " <> (T.pack $ show s))
                           (VS.take n $ VS.drop (n * s) bS)
        )
        [0 .. (numSlopes - 1)]
  liftIO $ reportMean0 "Residual" vEps
  let numberedLevels = zip [0 ..] (VB.toList levels)
  liftIO $ mapM_
    (\(lN, l) -> putStrLn ("Level " ++ show lN) >> levelReport l vb)
    numberedLevels
