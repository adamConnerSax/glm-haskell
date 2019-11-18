{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE Rank2Types          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
module Numeric.GLM.Bootstrap
  ( module Numeric.GLM.Bootstrap
  , module Polysemy.RandomFu
  , module Polysemy.Async
  )
where

import qualified Numeric.LinearAlgebra         as LA
import           Data.Maybe                     ( catMaybes )
import qualified Data.Text                     as T
import qualified Data.Vector.Storable          as VS
import qualified Data.Vector                   as VB

import           Control.Monad                  ( when )
import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Numeric.GLM.FunctionFamily    as GLM
import qualified Numeric.GLM.ModelTypes        as GLM
import qualified Numeric.GLM.MixedModel        as GLM
import qualified Numeric.LinearAlgebra.CHOLMOD.CholmodExtras
                                               as CH


import qualified Data.Random                   as DR
import qualified Data.Random.Distribution.Normal
                                               as DR
import qualified Data.Random.Distribution.Bernoulli
                                               as DR
import qualified Data.Random.Distribution.Binomial
                                               as DR
import qualified Data.Random.Distribution.Poisson
                                               as DR
import qualified Data.Random.Distribution.Gamma
                                               as DR

import qualified Knit.Effect.Logger            as P
import qualified Polysemy                      as P
import qualified Polysemy.Error                as P
import qualified Polysemy.RandomFu             as P
import           Polysemy.RandomFu              ( runRandomIO )
import qualified Polysemy.Async                as P
import           Polysemy.Async                 ( asyncToIO
                                                , asyncToIOFinal
                                                )
--import qualified Polysemy.ConstraintAbsorber.MonadRandom
--                                               as P

import           Control.Concurrent            as CC
import           Control.Concurrent.MVar       as CC

-- TODO: handle the sometimes requiring dev parameter in a cleaner way.
generateSample
  :: P.Member P.RandomFu r
  => GLM.ObservationDistribution
  -> Double
  -> Double
  -> P.Sem r Double
generateSample od mu dev = P.sampleRVar d
 where
  d = case od of
    GLM.DNormal    -> DR.normal mu dev
    GLM.DBernoulli -> DR.bernoulli mu
    GLM.DBinomial n ->
      fmap (\x -> realToFrac x / realToFrac n) $ DR.binomial n mu -- binomial observations are as fractions
    GLM.DPoisson -> DR.poisson mu
    GLM.DGamma   -> DR.gamma (1 / dev) (mu * dev) -- Shape & Scale as in deviance


-- NB: assumes the same deviance across all when deviance matters.  Not sure this is always correct. 
generateSamples
  :: P.Member P.RandomFu r
  => GLM.ObservationsDistribution
  -> VS.Vector Double
  -> Double
  -> P.Sem r (VS.Vector Double)
generateSamples osd vMu dev = do
  let n = VS.length vMu
      gs od mu = generateSample od mu dev
      vOD = case osd of
        GLM.Normal      -> VB.replicate n GLM.DNormal
        GLM.Bernoulli   -> VB.replicate n GLM.DBernoulli
        GLM.Binomial vN -> VB.map GLM.DBinomial (VB.convert vN) -- these are boxed
        GLM.Poisson     -> VB.replicate n GLM.DPoisson
        GLM.Gamma       -> VB.replicate n GLM.DGamma
  VS.convert <$> VB.zipWithM gs vOD (VB.convert vMu)


parametricBootstrap
  :: (GLM.EffectsIO r, P.Member P.RandomFu r, P.Member P.Async r)
  => GLM.MinimizeDevianceVerbosity
  -> GLM.DevianceType
  -> GLM.MixedModel b g -- we'll modify this, substuting bootstrapped data
  -> GLM.RandomEffectCalculated
  -> GLM.CholmodFactor
  -> GLM.CovarianceVec -- should be solution from a minimize run
  -> VS.Vector Double -- mu from solution
  -> Double -- deviance from solution, may be ignored
  -> Int -- number of times to bootstrap
  -> Bool
  -> P.Sem r [(GLM.BetaU, VS.Vector Double)] -- beta, u and b
parametricBootstrap mdv dt mm0 reCalc cf thSol vMuSol devSol n doConcurrently =
  do
    let (fpC, fpF, smP) = cf
        generateOne =
          generateSamples (GLM.observationDistribution mm0) vMuSol devSol
        solveOne (n, newMM) =
          P.wrapPrefix ("Bootstrap (" <> (T.pack $ show n) <> ")") $ do
            P.logLE P.Diagnostic "Starting..."
            fpF_Copy <- if doConcurrently
              then liftIO $ CH.copyFactor fpF fpC
              else return fpF
            (_, _, _, betaU, b, _) <- GLM.minimizeDevianceInner
              mdv
              dt
              newMM
              reCalc
              (fpC, fpF_Copy, smP)
              thSol
            P.logLE P.Diagnostic "Finished."
            return (betaU, b)
    newObservations <- mapM (const generateOne) $ replicate n ()
    let newMixedModels =
          fmap (\os -> GLM.changeMMObservations os mm0) newObservations
        seqF = if doConcurrently
          then fmap catMaybes . P.sequenceConcurrently
          else sequence
    res <- seqF $ fmap solveOne $ zip [1 ..] newMixedModels
    let m = length res
    when (m < n)
      $  P.throw
      $  GLM.OtherGLMError
      $  "parametric bootstrap attempted "
      <> (T.pack $ show n)
      <> " refits but only "
      <> (T.pack $ show m)
      <> " succeeded."
    return res

{-  
-- launch each action on its own thread, writing the result to an MVar.
-- as each MVar is written to,
-- place each result in the appropriate place in the structure.
-- return the new structure when all results are available.
sequenceConcurrently :: Traversable t => t (IO a) -> IO (t a)
sequenceConcurrently actions = do
  let f :: IO a -> IO (CC.MVar a, CC.ThreadId)
      f action = do
        mvar     <- CC.newEmptyMVar
        threadId <- forkIO $ (action >>= CC.putMVar mvar)
        return (mvar, threadId)
      g :: (CC.MVar a, CC.ThreadId) -> IO a
      g (mvar, _) = takeMVar mvar
  forked <- traverse f actions -- IO (t (MVar a, ThreadId))
  traverse g forked
-}
