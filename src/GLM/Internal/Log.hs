{-# LANGUAGE DataKinds        #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}
{-# LANGUAGE LambdaCase       #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators    #-}

module GLM.Internal.Log where

import qualified Data.Text                     as T
--import           Control.Monad.IO.Class         ( MonadIO(liftIO) )

import qualified Polysemy                      as P
import qualified Polysemy.Output               as P
import qualified Polysemy.Trace                as P

log :: P.Member P.Trace r => T.Text -> P.Sem r ()
log = P.trace . T.unpack

runLogIO :: P.Member (P.Lift IO) r => P.Sem (P.Trace ': r) a -> P.Sem r a
runLogIO = P.runTraceIO
