{-# LANGUAGE AllowAmbiguousTypes         #-}
{-# LANGUAGE FlexibleContexts            #-}
{-# LANGUAGE GADTs                       #-}
{-# LANGUAGE GeneralizedNewtypeDeriving  #-}
{-# LANGUAGE RankNTypes                  #-}
{-# LANGUAGE ScopedTypeVariables         #-}
{-# LANGUAGE TypeApplications            #-}
{-# LANGUAGE UndecidableInstances        #-}

module Polysemy.ConstraintAbsorber.MonadCatch
  (
    -- * Constraint Absorbers
    absorbMonadThrow
  , absorbMonadCatch
    -- * Re-exports
  , Exception(..)
  , SomeException
  )
where

import qualified Control.Monad.Catch           as C
import           Control.Monad.Catch            ( Exception(..)
                                                , SomeException
                                                )
import           Polysemy
import           Polysemy.ConstraintAbsorber
import qualified Polysemy.Error                as E


------------------------------------------------------------------------------
-- | Introduce a local 'S.MonadCatch' constraint on 'Sem' --- allowing it to
-- interop nicely with exceptions
--
-- @since 0.5.0.0
absorbMonadCatch
  :: Member (E.Error C.SomeException) r
  => (C.MonadCatch (Sem r) => Sem r a)
       -- ^ A computation that requires an instance of 'C.MonadCatch'
       -- or 'C.MonadThrow' for
       -- 'Sem'. This might be something with type @'C.MonadCatch' e m => m a@.
  -> Sem r a
absorbMonadCatch =
  absorbWithSem @C.MonadCatch @Action (CatchDict E.throw E.catch) (Sub Dict)
{-# INLINABLE absorbMonadCatch #-}

absorbMonadThrow
  :: Member (E.Error C.SomeException) r
  => (C.MonadThrow (Sem r) => Sem r a)
       -- ^ A computation that requires an instance of 'C.MonadCatch'
       -- or 'C.MonadThrow' for
       -- 'Sem'. This might be something with type @'C.MonadCatch' e m => m a@.
  -> Sem r a
absorbMonadThrow = absorbMonadCatch
{-# INLINABLE absorbMonadThrow #-}

------------------------------------------------------------------------------
-- | A dictionary of the functions we need to supply
-- to make an instance of Error
data CatchDict m = CatchDict
  { throwM_ :: forall a. C.SomeException -> m a
  , catch_ :: forall a. m a -> (C.SomeException -> m a) -> m a
  }


------------------------------------------------------------------------------
-- | Wrapper for a monadic action with phantom
-- type parameter for reflection.
-- Locally defined so that the instance we are going
-- to build with reflection must be coherent, that is
-- there cannot be orphans.
newtype Action m s' a = Action { action :: m a }
  deriving (Functor, Applicative, Monad)


------------------------------------------------------------------------------
-- | Given a reifiable mtl Error dictionary,
-- we can make an instance of @MonadError@ for the action
-- wrapped in @Action@.
instance ( Monad m
         , Reifies s' (CatchDict m)
         ) => C.MonadThrow (Action m s') where
  throwM e = Action $ throwM_ (reflect $ Proxy @s') (C.toException e)
  {-# INLINEABLE throwM #-}

instance ( Monad m
         , Reifies s' (CatchDict m)
         )  => C.MonadCatch (Action m s') where
  catch x f =
    let catchF = catch_ (reflect $ Proxy @s')
    in Action $ (action x) `catchF` \e -> case C.fromException e of
      Just e' -> action $ f e'
      _ -> throwM_ (reflect $ Proxy @s') (C.toException e)
  {-# INLINEABLE catch #-}


