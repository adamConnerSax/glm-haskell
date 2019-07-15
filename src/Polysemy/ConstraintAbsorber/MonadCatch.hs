{-# LANGUAGE AllowAmbiguousTypes         #-}
{-# LANGUAGE FlexibleContexts            #-}
{-# LANGUAGE GADTs                       #-}
{-# LANGUAGE GeneralizedNewtypeDeriving  #-}
{-# LANGUAGE MultiParamTypeClasses       #-}
{-# LANGUAGE RankNTypes                  #-}
{-# LANGUAGE ScopedTypeVariables         #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances        #-}

{-# LANGUAGE QuantifiedConstraints       #-}
module Polysemy.ConstraintAbsorber.MonadCatch
  ( absorbMonadCatch
  )
where

import qualified Control.Monad.Catch           as C
import           Polysemy
import           Polysemy.ConstraintAbsorber
import qualified Polysemy.Error                as E


------------------------------------------------------------------------------
-- | Introduce a local 'S.MonadCatch' constraint on 'Sem' --- allowing it to
-- interop nicely with exceptions
--
-- @since 0.5.0.0
absorbMonadCatch
  :: forall r a
   . (forall e. (Member (E.Error e) r
                ,C.Exception e)
     )
  => (C.MonadCatch (Sem r) => Sem r a)
       -- ^ A computation that requires an instance of 'C.MonadCatch'
       -- or 'C.MonadThrow' for
       -- 'Sem'. This might be something with type @'C.MonadCatch' e m => m a@.
  -> Sem r a
absorbMonadCatch = absorbWithSem @C.MonadCatch @Action
  (ThrowDict (E.throw) (E.catch))
  (Sub Dict)
{-# INLINABLE absorbMonadCatch #-}


------------------------------------------------------------------------------
-- | A dictionary of the functions we need to supply
-- to make an instance of Error
data ThrowDict m = ThrowDict
  { throwM_ :: forall a e. e -> m a
  , catch_ :: forall a e. m a -> (e -> m a) -> m a
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
         , forall e. (C.Exception e
                     , Reifies s' (ThrowDict m)
                     )
         ) => C.MonadThrow (Action m s') where
  throwM e = Action $ throwM_ (reflect $ Proxy @s') e
  {-# INLINEABLE throwM #-}

instance ( Monad m
         , forall e. (C.Exception e
                     , Reifies s' (ThrowDict m)
                     )
         ) => C.MonadCatch (Action m s') where
  catch x f = Action $ catch_ (reflect $ Proxy @s') (action x) (action . f)
  {-# INLINEABLE catch #-}

