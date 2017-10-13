{-# OPTIONS_GHC -fprint-explicit-kinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}

module Math.Linear.Vector.Dependent where

import Foundation
import Foundation.Array
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Class.Storable
import Foundation.Collection
import Foundation.Primitive

import GHC.TypeLits

import Prelude (fromIntegral)

newtype Vec a (n :: Nat) = V { vect :: UArray a } deriving (Eq, Ord, Show)

newtype MVec a (n :: Nat) s = MV { vect :: MUArray a s }

type instance Element (Vec a n) = a

instance PrimType a => IsList (Vec a n) where
  type Item (Vec a n) = a
  fromList = V . fromList
  fromListN n xs = V $ fromListN n xs
  toList V{..} = toList vect

deriving instance PrimType a => Monoid (Vec a n)

deriving instance NormalForm (Vec a n)

deriving instance PrimType a => Fold1able (Vec a n)

deriving instance PrimType a => Foldable (Vec a n)

deriving instance PrimType a => IndexedCollection (Vec a n)

deriving instance PrimType a => InnerFunctor (Vec a n)

deriving instance PrimType a => Copy (Vec a n)

deriving instance PrimType a => Collection (Vec a n)

deriving instance PrimType a => Sequential (Vec a n)

deriving instance PrimType a => Zippable (Vec a n)

instance (PrimType a, KnownNat n) => MutableCollection (MVec a n) where
  type MutableFreezed (MVec a n) = Vec a n
  type MutableKey (MVec a n) = Offset a
  type MutableValue (MVec a n) = a
  unsafeThaw V{..} = MV <$> unsafeThaw vect
  unsafeFreeze MV{..} = V <$> unsafeFreeze vect
  thaw V{..} = MV <$> thaw vect
  freeze MV{..} = V <$> freeze vect
  -- The given size argument would be ignored.
  mutNew ~_ = MV <$> mutNew nlen where nlen = CountOf (fromIntegral (natVal (Proxy :: Proxy n)))
  mutUnsafeWrite MV{..} loc = mutUnsafeWrite vect loc
  mutWrite MV{..} loc = mutWrite vect loc
  mutUnsafeRead MV{..} loc = mutUnsafeRead vect loc
  mutRead MV{..} loc = mutRead vect loc

withVPtr :: (PrimMonad monad, PrimType a) => Vec a n -> (Ptr a -> monad b) -> monad b
withVPtr V{..} = withPtr vect

withMutableVPtr :: (PrimMonad monad, PrimType a) => MVec a n (PrimState monad) -> (Ptr a -> monad b) -> monad b
withMutableVPtr MV{..} = withMutablePtr vect
