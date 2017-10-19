{-# OPTIONS_GHC -fprint-explicit-kinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}

module Math.Linear.Vector where

import Foundation
import Foundation.Array
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Class.Storable
import Foundation.Collection
import Foundation.Primitive

import GHC.TypeLits

import Prelude (fromIntegral)
import Control.Monad.ST (RealWorld)
import Foreign.Marshal.Alloc (alloca)
import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import Math.Linear.ElemWise

newtype Vec a (n :: Nat) = V { vect :: UArray a } deriving (Eq, Ord, Show)

newtype MVec a (n :: Nat) s = MV { vect :: MUArray a s }

type IOVec a (n :: Nat) = MVec a n RealWorld

type instance Element (Vec a n) = a

type instance Element (MVec a n s) = a

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

withVPtr' :: forall monad a b n. (PrimMonad monad, PrimType a, KnownNat n) => Vec a n -> (Ptr a -> Int32 -> monad b) -> monad b
withVPtr' V{..} f = withPtr vect $ \p -> f p (fromIntegral $ natVal (Proxy @n))

withMutableVPtr :: (PrimMonad monad, PrimType a) => MVec a n (PrimState monad) -> (Ptr a -> monad b) -> monad b
withMutableVPtr MV{..} = withMutablePtr vect

withMutableVPtr' :: forall monad a b n. (PrimMonad monad, PrimType a, KnownNat n) => MVec a n (PrimState monad) -> (Ptr a -> Int32 -> monad b) -> monad b
withMutableVPtr' MV{..} f = withMutablePtr vect $ \p -> f p (fromIntegral $ natVal (Proxy @n))

instance (I.Elem a, KnownNat n) => ElemWise (Vec a n) where
    -- * tensor and scalar
    shift x v = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' column ->
            withVPtr v $ \pv ->
                alloca $ \p -> do
                    poke p x
                    I.call $ I.shift pv' p pv 1 column
        unsafeFreeze v'
    times x v = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' column ->
            withVPtr v $ \pv ->
                alloca $ \p -> do
                    poke p x
                    I.call $ I.times pv' p pv 1 column
        unsafeFreeze v'
    -- * negative
    negative v = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v $ \pv ->
                I.call $ I.negative pv' pv 1 nlen
        unsafeFreeze v'
    -- * arithmetic
    add v1 v2 = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v1 $ \pv1 ->
                withVPtr v2 $ \pv2 ->
                    I.call $ I.add pv' 1 nlen 1 pv1 pv2
        unsafeFreeze v'
    minus v1 v2 = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v1 $ \pv1 ->
                withVPtr v2 $ \pv2 ->
                    I.call $ I.minus pv' 1 nlen 1 pv1 pv2
        unsafeFreeze v'
    mult v1 v2 = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v1 $ \pv1 ->
                withVPtr v2 $ \pv2 ->
                    I.call $ I.mult pv' 1 nlen 1 pv1 pv2
        unsafeFreeze v'
    division v1 v2 = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v1 $ \pv1 ->
                withVPtr v2 $ \pv2 ->
                    I.call $ I.division pv' 1 nlen 1 pv1 pv2
        unsafeFreeze v'
    -- * data generation
    constreplic x = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            alloca $ \p -> do
                poke p x
                I.call $ I.replicate pv' p nlen
        unsafeFreeze v'
    -- * extensions
    logistic v = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v $ \pv ->
                I.call $ I.logistic pv' pv 1 nlen
        unsafeFreeze v'
    logisticd v = unsafePerformIO $ do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            withVPtr v $ \pv ->
                I.call $ I.logisticd pv' pv 1 nlen
        unsafeFreeze v'

instance (I.Elem a, KnownNat n) => MutElemWise (MVec a n RealWorld) where
    -- * tensor and scalar
    shift' x v = do
        withMutableVPtr' v $ \pv nlen ->
            alloca $ \p -> do
                poke p x
                I.call $ I.shift pv p pv 1 nlen
        return ()
    times' x v = do
        withMutableVPtr' v $ \pv nlen ->
            alloca $ \p -> do
                poke p x
                I.call $ I.times pv p pv 1 nlen
        return ()
    -- * negative
    negative' v = do
        withMutableVPtr' v $ \pv nlen ->
            I.call $ I.negative pv pv 1 nlen
        return ()
    -- -- * arithmetic
    add' v1 v2 = do
        withMutableVPtr' v1 $ \pv1 nlen ->
            withMutableVPtr v2 $ \pv2 ->
                I.call $ I.add pv1 1 nlen 1 pv1 pv2
        return ()
    minus' v1 v2 = do
        withMutableVPtr' v1 $ \pv1 nlen ->
            withMutableVPtr v2 $ \pv2 ->
                I.call $ I.minus pv1 1 nlen 1 pv1 pv2
        return ()
    mult' v1 v2 = do
        withMutableVPtr' v1 $ \pv1 nlen ->
            withMutableVPtr v2 $ \pv2 ->
                I.call $ I.mult pv1 1 nlen 1 pv1 pv2
        return ()
    division' v1 v2 = do
        withMutableVPtr' v1 $ \pv1 nlen ->
            withMutableVPtr v2 $ \pv2 ->
                I.call $ I.division pv1 1 nlen 1 pv1 pv2
        return ()
    -- * data generation
    constreplic' x = do
        v' <- mutNew undefined
        withMutableVPtr' v' $ \pv' nlen ->
            alloca $ \p -> do
                poke p x
                I.call $ I.replicate pv' p nlen
        return v'
    -- * extensions
    logistic' v = do
        withMutableVPtr' v $ \pv nlen ->
            I.call $ I.logistic pv pv 1 nlen
        return ()
    logisticd' v = do
        withMutableVPtr' v $ \pv nlen ->
            I.call $ I.logisticd pv pv 1 nlen
        return ()
