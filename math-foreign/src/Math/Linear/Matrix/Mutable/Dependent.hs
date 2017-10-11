-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Matrix.Mutable.Dependent
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Multiple dimensions matrices repersentation in Haskell.
--
{-# OPTIONS_GHC -fprint-explicit-kinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}

module Math.Linear.Matrix.Mutable.Dependent where

import Foundation
import Foundation.Array
import Foundation.Array.Internal (withMutablePtr)
import Foundation.Class.Storable
import Foundation.Collection
import Foundation.Primitive

import GHC.TypeLits
import GHC.Num (Num)
import Control.Monad.ST (RealWorld)
import Foreign.Marshal.Alloc (alloca)

import qualified Math.Linear.Internal as I

data MMat a (m :: Nat) (n :: Nat) s = M
    { vect :: MUArray a s -- ^ data in plain vector.
    }

type IOMat a m n = MMat a m n RealWorld

type STMat a m n s = MMat a m n s

-- | If the matrix is empty.
null :: forall a m n s. (KnownNat m, KnownNat n)
    => MMat a m n s -> Bool
null _ = isJust (sameNat (Proxy :: Proxy m) (Proxy :: Proxy 0)) && isJust (sameNat (Proxy :: Proxy n) (Proxy :: Proxy 0))

{-# INLINE null #-}

-- | If the matrix is a square matrix.
square :: forall a m n s. (KnownNat m, KnownNat n)
    => MMat a m n s -> Bool
square _ = isJust (sameNat (Proxy :: Proxy m) (Proxy :: Proxy n))

{-# INLINE square #-}

-- | Construct a new matrix without initialisation.
new :: (PrimMonad monad, PrimType a, KnownNat m, KnownNat n)
    => Proxy m -> Proxy n -> monad (MMat a m n (PrimState monad))
new m n = M <$> mutNew (integralCast (r * c))
    where r = natVal m
          c = natVal n

{-# INLINE new #-}

-- | Construct a matrix with all zeros.
zeros :: (PrimMonad monad, PrimType a, Num a, KnownNat m, KnownNat n)
    => Proxy m -> Proxy n -> monad (MMat a m n (PrimState monad))
zeros r c = replicate' r c 0

{-# INLINE zeros #-}

-- | Construct a matrix with all ones.
ones :: (PrimMonad monad, PrimType a, Num a, KnownNat m, KnownNat n)
    => Proxy m -> Proxy n -> monad (MMat a m n (PrimState monad))
ones r c = replicate' r c 1

{-# INLINE ones #-}

-- | Construct a identity matrix, square is not required.
identity :: (I.Elem a, KnownNat m, KnownNat n)
    => Proxy m -> Proxy n -> IO (IOMat a m n)
identity r c = do
    m <- new r c
    unsafeWith m $ \xs r' c' -> I.call $ I.identity xs r' c'
    return m

{-# INLINE identity #-}

-- | Construct a identity matrix, square is not required.
random :: (I.Elem a, KnownNat m, KnownNat n)
    => Proxy m -> Proxy n -> IO (IOMat a m n)
random r c = do
    m <- new r c
    unsafeWith m $ \xs r' c' -> I.call' $ I.random_ xs r' c'
    return m

{-# INLINE random #-}

-- | Construct a matrix with all given constants.
replicate' :: (PrimMonad monad, PrimType a, KnownNat m, KnownNat n)
    => Proxy m -> Proxy n -> a -> monad (MMat a m n (PrimState monad))
replicate' m n v = M <$> thaw (replicate (integralCast (r * c)) v)
    where r = natVal m
          c = natVal n

{-# INLINE replicate' #-}

-- | Get specified element from matrix unsafely.
(!) :: (PrimMonad monad, PrimType a, KnownNat m, KnownNat n)
    => MMat a m n (PrimState monad) -> (Proxy m, Proxy n) -> monad a
(!) = uncurry . unsafeRead

{-# INLINE (!) #-}

-- | Add scalar value to every elements in matrix.
shift :: (I.Elem a, KnownNat m, KnownNat n)
    => a -> IOMat a m n -> IO ()
shift x m = do
    unsafeWith m $
        \xs r c ->
             alloca $
             \p -> do
                 poke p x
                 I.call $ I.shift xs p xs r c
    return ()

{-# INLINE shift #-}

-- | Multiply scalar value to every elements in matrix.
times :: (I.Elem a, KnownNat m, KnownNat n)
    => a -> IOMat a m n -> IO ()
times x m = do
    unsafeWith m $
        \xs r c ->
             alloca $
             \p -> do
                 poke p x
                 I.call $ I.times xs p xs r c
    return ()

{-# INLINE times #-}

-- | Read value from matrix.
read :: forall monad a m n u v. (PrimMonad monad, PrimType a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, u <= m, v <= n)
    => MMat a m n (PrimState monad) -> Proxy u -> Proxy v -> monad a
read M{..} u v = mutRead vect (integralCast (r * column + c))
    where r = natVal u
          c = natVal v
          column = natVal (Proxy :: Proxy m)

{-# INLINE read #-}

-- | Write value to matrix.
write :: forall monad a m n u v. (PrimMonad monad, PrimType a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, u <= m, v <= n)
    => MMat a m n (PrimState monad) -> Proxy u -> Proxy v -> a -> monad ()
write M{..} u v = mutWrite vect (integralCast (r * column + c))
    where r = natVal u
          c = natVal v
          column = natVal (Proxy :: Proxy m)

{-# INLINE write #-}

-- -- | Copy one matrix to another.
-- copy
--     :: (PrimMonad m, I.Elem a)
--     => MMat (PrimState m) a -> MMat (PrimState m) a -> m ()
-- copy m1 m2 = V.copy (vect m1) (vect m2)

-- {-# INLINE copy #-}

-- -- | Modify element in matrix using given function.
-- modify
--     :: (PrimMonad m, I.Elem a)
--     => MMat (PrimState m) a -> (a -> a) -> Int -> Int -> m ()
-- modify M {..} f r c = V.modify vect f (r * column + c)

-- {-# INLINE modify #-}

unsafeRead :: forall monad a m n u v. (PrimMonad monad, PrimType a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, u <= m, v <= n)
    => MMat a m n (PrimState monad) -> Proxy u -> Proxy v -> monad a
unsafeRead M{..} u v = mutUnsafeRead vect (integralCast (r * column + c))
    where r = natVal u
          c = natVal v
          column = natVal (Proxy :: Proxy m)

{-# INLINE unsafeRead #-}

unsafeWrite :: forall monad a m n u v. (PrimMonad monad, PrimType a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, u <= m, v <= n)
    => MMat a m n (PrimState monad) -> Proxy u -> Proxy v -> a -> monad ()
unsafeWrite M{..} u v = mutUnsafeWrite vect (integralCast (r * column + c))
    where r = natVal u
          c = natVal v
          column = natVal (Proxy :: Proxy m)

{-# INLINE unsafeWrite #-}

-- unsafeCopy
--     :: (PrimMonad m, I.Elem a)
--     => MMat (PrimState m) a -> MMat (PrimState m) a -> m ()
-- unsafeCopy m1 m2 = V.unsafeCopy (vect m1) (vect m2)

-- {-# INLINE unsafeCopy #-}

-- unsafeModify
--     :: (PrimMonad m, I.Elem a)
--     => MMat (PrimState m) a -> (a -> a) -> Int -> Int -> m ()
-- unsafeModify M {..} f r c = V.unsafeModify vect f (r * column + c)

-- {-# INLINE unsafeModify #-}

unsafeWith :: forall monad a b m n. (PrimMonad monad, PrimType a, KnownNat m, KnownNat n)
    => MMat a m n (PrimState monad) -> (Ptr a -> Int32 -> Int32 -> monad b) -> monad b
unsafeWith M{..} f = withMutablePtr vect $ \p -> f p row column
    where row = integralDownsize $ natVal (Proxy :: Proxy m)
          column = integralDownsize $ natVal (Proxy :: Proxy n)

{-# INLINE unsafeWith #-}
