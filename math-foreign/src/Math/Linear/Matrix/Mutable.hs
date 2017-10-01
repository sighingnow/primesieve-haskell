-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Matrix.Mutable
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Multiple dimensions matrices repersentation in Haskell.
--
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeSynonymInstances #-}

module Math.Linear.Matrix.Mutable where
  -- ( MMat(..)
  -- , IOMat
  -- , STMat
  -- , null
  -- , square
  -- , new
  -- , zeros
  -- , ones
  -- , identity
  -- , random
  -- , replicate
  -- , (!)
  -- , shift
  -- , times
  -- , read
  -- , write
  -- , copy
  -- , modify
  -- , unsafeRead
  -- , unsafeWrite
  -- , unsafeCopy
  -- , unsafeModify
  -- , unsafeWith)
  -- where

import Foundation
import Foundation.Array
import Foundation.Array.Internal (withMutablePtr)
import Foundation.Class.Storable
import Foundation.Collection
import Foundation.Primitive

import GHC.Num (Num)
import Control.Monad.ST (RealWorld)
import Foreign.Marshal.Alloc (alloca)

import qualified Math.Linear.Internal as I

data MMat a s = M
    { row :: {-# UNPACK #-}!Int -- ^ rows
    , column :: {-# UNPACK #-}!Int -- ^ columns
    , vect :: MUArray a s -- ^ data in plain vector.
    }

type IOMat a = MMat a RealWorld

type STMat a s = MMat a s

-- | If the matrix is empty.
null :: MMat a s -> Bool
null M{..} = row == 0 && column == 0

{-# INLINE null #-}

-- | If the matrix is a square matrix.
square :: MMat a s -> Bool
square M{..} = row == column

{-# INLINE square #-}

-- | Construct a new matrix without initialisation.
new :: (PrimMonad m, PrimType a)
    => Int -> Int -> m (MMat a (PrimState m))
new r c = M r c <$> mutNew (CountOf (r * c))

{-# INLINE new #-}

-- | Construct a matrix with all zeros.
zeros :: (PrimMonad m, PrimType a, Num a)
    => Int -> Int -> m (MMat a (PrimState m))
zeros r c = replicate' r c 0

{-# INLINE zeros #-}

-- | Construct a matrix with all ones.
ones :: (PrimMonad m, PrimType a, Num a)
    => Int -> Int -> m (MMat a (PrimState m))
ones r c = replicate' r c 1

{-# INLINE ones #-}

-- | Construct a identity matrix, square is not required.
identity :: I.Elem a
    => Int -> Int -> IO (IOMat a)
identity r c = do
    m <- new r c
    unsafeWith m $ \xs r' c' -> I.call $ I.identity xs r' c'
    return m

{-# INLINE identity #-}

-- | Construct a identity matrix, square is not required.
random :: I.Elem a
    => Int -> Int -> IO (IOMat a)
random r c = do
    m <- new r c
    unsafeWith m $ \xs r' c' -> I.call' $ I.random_ xs r' c'
    return m

{-# INLINE random #-}

-- | Construct a matrix with all given constants.
replicate' :: (PrimMonad m, PrimType a)
    => Int -> Int -> a -> m (MMat a (PrimState m))
replicate' r c v = M r c <$> unsafeThaw (replicate (CountOf (r * c)) v)

{-# INLINE replicate' #-}

-- | Get specified element from matrix unsafely.
(!) :: (PrimMonad m, PrimType a)
    => MMat a (PrimState m) -> (Int, Int) -> m a
(!) m = uncurry (unsafeRead m)

{-# INLINE (!) #-}

-- | Add scalar value to every elements in matrix.
shift :: I.Elem a
    => a -> IOMat a -> IO ()
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
times :: I.Elem a
    => a -> IOMat a -> IO ()
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
read :: (PrimMonad m, PrimType a)
    => MMat a (PrimState m) -> Int -> Int -> m a
read M{..} r c = mutRead vect (Offset (r * column + c))

{-# INLINE read #-}

-- | Write value to matrix.
write :: (PrimMonad m, PrimType a)
    => MMat a (PrimState m) -> Int -> Int -> a -> m ()
write M{..} r c = mutWrite vect (Offset (r * column + c))

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

unsafeRead :: (PrimMonad m, PrimType a)
    => MMat a (PrimState m) -> Int -> Int -> m a
unsafeRead M{..} r c = mutUnsafeRead vect (Offset (r * column + c))

{-# INLINE unsafeRead #-}

unsafeWrite :: (PrimMonad m, PrimType a)
    => MMat a (PrimState m) -> Int -> Int -> a -> m ()
unsafeWrite M{..} r c = mutUnsafeWrite vect (Offset (r * column + c))

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

unsafeWith :: (PrimMonad m, PrimType a)
    => MMat a (PrimState m) -> (Ptr a -> Int -> Int -> m b) -> m b
unsafeWith M{..} f = withMutablePtr vect $ \p -> f p row column

{-# INLINE unsafeWith #-}
