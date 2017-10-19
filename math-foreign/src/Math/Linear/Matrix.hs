-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Matrix
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Matrices representation in Haskell.
--
{-# OPTIONS_GHC -fprint-explicit-kinds #-}
{-# OPTIONS_GHC -Wno-orphans #-}
{-# OPTIONS_GHC -Wno-inline-rule-shadowing #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}

module Math.Linear.Matrix
    ( module Math.Linear.Matrix
    , ElemWise (..)
    ) where

import Foundation
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Class.Storable
import Foundation.Collection
import Foundation.Primitive

import GHC.TypeLits

import Prelude (fromIntegral)
import GHC.Num (Num)
import Foreign.Marshal.Alloc (alloca)
import qualified Foreign.Storable as Foreign.Storable
import System.IO.Unsafe (unsafePerformIO)

import Math.Linear.ElemWise
import qualified Math.Linear.Internal as I
import Math.Linear.Vector
import qualified Math.Linear.Matrix.Mutable as Mutable

newtype Mat a (m :: Nat) (n :: Nat) = M
    { vect :: Vec a (m * n) -- ^ data in plain vector.
    } deriving (Eq, Show)

type instance Element (Mat a m n) = a

-- | Construct an empty matrix.
empty :: PrimType a => Mat a 0 0
empty = M mempty

{-# INLINE empty #-}

-- | Verify matrix dimensions and memory layout
valid :: forall a m n. (KnownNat m, KnownNat n)
    => PrimType a => Mat a m n -> Bool
valid M{..} = row >= 0 && column >= 0 && length vect == integralCast (row * column)
    where row = natVal (Proxy :: Proxy m)
          column = natVal (Proxy :: Proxy n)

{-# INLINE valid #-}

-- | If the matrix is a square matrix.
square :: forall a m n. (KnownNat m, KnownNat n)
    => Mat a m n -> Bool
square M{..} = isJust (sameNat (Proxy :: Proxy m) (Proxy :: Proxy n))

{-# INLINE square #-}

-- | Construct a matrix with all zeros.
zeros :: (PrimType a, Num a, KnownNat m, KnownNat n) => proxy m -> proxy n -> Mat a m n
zeros r c = replicate' r c 0

{-# INLINE zeros #-}

-- | Construct a matrix with all ones.
ones :: (PrimType a, Num a, KnownNat m, KnownNat n) => proxy m -> proxy n -> Mat a m n
ones r c = replicate' r c 1

{-# INLINE ones #-}

-- | Construct a identity matrix, square is not required.
identity :: (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n))
    => proxy m -> proxy n -> Mat a m n
identity r c = unsafePerformIO $ do
    m <- Mutable.zeros r c
    Mutable.unsafeWith m $ \xs r' c' ->
        I.call $ I.identity xs r' c'
    unsafeFreeze m

{-# INLINE identity #-}

-- | Construct a random matrix.
random :: (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n))
    => proxy m -> proxy n -> IO (Mat a m n)
random r c = do
    m <- Mutable.new r c
    Mutable.unsafeWith m $ \xs r' c' ->
        I.call' $ I.random_ xs r' c'
    unsafeFreeze m

{-# INLINE random #-}

-- -- | Permutation matrix, assume that /P/ is the permutation matrix, then,
-- --
-- -- * exchange two columns: /x' = x A/.
-- --
-- -- * exchange two rows: /x' = A x/.
-- permute :: I.Elem a
--     => Int32 -- ^ size of matrix, /n x n/.
--     -> Int32
--     -> Int32
--     -> Mat a
-- permute n i j = unsafePerformIO $ do
--     m <- Mutable.identity n n
--     Mutable.unsafeWrite m i i 0
--     Mutable.unsafeWrite m j j 0
--     Mutable.unsafeWrite m i j 1
--     Mutable.unsafeWrite m j i 1
--     unsafeFreeze m

-- {-# INLINE permute #-}

-- | Construct a matrix with given constant.
replicate' :: (PrimType a, KnownNat m, KnownNat n)
    => proxy m -> proxy n -> a -> Mat a m n
replicate' m n = M . replicate (integralCast (r * c))
    where r = natVal m
          c = natVal n

{-# INLINE replicate' #-}

-- | Construct matrix with given generate function.
matrix :: (PrimType a, KnownNat m, KnownNat n)
    => proxy m -> proxy n -> (Int32 -> Int32 -> a) -> Mat a m n
matrix m n func = fromList' r c [ func i j
                                | j <- [0 .. c - 1]
                                , i <- [0 .. r - 1] ]
    where r = integralDownsize (natVal m)
          c = integralDownsize (natVal n)

{-# INLINE matrix #-}

-- | Construct a diagonal matrix from given vector, elements in vector will be diagonal elements of the matrix.
diag :: (I.Elem a, KnownNat n, KnownNat (n * n))
    => Proxy n -> UArray a -> Mat a n n
diag nlen xs = unsafePerformIO $ do
    m <- Mutable.new nlen nlen
    Mutable.unsafeWith m $ \vect r c ->
        withPtr xs $ \p ->
            I.call $ I.diag vect r c p
    unsafeFreeze m

{-# INLINE [1] diag #-}

-- | Construct matrix from given list and size, assume that /row * column == length xs/.
fromList' :: PrimType a
    => Int32 -- ^ rows
    -> Int32 -- ^ columns
    -> [a] -- ^ values
    -> Mat a m n
fromList' r c = M . fromListN (integralUpsize $ r * c)

{-# INLINE [1] fromList' #-}

-- {-# RULES
-- "fromlist'/tolist" forall a b m. fromList' a b (toList m) = m
--  #-}

-- # RULES
-- "tolist/fromlist'" forall a b l. toList (fromList' a b l) = l
--  # -- assume /a * b == length l/

-- | If the vector contains a specific element that statisfy some predicate.
find' :: PrimType a => (a -> Bool) -> Mat a m n -> Maybe a
find' predicate M{..} = find predicate vect

{-# INLINE find' #-}

-- -- | Left fold strictly, /O(n)/.
-- foldl'
--     :: Storable b
--     => (a -> b -> a) -> a -> Mat b -> a
-- foldl' f x = V.foldl' f x . vect

-- {-# INLINE foldl' #-}

-- | Get specified element from matrix.
at :: forall proxy a m n u v. (PrimType a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, u <= m, v <= n)
    => Mat a m n -> (proxy u, proxy v) -> Maybe a
at M{..} (u, v) = vect ! (integralCast (i * column + j))
    where i = natVal u
          j = natVal v
          column = natVal (Proxy :: Proxy m)

{-# INLINE at #-}

-- | Unsafely get specified element from matrix.
at' :: forall a m n u v. (PrimType a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, u <= m, v <= n)
    => Mat a m n -> (Proxy u, Proxy v) -> a
at' M{..} (u, v) = case vect ! (integralCast (i * column + j)) of
                     Just v' -> v'
                     Nothing -> error "Matrix.at': no such element."
    where i = natVal u
          j = natVal v
          column = natVal (Proxy :: Proxy m)

{-# INLINE at' #-}

-- | Get elements in diagonal positions, square matrix is not necessary.
diagonal :: forall a m n. (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n))
    => Mat a m n -> UArray a
diagonal m@M{..} = unsafePerformIO $ do
    v <- mutNew (integralCast (min row column))
    unsafeWith m $ \xs r c ->
        withMutablePtr v $ \p ->
            I.call $ I.diagonal p xs r c
    unsafeFreeze v
  where row = natVal (Proxy :: Proxy m)
        column = natVal (Proxy :: Proxy m)

{-# INLINE [1] diagonal #-}

-- # RULES
-- "diagonal/diag" forall v. diagonal (diag v) = v
--  #

-- | Get deminsion of the matrix.
shape :: forall a m n. (KnownNat m, KnownNat n)
    => Mat a m n -> (Int32, Int32)
shape M{..} = (integralDownsize (natVal (Proxy :: Proxy m)), integralDownsize (natVal (Proxy :: Proxy n)))

{-# INLINE [1] shape #-}

-- | Reshape a matrix.
reshape :: forall a m0 n0 m1 n1. (KnownNat m0, KnownNat n0, KnownNat m1, KnownNat n1, (m0 * n0) ~ (m1 * n1))
    => Mat a m0 n0 -> (Proxy m1, Proxy n1) -> Mat a m1 n1
reshape M{..} _ = M vect

{-# INLINE [1] reshape #-}

-- | Reshape a matrix.
reshape' :: forall a m0 n0 m1 n1. (KnownNat m0, KnownNat n0, KnownNat m1, KnownNat n1, (m0 * n0) ~ (m1 * n1))
    => Mat a m0 n0 -> Mat a m1 n1
reshape' M{..} = M vect

{-# INLINE reshape' #-}

-- # RULES
-- "reshape/reshape" forall a b m. reshape (reshape m a) b = reshape m b
-- "reshape/shape" forall m. reshape m (shape m) = m
--  #

-- | Infix operator of reshape.
(==>) :: forall a m0 n0 m1 n1. (KnownNat m0, KnownNat n0, KnownNat m1, KnownNat n1, (m0 * n0) ~ (m1 * n1))
    => Mat a m0 n0 -> (Proxy m1, Proxy n1) -> Mat a m1 n1
(==>) = reshape

infix 8 ==>

{-# INLINE (==>) #-}

-- -- | Lines of a matrix.
-- rows :: Storable a => Mat a -> [V.Vector a]
-- rows m@M {..} = Prelude.map (rowAt m) [0 .. row - 1]

-- {-# INLINE rows #-}

-- -- | Get specific row.
-- rowAt
--     :: Storable a
--     => Mat a -> Int -> V.Vector a
-- rowAt M {..} r = V.slice (r * column) column vect

-- {-# INLINE rowAt #-}

-- -- | Columns of a matrix.
-- columns
--     :: Storable a
--     => Mat a -> [V.Vector a]
-- columns m@M {..} = Prelude.map (columnAt m) [0 .. column - 1]

-- {-# INLINE columns #-}

-- -- | Get specific column.
-- columnAt
--     :: Storable a
--     => Mat a -> Int -> V.Vector a
-- columnAt M {..} c = V.generate row $ \i -> vect V.! (i * column + c)

-- {-# INLINE columnAt #-}

-- -- | Map an unary function over a matrix.
-- map
--     :: (Storable a, Storable b)
--     => (a -> b) -> Mat a -> Mat b
-- map f M {..} = M row column $ V.map f vect

-- {-# INLINE [1] map #-}

-- {-# RULES
-- "map/map" forall f g m. map f (map g m) = map (f . g) m
--  #-}

-- -- | Filter elements from a matrix ans accumulate these elements into a vector.
-- filter
--     :: Storable a
--     => (a -> Bool) -> Mat a -> V.Vector a
-- filter f = V.filter f . vect

-- {-# INLINE [1] filter #-}

-- {-# RULES
-- "filter/map" forall f g m. filter f (map g m) = filter (f . g) m
--  #-}

-- | Matrix transpose.
transpose :: forall a m n. (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n), KnownNat (n * m))
    => Mat a m n -> Mat a n m
transpose m@M{..} = unsafeUnaryOp (Proxy :: Proxy n) (Proxy :: Proxy m) I.transpose m -- matrix column row $ \i j -> m ! (j, i)

{-# INLINE [1] transpose #-}

-- # RULES
-- "transpose/transpose" forall m. transpose (transpose m) = m
--  #

-- | Lower triangularize the given matrix strictly.
lower :: forall a m n. (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n))
    => Mat a m n -> Mat a m n
lower m@M{..} = unsafeUnaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.lower m

{-# INLINE lower #-}

-- | Upper triangularize the given matrix strictly.
upper :: forall a m n. (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n))
    => Mat a m n -> Mat a m n
upper m@M{..} = unsafeUnaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.upper m

{-# INLINE upper #-}

-- -- | Concatenate mutiple matrices into one with row as axis, columns of these operands must be equal.
-- --
-- -- > concatenate [[1,2],[3,4]] [[5,6],[7,8]]
-- -- >    1 2
-- -- >    3 4
-- -- >    5 6
-- -- >    7 8
-- concatenate
--     :: Storable a
--     => [Mat a] -> Mat a
-- concatenate [] = error "Matrix.concatenate: no operand found."
-- concatenate ms
--     | length csc <= 1 = M rs (List.head csc) vects
--     | otherwise = error "Matrix.concatenate: all matrice must have same column size."
--   where
--     rs = Prelude.sum . Prelude.map row $ ms
--     csc = List.concat . List.group . Prelude.map column $ ms
--     vects = V.concat . Prelude.map vect $ ms

-- {-# INLINE concatenate #-}

-- -- | Horizontally join two matrices.
-- --
-- -- > ( A ) <|> ( B ) = ( A | B )
-- --
-- -- Where both matrices /A/ and /B/ have the same number of rows, but no bounds check will be performed.
-- (<|>)
--     :: Storable a
--     => Mat a -> Mat a -> Mat a
-- m1 <|> m2 =
--     matrix (row m1) (c + column m2) $
--     \i j ->
--          if j <= c
--              then m1 ! (i, j)
--              else m2 ! (i, j - c)
--   where
--     c = column m1

-- -- {-# INLINE (<|>) #-}

-- -- | Vertically join two matrices.
-- --
-- -- >                   ( A )
-- -- > ( A ) <-> ( B ) = ( - )
-- -- >                   ( B )
-- --
-- -- Where both matrices /A/ and /B/ have the same number of columns, but no bounds check will be performed.
-- (<->)
--     :: Storable a
--     => Mat a -> Mat a -> Mat a
-- m1 <-> m2 =
--     matrix (r + row m2) (column m1) $
--     \i j ->
--          if i <= r
--              then m1 ! (i, j)
--              else m2 ! (i - r, j)
--   where
--     r = row m1

-- {-# INLINE (<->) #-}

-- -- | Zip two matrix with given function and truncate extra elements both row and column.
-- zipWith
--     :: (Storable a, Storable b, Storable c)
--     => (a -> b -> c) -> Mat a -> Mat b -> Mat c
-- zipWith f m1@(M r1 c1 _) m2@(M r2 c2 _)
--     | r1 == r1 && c1 == c2 = M r1 c1 $ V.zipWith f (vect m1) (vect m2)
--     | r1 <= r2 && c1 <= c2 = matrix r1 c1 $ \i j -> f (m1 ! (i, j)) (m2 ! (i, j))
--     | r1 >= r2 && c1 >= c2 = matrix r2 c2 $ \i j -> f (m1 ! (i, j)) (m2 ! (i, j))
--     | otherwise = error "Mat.zipWith: uncompatible shape of two operands."

-- | The sum of all coefficients of the matrix
sum :: (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n)) => Mat a m n -> a
sum = unsafeOp I.sum

{-# INLINE sum #-}

-- | The product of all coefficients of the matrix
product :: (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n)) => Mat a m n -> a
product = unsafeOp I.product -- V.product . vect

{-# INLINE product #-}

-- | The mean of all coefficients of the matrix
mean :: (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n)) => Mat a m n -> a
mean = unsafeOp I.mean

{-# INLINE mean #-}

-- {-# INLINE [1] shift #-}

-- # RULES
-- "shift/shift" forall a b m. shift a (shift b m) = shift (a + b) m
--  #

-- {-# INLINE [1] times #-}

-- # RULES
-- "times/times" forall a b m. times a (times b m) = times (a * b) m
--  #

instance (I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n)) => ElemWise (Mat a m n) where
    -- * tensor and scalar
    shift x m = unsafePerformIO $ do
        m' <- Mutable.new (Proxy :: Proxy m) (Proxy :: Proxy n)
        Mutable.unsafeWith m' $ \v' _ _ ->
            unsafeWith m $ \v row column ->
                alloca $ \p -> do
                    poke p x
                    I.call $ I.shift v' p v row column
        unsafeFreeze m'
    times x m = unsafePerformIO $ do
        m' <- Mutable.new (Proxy :: Proxy m) (Proxy :: Proxy n)
        Mutable.unsafeWith m' $ \v' _ _ ->
            unsafeWith m $ \v row column ->
                alloca $ \p -> do
                    poke p x
                    I.call $ I.times v' p v row column
        unsafeFreeze m'
    -- * negative
    negative = unsafeUnaryOp (Proxy @m) (Proxy @n) I.negative
    -- * arithmetic
    add = unsafeBinaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.add
    minus = unsafeBinaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.minus
    mult = unsafeBinaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.mult
    division = unsafeBinaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.division
    -- * data generation
    constreplic x = unsafePerformIO $ do
        m' <- Mutable.new (Proxy @m) (Proxy @n)
        Mutable.unsafeWith m' $ \pm' row column ->
            alloca $ \p -> do
                poke p x
                I.call $ I.replicate pm' p (row * column)
        unsafeFreeze m'
    -- * extensions
    logistic = unsafeUnaryOp (Proxy @m) (Proxy @n) I.logistic
    logisticd = unsafeUnaryOp (Proxy @m) (Proxy @n) I.logisticd

-- | Matrix multiplication WITHOUT bounds check.
dot :: forall a m k n. (I.Elem a, KnownNat m, KnownNat k, KnownNat n, KnownNat (m * k), KnownNat (k * n), KnownNat (m * n))
    => Mat a m k -> Mat a k n -> Mat a m n
dot m1 m2 = unsafeBinaryOp (Proxy :: Proxy m) (Proxy :: Proxy n) I.dot m1 m2

{-# INLINE dot #-}

-- | Infix operator of dot.
(.*) :: (I.Elem a, KnownNat m, KnownNat k, KnownNat n, KnownNat (m * k), KnownNat (k * n), KnownNat (m * n))
    => Mat a m k -> Mat a k n -> Mat a m n
(.*) = dot

infix 8 .*

{-# INLINE (.*) #-}

-- -- | Matrix multiplication in pure Haskell without FFI.
-- dot'' :: I.Elem a => Mat a -> Mat a -> Mat a
-- dot'' m1 m2 = matrix (row m1) (column m2) $ \i j -> V.sum $ V.zipWith (*) (rowAt m1 i) (columnAt m2 j)

-- {-# INLINE dot'' #-}

-- | Matrix power with matrix multiplication (dot) as multiply operator, the given matrix must
-- be a square matrix.
pow :: forall a b n. (I.Elem a, Integral b, Ord b, Num b, IDivisible b, IsNatural b, KnownNat n, KnownNat (n * n))
    => Mat a n n -> b -> Mat a n n
pow m@M {..} k
    | otherwise = go (diag (Proxy :: Proxy n) $ replicate (integralCast nlen) 1) m k
  where
    nlen = natVal (Proxy :: Proxy n)

    go r _ 0 = r
    go r a n = go r' a' n'
      where
        r' = if n `mod` 2 == 0
                then r
                else r `dot` a
        a' = a `dot` a
        n' = n `div` 2

{-# INLINE pow #-}

-- | inner product of two specific columns (column i and column j) WITHOUT bounds check.
--
--    * For real vector, r = x^T y
--    * For complex vector, r = x^H y
inner :: (I.Elem a, KnownNat m, KnownNat n, KnownNat u, KnownNat v, KnownNat (m * n), u <= m, v <= n)
    => Mat a m n -> Proxy u -> Proxy v -> a
inner m u v = unsafePerformIO $
    unsafeWith m $ \xs r c ->
        alloca $ \p -> do
            I.call $ I.inner p r xs c i xs c j
            peek p
  where
    i = integralDownsize $ natVal u
    j = integralDownsize $ natVal v

{-# INLINE inner #-}

-- | inner product of two vectors with same length, but no bounds check will be performed.
inner' :: I.Elem a
    => Vec a n -> Vec a n -> a
inner' (V v1) (V v2) = unsafePerformIO $
    withPtr v1 $ \xs1 ->
        withPtr v2 $ \xs2 ->
            alloca $ \p -> do
                I.call $ I.inner p n xs1 1 0 xs2 1 0
                peek p
  where
    n = integralDownsize $ min (length v1) (length v2)

{-# INLINE inner' #-}

-- slice :: I.Elem a
    -- => Mat a m n -> Proxy u -> Proxy v -> Mat a m n

-- -- | Copy the matrix data and drop extra memory.
-- force :: forall a m n. (PrimType a, KnownNat m, KnownNat n)
--     => Mat a m n -> Mat a m n
-- force m@M{..} = matrix' row column $ curry (m `at'`)
--     where row = integralCast $ natVal (Proxy :: Proxy m)
--           column = integralCast $ natVal (Proxy :: Proxy n)

-- {-# INLINE force #-}

-- | Pass a pointer to the matrix's data to the IO action. The data may not be modified through the pointer.
unsafeWith :: forall monad a c m n. (PrimMonad monad, I.Elem a, KnownNat m, KnownNat n, KnownNat (m * n))
    => Mat a m n -> (Ptr a -> Int32 -> Int32 -> monad c) -> monad c
unsafeWith M{..} f = withVPtr vect $ \p -> f p row column
    where row = integralDownsize $ natVal (Proxy :: Proxy m)
          column = integralDownsize $ natVal (Proxy :: Proxy n)

{-# INLINE unsafeWith #-}

withMPtr :: (PrimMonad monad, PrimType a) => Mat a m n -> (Ptr a -> monad b) -> monad b
withMPtr M{..} = withVPtr vect

{-# INLINE withMPtr #-}

-- | Take one matrix as argument and only return a scalar, like sum, mean, max and min.
unsafeOp :: (I.Elem a, Foreign.Storable.Storable b, Storable b, KnownNat m0, KnownNat n0, KnownNat (m0 * n0))
    => (Ptr b -- result ptr, a scalar value.
        -> Ptr a -- matrix operand
        -> Int32 -- rows of matrix operand
        -> Int32 -- columns of matrix operand
        -> Int32 -- CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a m0 n0 -- ^ operand
    -> b
unsafeOp f m = unsafePerformIO $
    unsafeWith m $ \xs r c ->
        alloca $ \p -> do
            I.call $ f p xs r c
            peek p

{-# INLINE unsafeOp #-}

-- | Take one matrix as it's argument and return a matrix as result, like transpose and inverse.
unsafeUnaryOp :: (I.Elem a, KnownNat m0, KnownNat n0, KnownNat m1, KnownNat n1, KnownNat (m0 * n0), KnownNat (m1 * n1))
    => Proxy m1 -- ^ target rows
    -> Proxy n1 -- ^ target columns
    -> (Ptr a -- ^ result ptr, a matrix
        -> Ptr a -- ^ matrix operand ptr
        -> Int32 -- ^ rows of matrix operand
        -> Int32 -- ^ columns of matrix operand
        -> Int32 -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a m0 n0 -- ^ matrix operand
    -> Mat a m1 n1
unsafeUnaryOp r c f m = unsafePerformIO $ do
    m' <- Mutable.new r c
    Mutable.unsafeWith m' $ \xs' _ _ ->
        unsafeWith m $ \xs u v ->
            I.call $ f xs' xs u v
    unsafeFreeze m'

{-# INLINE unsafeUnaryOp #-}

-- | Take two matrices as arguments and return a matrix as result, like dot, add, multiplication.
unsafeBinaryOp :: (I.Elem a, KnownNat m0, KnownNat n0, KnownNat ml, KnownNat nl, KnownNat mr, KnownNat nr, KnownNat (m0 * n0), KnownNat (ml * nl), KnownNat (mr * nr))
    => Proxy m0 -- ^ target rows
    -> Proxy n0 -- ^ target columns
    -> (Ptr a -- ^ result ptr, a matrix, with shape /m x n/
        -> Int32 -- ^ m
        -> Int32 -- ^ n
        -> Int32 -- ^ k
        -> Ptr a -- ^ left matrix operand, with shape /m x k/
        -> Ptr a -- ^ right matrix operand, with shape /k x n/
        -> Int32 -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a ml nl -- ^ left operand
    -> Mat a mr nr -- ^ right operand
    -> Mat a m0 n0
unsafeBinaryOp r c f m1 m2 = unsafePerformIO $ do
  m0 <- Mutable.new r c
  Mutable.unsafeWith m0 $
    \vect0 m n ->
      unsafeWith m1 $
        \vect1 _ k -> unsafeWith m2 $ \vect2 _ _ -> I.call $ f vect0 m n k vect1 vect2
  unsafeFreeze m0

{-# INLINE unsafeBinaryOp #-}

-- | Take one matrix as argument and two matrices as result, like QR factorization.
unsafeFactorizeOp
    :: (I.Elem a, I.Elem b, I.Elem c, KnownNat m0, KnownNat n0, KnownNat m1, KnownNat n1, KnownNat m2, KnownNat n2, KnownNat (m0 * n0), KnownNat (m1 * n1), KnownNat (m2 * n2))
    => Mat a m0 n0 -- ^ matrix to decompose
    -> (Proxy m1, Proxy n1) -- ^ size of result 1
    -> (Proxy m2, Proxy n2) -- ^ size of result 2
    -> (Ptr a -- ^ matrix operand
        -> Int32 -> Int32 -- ^ size of operand
        -> Ptr b -- ^ result ptr 1, a matrix
        -> Int32 -> Int32 -- ^ size of result 1
        -> Ptr c -- ^ result ptr 2, a matrix
        -> Int32 -> Int32 -- ^ size of result 2
        -> Int32 -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> (Mat b m1 n1, Mat c m2 n2)
unsafeFactorizeOp m0 (r1, c1) (r2, c2) f = unsafePerformIO $ do
    m1 <- Mutable.new r1 c1
    m2 <- Mutable.new r2 c2
    Mutable.unsafeWith m1 $ \vect1 r1' c1' ->
        Mutable.unsafeWith m2 $ \vect2 r2' c2' ->
            unsafeWith m0 $ \vect0 r0 c0 -> I.call $ f vect0 r0 c0 vect1 r1' c1' vect2 r2' c2'
    m1' <- unsafeFreeze m1
    m2' <- unsafeFreeze m2
    return (m1', m2')

{-# INLINE unsafeFactorizeOp #-}

-- | Take one matrix as argument and three matrices as result, like LU decomposition and SVD decomposition.
unsafeDecomposeOp
    :: (I.Elem a, I.Elem b, I.Elem c, I.Elem d, KnownNat m0, KnownNat n0, KnownNat m1, KnownNat n1, KnownNat m2, KnownNat n2, KnownNat m3, KnownNat n3, KnownNat (m0 * n0), KnownNat (m1 * n1), KnownNat (m2 * n2), KnownNat (m3 * n3))
    => Mat a m0 n0 -- ^ matrix to decompose
    -> (Proxy m1, Proxy n1) -- ^ size of result 1
    -> (Proxy m2, Proxy n2) -- ^ size of result 2
    -> (Proxy m3, Proxy n3) -- ^ size of result 3
    -> (Ptr a -- ^ matrix operand
        -> Int32 -> Int32 -- ^ size of operand
        -> Ptr b -- ^ result ptr 1, a matrix
        -> Int32 -> Int32 -- ^ size of result 1
        -> Ptr c -- ^ result ptr 2, a matrix
        -> Int32 -> Int32 -- ^ size of result 2
        -> Ptr d -- ^ result ptr 3, a matrix
        -> Int32 -> Int32 -- ^ size of result 3
        -> Int32 -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> (Mat b m1 n1, Mat c m2 n2, Mat d m3 n3)
unsafeDecomposeOp m0 (r1, c1) (r2, c2) (r3, c3) f = unsafePerformIO $ do
    m1 <- Mutable.new r1 c1
    m2 <- Mutable.new r2 c2
    m3 <- Mutable.new r3 c3
    Mutable.unsafeWith m1 $ \vect1 r1' c1' ->
        Mutable.unsafeWith m2 $ \vect2 r2' c2' ->
            Mutable.unsafeWith m3 $ \vect3 r3' c3' ->
                unsafeWith m0 $ \vect0 r0 c0 -> I.call $ f vect0 r0 c0 vect1 r1' c1' vect2 r2' c2' vect3 r3' c3'
    m1' <- unsafeFreeze m1
    m2' <- unsafeFreeze m2
    m3' <- unsafeFreeze m3
    return (m1', m2', m3')

{-# INLINE unsafeDecomposeOp #-}

instance (PrimType a, KnownNat m, KnownNat n) => IsList (Mat a m n) where
    type Item (Mat a m n) = a
    fromList xs
        | m' * n' == nlen = M (fromList xs)
        | otherwise = error "Matrix.fromList: the give size doesn't match"
        where m' = natVal (Proxy :: Proxy m)
              n' = natVal (Proxy :: Proxy n)
              nlen = let CountOf x = length xs in integralUpsize x
    fromListN n xs
        | m' * n' == nlen = M (fromListN n xs)
        | otherwise = error "Matrix.fromListN: the give size doesn't match"
        where m' = natVal (Proxy :: Proxy m)
              n' = natVal (Proxy :: Proxy n)
              nlen = let CountOf x = length xs in integralUpsize x
    toList M{..} = toList vect

deriving instance NormalForm (Mat a m n)

deriving instance PrimType a => Fold1able (Mat a m n)

deriving instance PrimType a => Foldable (Mat a m n)

deriving instance PrimType a => IndexedCollection (Mat a m n)

deriving instance PrimType a => InnerFunctor (Mat a m n)

deriving instance PrimType a => Copy (Mat a m n)

deriving instance (PrimType a, KnownNat m, KnownNat n) => Collection (Mat a m n)

instance (PrimType a, KnownNat m, KnownNat n, KnownNat (m * n)) => MutableCollection (Mutable.MMat a m n) where
    type MutableFreezed (Mutable.MMat a m n) = Mat a m n
    type MutableKey (Mutable.MMat a m n) = (Offset a, Offset a)
    type MutableValue (Mutable.MMat a m n) = a
    unsafeThaw (M v) = Mutable.MM <$> unsafeThaw v
    unsafeFreeze (Mutable.MM v) = M <$> unsafeFreeze v
    thaw (M v) = Mutable.MM <$> thaw v
    freeze (Mutable.MM v) = M <$> freeze v
    -- The given size argument would be ignored.
    mutNew :: forall monad. PrimMonad monad => CountOf a -> monad (Mutable.MMat a m n (PrimState monad))
    mutNew ~_ = Mutable.MM <$> mutNew undefined :: monad (Mutable.MMat a m n (PrimState monad))
    mutUnsafeWrite (Mutable.MM v) (Offset i, Offset j) =
        mutUnsafeWrite v (Offset (i * c + j))
            where c = fromIntegral $ natVal (Proxy :: Proxy n)
    mutWrite (Mutable.MM v) (Offset i, Offset j) =
        mutWrite v (Offset (i * c + j))
            where c = fromIntegral $ natVal (Proxy :: Proxy n)
    mutUnsafeRead (Mutable.MM v) (Offset i, Offset j) =
        mutUnsafeRead v (Offset (i * c + j))
            where c = fromIntegral $ natVal (Proxy :: Proxy n)
    mutRead (Mutable.MM v) (Offset i, Offset j) =
        mutRead v (Offset (i * c + j))
            where c = fromIntegral $ natVal (Proxy :: Proxy n)
