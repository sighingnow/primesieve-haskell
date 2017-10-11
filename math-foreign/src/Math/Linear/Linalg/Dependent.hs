-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Linalg.Dependent
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Linear algebra.
--
{-# OPTIONS_GHC -fprint-explicit-kinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Math.Linear.Linalg.Dependent where

import Foundation
import Foundation.Collection
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Foreign
import Foundation.Primitive

import GHC.TypeLits
import Foreign.C.String (castCharToCChar)
import Foreign.Ptr (nullPtr)
import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import qualified Math.Linear.Matrix.Mutable.Dependent as Mutable
import Math.Linear.Matrix.Dependent (Mat(..), unsafeWith, unsafeOp, unsafeUnaryOp, unsafeFactorizeOp, unsafeDecomposeOp)

-- | Matrix determinant.
det :: (I.Elem a, KnownNat m, KnownNat n)
    => Mat a m n -> a
det = unsafeOp I.det

{-# INLINE det #-}

-- | Trace of a matrix.
trace :: (I.Elem a, KnownNat m, KnownNat n)
    => Mat a m n -> a
trace = unsafeOp I.trace

{-# INLINE trace #-}

-- | Rank of a matrix.
rank :: (I.Elem a, KnownNat m, KnownNat n)
    => Mat a m n -> Int32
rank = unsafeOp I.rank

{-# INLINE rank #-}

-- | .
norm :: (I.Elem a, KnownNat m, KnownNat n)
    => Mat a m n -> a
norm = unsafeOp I.norm

{-# INLINE norm #-}

-- | Inverse matrix.
inverse :: forall a n. (I.Elem a, KnownNat n)
    => Mat a n n -> Mat a n n
inverse m@M{..} = unsafeUnaryOp nlen nlen I.inverse m
    where nlen = Proxy :: Proxy n

{-# INLINE inverse #-}

-- | Compute the eigenvalues and right eigenvectors of a square matrix.
--
--    + For every column vector x in right eigenvectors V, we have Ax = \lambda x.
--    + For every column vector x in left eigenvectors U, we have x^H * A = \lambda x^H.
--
-- Noticing that eigenvectors are stored as column vectors in U^T and V.
eigen :: forall a n. (I.Elem a, I.Elem (I.ComplexType a), KnownNat n)
    => Mat a n n -- ^ square matrix.
    -> (UArray (I.ComplexType a), Mat (I.ComplexType a) n n, Mat (I.ComplexType a) n n) -- ^ (\Lambda, UT (Hermitian of left eigenvectors), V (right eigenvectors))
eigen m@M{..} = (lam, ut, v)
    where nlen = Proxy :: Proxy n
          ((M lam), ut, v) = unsafeDecomposeOp m (Proxy :: Proxy 1, nlen) (nlen, nlen) (nlen, nlen) I.eigen

{-# INLINE eigen #-}

-- | Compute the eigenvalues and right eigenvectors of a Hermitian or symmetric matrix.
--
--  A = V \Lambda V^T (or V^H), V is an orthogonal matrix whose columns are the eigenvectors of A.
eigenh :: forall a n. (I.Elem a, I.Elem (I.RealType a), KnownNat n)
    => Mat a n n -- ^ symmetric matrix A, the upper triangular part of A is used when call lapack routines.
    -> (UArray (I.RealType a), Mat a n n) -- ^ The eigenvalues is stored in ascending order,
                                      -- eigenvalues of symmetric complex matrix is real number.
eigenh m@M{..} = (lam, v)
    where nlen = Proxy :: Proxy n
          ((M lam), v) = unsafeFactorizeOp m (Proxy :: Proxy 1, nlen) (nlen, nlen) I.eigenh

{-# INLINE eigenh #-}

-- | Compute the eignevalues of a square matrix.
eigenvals :: forall a n. (I.Elem a, I.Elem (I.ComplexType a), KnownNat n)
    => Mat a n n
    -> UArray (I.ComplexType a)
eigenvals m@M{..} = lam
    where nlen = Proxy :: Proxy n
          ((M lam), _, _) = unsafeDecomposeOp m (Proxy :: Proxy 1, nlen) (Proxy :: Proxy 0, nlen) (Proxy :: Proxy 0, nlen) I.eigen

{-# INLINE eigenvals #-}

-- | Compute the eignevalues of a Hermitian or symmetric matrix.
eigenvalsh :: forall a n. (I.Elem a, I.Elem (I.RealType a), KnownNat n)
    => Mat a n n
    -> UArray (I.RealType a) -- ^ The eigenvalues is stored in ascending order,
                             -- eigenvalues of symmetric complex matrix is real number.
eigenvalsh m@M{..} = lam
    where nlen = Proxy :: Proxy n
          ((M lam), _) = unsafeFactorizeOp m (Proxy :: Proxy 1, nlen) (Proxy :: Proxy 0, Proxy :: Proxy 1) I.eigenh

{-# INLINE eigenvalsh #-}

-- | LU decomposition, for matrix A, return three matrices L, U and P as the result, where
-- L is a lower triangular matrix, U is a upper triangular matrix and P is the permutation matirx.
-- We have /AP = LU/.
lu :: I.Elem a
    => Mat a m n -- ^ matrix A
    -> (Mat a m n, Mat a m n, Mat a m n) -- ^ matrices (L, U, P)
lu m@M{..} = undefined
    -- | row < column = unsafeDecomposeOp (row, row) (row, column) (column, column) I.lu m
    -- | otherwise = unsafeDecomposeOp (row, column) (column, column) (column, column) I.lu m

{-# INLINE lu #-}

-- | Reduced QR decomposition: A = QR where A is a m-by-n matrix, Q is a m-by-min(m, n)
-- orthogonal matrix (column vectors of Q are orthonormal), R is a min(m, n)-by-n upper trapezoidal matrix.
--
-- Equals to numpy's numpy.linalg.qr(src, mode='reduced')
-- qr :: I.Elem a
--     => Mat a m n -- ^ matrix A
--     -> (Mat a m m, Mat a m n) -- ^ matrices (Q, R)
-- qr m@M{..} = unsafeFactorizeOp m (row, k) (k, column) I.qr
--     where k = min row column

-- {-# INLINE qr #-}

-- | Complete QR decomposition: A = QR where A is a m-by-n matrix, Q is a m-by-m
-- orthogonal matrix (column vectors of Q are orthonormal), R is a m-by-n upper trapezoidal matrix.
--
-- Equals to numpy's numpy.linalg.qr(src, mode='complete')
qr' :: forall a m n. (I.Elem a, KnownNat m, KnownNat n)
    => Mat a m n -- ^ matrix A
    -> (Mat a m m, Mat a m n) -- ^ matrices (Q, R)
qr' m@M{..} = unsafeFactorizeOp m (row, row) (row, column) I.qr
    where row = Proxy :: Proxy m
          column = Proxy :: Proxy n

{-# INLINE qr' #-}

-- -- | Singular value decomposition: A = U \Sigma V^T (or V^H)
-- svd :: (I.Elem a, I.Elem (I.RealType a))
--     => Mat a m n -- ^ matrix A
--     -> (Mat a m m, UArray (I.RealType a), Mat a n n) -- ^ matrices (U, \sigma, V^T (or V^H)), where \sigma is
--                                              -- diagonal elements of matrix \Sigma.
-- svd m@M{..} = (u, sigma, vt)
--     where k = min row column
--           (u, (M _ _ sigma), vt) = unsafeDecomposeOp m (row, row) (k, 1) (column, column) I.svd

-- {-# INLINE svd #-}

-- -- | Compact singular value decomposition: A = U \Sigma V^T (or V^H). Only first min(m, n) columns of U and
-- -- first min(m, n) rows of V^H is computed, and \Sigma is a min(m, n)-by-min(m, n) square matrix.
-- svd' :: (I.Elem a, I.Elem (I.RealType a))
--     => Mat a -- ^ matrix A
--     -> (Mat a, UArray (I.RealType a), Mat a) -- ^ matrices (U, \sigma, V^T (or V^H)), where \sigma is
--                                              -- diagonal elements of matrix \Sigma.
-- svd' m@M{..} = (u, sigma, vt)
--     where k = min row column
--           (u, (M _ _ sigma), vt) = unsafeDecomposeOp m (row, k) (k, 1) (k, column) I.svd

-- {-# INLINE svd' #-}

-- | Singular value decomposition, only result singular values, without compute U and V^T.
singular :: forall a m n. (I.Elem a, I.Elem (I.RealType a), KnownNat m, KnownNat n)
    => Mat a m n -- ^ matrix A
    -> UArray (I.RealType a) -- ^ singular value array.
singular M{..} = unsafePerformIO $ do
    arr <- mutNew (integralCast (min row column))
    withMutablePtr arr $ \pr ->
        withPtr vect $ \px ->
            I.call $ I.svd px row column nullPtr 0 1 pr 1 (min row column) nullPtr 0 column
    unsafeFreeze arr
  where
    row = integralDownsize $ natVal (Proxy :: Proxy m)
    column = integralDownsize $ natVal (Proxy :: Proxy n)

{-# INLINE singular #-}

-- | Jordan decomposition.
jordan :: I.Elem a
    => Mat a m n -> (Mat a m n, Mat a m n, Mat a m n)
jordan m@M{..} = undefined

{-# INLINE jordan #-}

-- | Cholesky decomposition: A = U^T U for real data and A = U^H U for complex data, where U is a upper/lower triangular matrix.
-- A must be a symmetric positive-definite matrix (but not checked).
cholesky :: forall a n. (I.Elem a, KnownNat n)
    => Mat a n n
    -> Char       -- ^ 'U' for upper triangular matrix and 'L' for lower triangular matrix.
    -> Mat a n n  -- ^ an upper/lower triangular matrix.
cholesky m@M{..} uplo = unsafePerformIO $ do
    mu <- Mutable.new (Proxy :: Proxy n) (Proxy :: Proxy n)
    Mutable.unsafeWith mu $ \pr r1 c1 ->
        unsafeWith m $ \pm r0 c0 ->
            I.call $ I.cholesky (castCharToCChar uplo) pm r0 c0 pr r1 c1
    unsafeFreeze mu

{-# INLINE cholesky #-}

-- | Schur decomposition.
schur :: I.Elem a
    => Mat a m n -> (Mat a m n, Mat a m n, Mat a m n)
schur m@M {..} = undefined

{-# INLINE schur #-}

-- | Linear transformation, /transform A x = Ax/.
transform :: forall a m n. (I.Elem a, KnownNat m, KnownNat n)
    => Mat a m n -> UArray a -> UArray a
transform m@M{..} v
    | row * column /= integralUpsize nlen' = error "Linalg.transform: the given size doesn't match"
    | otherwise = unsafePerformIO $ do
        v' <- mutNew nlen
        withMutablePtr v' $ \xr ->
            unsafeWith m $ \xs r c ->
                withPtr v $ \xv ->
                    I.call $ I.transform xr r c xs xv
        unsafeFreeze v'
  where row = natVal (Proxy :: Proxy m)
        column = natVal (Proxy :: Proxy n)
        nlen = length v
        nlen' = let CountOf x = nlen in x
