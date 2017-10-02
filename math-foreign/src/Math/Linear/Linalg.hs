-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Linalg
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Linear algebra.
--
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}

module Math.Linear.Linalg
  ( Mat(..)
   -- * Matrix Properties
  , det
  , trace
  , rank
  , norm
   -- * Matrix inverse
  , inverse
   -- * Eigen system
  , eigen
  , eigenh
  , eigenvals
  , eigenvalsh
   -- * Matrix decomposition
  , lu
  , qr
  , qr'
  , svd
  , svd'
  , singular
  , jordan
  , cholesky
  , schur
   -- * Linear transformations
  , transform)
  where

import Foundation
import Foundation.Collection
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Foreign
import Foundation.Primitive

import Foreign.C.String (castCharToCChar)
import Foreign.Ptr (nullPtr)
import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import qualified Math.Linear.Matrix.Mutable as Mutable
import Math.Linear.Matrix (Mat(..), unsafeWith, unsafeOp, unsafeUnaryOp, unsafeFactorizeOp, unsafeDecomposeOp)

-- | Matrix determinant.
det :: I.Elem a
    => Mat a -> a
det = unsafeOp I.det

{-# INLINE det #-}

-- | Trace of a matrix.
trace :: I.Elem a
    => Mat a -> a
trace = unsafeOp I.trace

{-# INLINE trace #-}

-- | Rank of a matrix.
rank :: I.Elem a
    => Mat a -> Int32
rank = unsafeOp I.rank

{-# INLINE rank #-}

-- | .
norm :: I.Elem a
    => Mat a -> a
norm = unsafeOp I.norm

{-# INLINE norm #-}

-- | Inverse matrix.
inverse :: I.Elem a
    => Mat a -> Mat a
inverse m@M{..} = unsafeUnaryOp row column I.inverse m

{-# INLINE inverse #-}

-- | Compute the eigenvalues and right eigenvectors of a square matrix.
eigen :: I.Elem a
    => Mat a -> Mat a
eigen m@M{..} = undefined

{-# INLINE eigen #-}

-- | Compute the eigenvalues and right eigenvectors of a Hermitian or symmetric matrix.
eigenh :: I.Elem a
    => Mat a -> Mat a
eigenh m@M{..} = undefined

{-# INLINE eigenh #-}

-- | Compute the eignevalues of a square matrix.
eigenvals :: I.Elem a
    => Mat a -> UArray a
eigenvals m@M{..} = undefined

{-# INLINE eigenvals #-}

-- | Compute the eignevalues of a Hermitian or symmetric matrix.
eigenvalsh :: I.Elem a
    => Mat a -> UArray a
eigenvalsh m@M{..} = undefined

{-# INLINE eigenvalsh #-}

-- | LU decomposition, for matrix A, return three matrices L, U and P as the result, where
-- L is a lower triangular matrix, U is a upper triangular matrix and P is the permutation matirx.
-- We have /AP = LU/.
lu :: I.Elem a
    => Mat a -- ^ matrix A
    -> (Mat a, Mat a, Mat a) -- ^ matrices (L, U, P)
lu m@M {..} = undefined
    -- | row < column = unsafeDecomposeOp (row, row) (row, column) (column, column) I.lu m
    -- | otherwise = unsafeDecomposeOp (row, column) (column, column) (column, column) I.lu m

{-# INLINE lu #-}

-- | Reduced QR decomposition: A = QR where A is a m-by-n matrix, Q is a m-by-m
-- orthogonal matrix (column vectors of Q are orthonormal), R is a m-by-n upper trapezoidal matrix.
--
-- Equals to numpy's numpy.linalg.qr(src, mode='reduced')
qr :: I.Elem a
    => Mat a -- ^ matrix A
    -> (Mat a, Mat a) -- ^ matrices (Q, R)
qr m@M{..} = unsafeFactorizeOp m (row, k) (k, column) I.qr
    where k = min row column

{-# INLINE qr #-}

-- | Complete QR decomposition: A = QR where A is a m-by-n matrix, Q is a m-by-min(m, n)
-- orthogonal matrix (column vectors of Q are orthonormal), R is a min(m, n)-by-n upper trapezoidal matrix.
--
-- Equals to numpy's numpy.linalg.qr(src, mode='complete')
qr' :: I.Elem a
    => Mat a -- ^ matrix A
    -> (Mat a, Mat a) -- ^ matrices (Q, R)
qr' m@M{..} = unsafeFactorizeOp m (row, row) (row, column) I.qr

{-# INLINE qr' #-}

-- | Singular value decomposition: A = U \Sigma V^T (or V^H)
svd :: (I.Elem a, I.Elem (I.RealType a))
    => Mat a -- ^ matrix A
    -> (Mat a, UArray (I.RealType a), Mat a) -- ^ matrices (U, \sigma, V^T (or V^H)), where \sigma is
                                             -- diagonal elements of matrix \Sigma.
svd m@M{..} = (u, sigma, vt)
    where k = min row column
          (u, (M _ _ sigma), vt) = unsafeDecomposeOp m (row, row) (k, 1) (column, column) I.svd

{-# INLINE svd #-}

-- | Compact singular value decomposition: A = U \Sigma V^T (or V^H). Only first min(m, n) columns of U and
-- first min(m, n) rows of V^H is computed, and \Sigma is a min(m, n)-by-min(m, n) square matrix.
svd' :: (I.Elem a, I.Elem (I.RealType a))
    => Mat a -- ^ matrix A
    -> (Mat a, UArray (I.RealType a), Mat a) -- ^ matrices (U, \sigma, V^T (or V^H)), where \sigma is
                                             -- diagonal elements of matrix \Sigma.
svd' m@M{..} = (u, sigma, vt)
    where k = min row column
          (u, (M _ _ sigma), vt) = unsafeDecomposeOp m (row, k) (k, 1) (k, column) I.svd

{-# INLINE svd' #-}

-- | Singular value decomposition, only result singular values, without compute U and V^T.
singular :: (I.Elem a, I.Elem (I.RealType a))
    => Mat a -- ^ matrix A
    -> UArray (I.RealType a) -- ^ singular value array.
singular M{..} = unsafePerformIO $ do
    arr <- mutNew (integralCast (min row column))
    withMutablePtr arr $ \pr ->
        withPtr vect $ \px ->
            I.call $ I.svd px row column nullPtr 0 1 pr 1 (min row column) nullPtr 0 2 -- TODO why the last argument `c3` can be set as 1 ?
    unsafeFreeze arr

{-# INLINE singular #-}

-- | Jordan decomposition.
jordan :: I.Elem a
    => Mat a -> (Mat a, Mat a, Mat a)
jordan m@M {..} = undefined

{-# INLINE jordan #-}

-- | Cholesky decomposition: A = U^T U for real data and A = U^H U for complex data, where U is a upper/lower triangular matrix.
-- A must be a symmetric positive-definite matrix (but not checked).
cholesky :: I.Elem a
    => Mat a
    -> Char   -- ^ 'U' for upper triangular matrix and 'L' for lower triangular matrix.
    -> Mat a  -- ^ an upper/lower triangular matrix.
cholesky m@M{..} uplo
    | row /= column = error "Cholesky decomposition can only be applied to a symmetric positive-definite matrix."
    | otherwise = unsafePerformIO $ do
        mu <- Mutable.new row column
        Mutable.unsafeWith mr $ \pr _ _ ->
            unsafeWith m $ \pm ->
                I.call $ I.cholesky (castCharToCChar uplo) pm row column pr row column
        unsafeFreeze mu

{-# INLINE cholesky #-}

-- | Schur decomposition.
schur :: I.Elem a
    => Mat a -> (Mat a, Mat a, Mat a)
schur m@M {..} = undefined

{-# INLINE schur #-}

-- | Linear transformation, /transform A x = Ax/.
transform :: I.Elem a
    => Mat a -> UArray a -> UArray a
transform m@M {..} v = unsafePerformIO $ do
    v' <- mutNew (integralCast row)
    withMutablePtr v' $ \xr ->
        unsafeWith m $ \xs r c ->
            withPtr v $ \xv ->
                I.call $ I.transform xr r c xs xv
    unsafeFreeze v'
