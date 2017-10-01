-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Linalg
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Linear algebra.
--
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
  , svd
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

import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import qualified Math.Linear.Matrix as M
import Math.Linear.Matrix (Mat(..), unsafeDecomposeOp, unsafeOp, unsafeUnaryOp)

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
lu m@M {..}
    | row < column = unsafeDecomposeOp (row, row) (row, column) (column, column) I.lu m
    | otherwise = unsafeDecomposeOp (row, column) (column, column) (column, column) I.lu m

{-# INLINE lu #-}

-- | QR decomposition.
qr :: I.Elem a
    => Mat a -> (Mat a, Mat a, Mat a)
qr m@M {..} = undefined

{-# INLINE qr #-}

-- | Singular value decomposition.
svd :: I.Elem a
    => Mat a -> (Mat a, Mat a, Mat a)
svd m@M {..} = undefined

{-# INLINE svd #-}

-- | Jordan decomposition.
jordan :: I.Elem a
    => Mat a -> (Mat a, Mat a, Mat a)
jordan m@M {..} = undefined

{-# INLINE jordan #-}

-- | Cholesky decomposition.
cholesky :: I.Elem a
    => Mat a -> (Mat a, Mat a, Mat a)
cholesky m@M {..} = undefined

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
        M.unsafeWith m $ \xs r c ->
            withPtr v $ \xv ->
                I.call $ I.transform xr r c xs xv
    unsafeFreeze v'
