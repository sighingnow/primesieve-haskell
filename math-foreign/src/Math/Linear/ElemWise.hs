-- {-# LANGUAGE TypeFamilies #-}

module Math.Linear.ElemWise where

import Foundation
import Foundation.Collection

-- | Types that support elementwise operations.
class ElemWise tensor where
  -- * tensor and scalar
  shift :: Element tensor -> tensor -> tensor
  times :: Element tensor -> tensor -> tensor
  -- * negative
  negative :: tensor -> tensor
  -- * arithmetic
  add :: tensor -> tensor -> tensor
  minus :: tensor -> tensor -> tensor
  mult :: tensor -> tensor -> tensor
  division :: tensor -> tensor -> tensor
  -- * extensions
  logistic :: tensor -> tensor
  logisticd :: tensor -> tensor

class MutElemWise tensor where
  -- * tensor and scalar
  shift' :: Element tensor -> tensor -> IO ()
  times' :: Element tensor -> tensor -> IO ()
  -- * negative
  negative' :: tensor -> IO ()
  -- * arithmetic
  add' :: tensor -> tensor -> IO ()
  minus' :: tensor -> tensor -> IO ()
  mult' :: tensor -> tensor -> IO ()
  division' :: tensor -> tensor -> IO ()
  -- * extensions
  logistic' :: tensor -> IO ()
  logisticd' :: tensor -> IO ()
