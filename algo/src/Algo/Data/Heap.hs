-------------------------------------------------------------------------------
-- |
-- module:          Algo.Data.Heap
--
-- Heap type constraints.

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}

module Algo.Data.Heap where

import Foundation
import GHC.TypeLits

class Heap (heap :: * -> Nat -> *) where
    -- * Properties
    nil :: Ord a => heap a 0

    toHeap :: Ord a => a -> heap a 1

    size :: heap a s -> Int

    isEmpty :: heap a s -> Bool

    minElem :: heap a s -> a

    -- * Operations

    deleteMin :: Ord a => heap a s -> heap a (s - 1)

    merge :: Ord a => heap a m -> heap a n -> heap a (m + n)

    -- * Default operations

    insert :: Ord a => a -> heap a s -> heap a (1 + s)
    insert x = merge (toHeap x)
    {-# INLINE insert #-}
