-------------------------------------------------------------------------------
-- |
-- module:          Algo.Data.Heap.Maxiphobic
--
-- Maxiphobic heaps from <Fun with binary heap trees> by Chris Okasaki.

{-# OPTIONS_GHC -fplugin=Plugin.Intuition #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE RankNTypes #-}

module Algo.Data.Heap.Maxiphobic where

import Foundation
import GHC.TypeLits

import Algo.Data.Heap

data MaxiphobicHeap a (s :: Nat) where
    Nil :: Ord a => MaxiphobicHeap a 0
    Fork :: Ord a
         => Int -- ^ keep track of the size of each tree
         -> a -> MaxiphobicHeap a m -> MaxiphobicHeap a n -> MaxiphobicHeap a (1 + m + n)

-- TODO
-- deriving instance Eq a => Eq (MaxiphobicHeap a s)

deriving instance Show a => Show (MaxiphobicHeap a s)

-- deleteMin (Fork n x l r) = merge l r

-- {-# INLINE deleteMin #-}

instance Heap MaxiphobicHeap where
    nil = Nil
    {-# INLINE nil #-}

    toHeap e = Fork 1 e Nil Nil
    {-# INLINE toHeap #-}

    size Nil = 0
    size (Fork n _ _ _) = n
    {-# INLINE size #-}

    isEmpty Nil = True
    isEmpty (Fork _ _ _ _) = False
    {-# INLINE isEmpty #-}

    minElem (Fork _ x _ _) = x
    {-# INLINE minElem #-}

    deleteMin (Fork _ _ l r) = merge l r

    -- | merge: keep the maxiphobic ("biggest avoiding") property.
    merge Nil Nil = Nil
    merge Nil r@(Fork _ _ _ _) = r
    merge l@(Fork _ _ _ _) Nil = l
    merge l@(Fork _ v1 _ _) r@(Fork _ v2 _ _)
        | v1 <= v2 = join l r
        | otherwise = join r l
        where
            -- explicit signature for `join` is necessary.
            join :: forall a (m :: Nat) (n :: Nat). MaxiphobicHeap a m -> MaxiphobicHeap a n -> MaxiphobicHeap a (m + n)
            join (Fork n x ll lr) r
                | size ll <= size lr && size ll <= size r
                    = Fork (n + size r) x ll (merge lr r)
                | size lr <= size ll && size lr <= size r
                    = Fork (n + size r) x lr (merge ll r)
                | size r <= size ll && size r <= size lr
                    = Fork (n + size r) x r (merge ll lr)
            -- TODO needs to improve the type constraint: m > 0.
