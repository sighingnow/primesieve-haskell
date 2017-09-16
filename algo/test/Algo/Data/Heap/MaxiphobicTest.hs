{-# LANGUAGE DataKinds #-}

module Algo.Data.Heap.MaxiphobicTest
  ( maxiphobicTest
  ) where

import Foundation
import Foundation.Check

import Algo.Data.Heap
import Algo.Data.Heap.Maxiphobic

maxiphobicTest :: Test
maxiphobicTest = Group "maxiphobicTest"
  [ testSize
  ]

testSize :: Test
testSize = Group "testSize"
  [ Property "size of nil is 0: " $ size heap1 === 0
  , Property "size of singleton is 1: " $ size heap2 === 1
  , Property "size of merge: " $ size heap3 === size heap2 + size heap2
  ]

heap1 :: MaxiphobicHeap Int 0
heap1 = Nil

heap2 :: MaxiphobicHeap Int 1
heap2 = toHeap 1

heap3 :: MaxiphobicHeap Int 2
heap3 = merge heap2 heap2
