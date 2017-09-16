-- test for Algo.Data

module Main where

import Foundation
import Foundation.Check
import Foundation.Check.Main

import Algo.Data.Heap.MaxiphobicTest

main :: IO ()
main = defaultMain $ Group "Test: Algo.Data"
  [ maxiphobicTest
  ]
