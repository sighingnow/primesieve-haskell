module Main where

import Foundation

import Math.Prime.FastSieve

main :: IO ()
main = do
    primesieveVersion >>= putStrLn . ("Version of primesieve: " <>)
    
    ns <- generateNPrimes 10 0 :: IO (UArray Int32)
    putStrLn $ "The first 10 primes: " <> show ns

    putStrLn $ "Count of primes between 100 and 1000000: " <> show (countPrimes 100 1000000)

    putStrLn $ "The 12345678th prime number (from 0): " <> show (nthPrime 12345678 0)
