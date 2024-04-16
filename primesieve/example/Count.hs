import Math.Prime.FastSieve

-- Test the results of prime count.

testEq :: (Eq a, Show a) => String -> a -> a -> IO ()
testEq message a b = do
    putStr $ message <> show a
    if a /= b
        then putStrLn $ ", but should be: " <> show b
        else putStrLn ""

main :: IO ()
main = do
    putStrLn "Test count of prime numbers"

    testEq "count within 10^8: " (countPrimes 0 100000000) 5761455
    testEq "count within 10^9: " (countPrimes 0 1000000000) 50847534
    testEq "count within 2^32: " (countPrimes 0 4294967296) 203280221
    testEq "count within 10^10: " (countPrimes 0 10000000000) 455052511
    testEq "count within 10^11: " (countPrimes 0 100000000000) 4118054813

