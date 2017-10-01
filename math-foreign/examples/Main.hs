module Main where

import Foundation
import Foundation.Class.Storable

import Foreign.Storable (poke)
import Foreign.Marshal.Alloc (alloca)

-- import qualified Math.Linear.Matrix.Mutable as MM
import Math.Complex
import Math.Linear.Matrix
import Math.Linear.Linalg

-- app1 :: IO ()
-- app1 =
--     print $
--     simpsonWithConfig
--         (defaultConf
--          { iteration = 20
--          })
--         (recip . sqrt . sin)
--         (1e-8, 1)

-- app2 :: IO ()
-- app2 = do
--     m <- MM.ones 3 4 :: IO (MM.IOMat Float)
--     MM.shift 4 m
--     MM.times 3 m
--     print m

-- app3 :: IO ()
-- app3 = do
--     let a = M.transpose . M.transpose $ (M.ones 4 3 :: M.Mat Float)
--     print (M.map (+ 1) (M.map (+ 2) (M.ones 3 3 :: M.Mat Float)))
--     print (M.times 3 (M.times 4 a))

-- app4 :: IO ()
-- app4 = do
--     let a = M.identity 4 4 :: M.Mat Float
--     print a
--     print a
--     print a
--     print (a + a)
--     let d = M.dot a a
--     print d
--     print d
--     print d
--     print (M.dot a a)
--     print (M.dot a a)


app5 :: IO ()
app5 = do
    -- x <- ones 4 3 `dot` ones 4 3
    -- let x = ones 4 4 :: Mat (Complex Float)
    --     y = ones 4 4 :: Mat (Complex Float)
    -- putStrLn $ show $ x `dot` y
    -- putStrLn $ show $ det x
    -- putStrLn $ show $ det (x `dot` y)
    let x = fromList' 3 3 [1.0, 1.0, 1.0, 3.0, 2.0, 3.0, 4.0, 5.0, 7.0] :: Mat (Complex Float)
    putStrLn $ show $ (x `at` (1, 2))
    putStrLn $ show $ det x

main :: IO ()
main = do
    app5
