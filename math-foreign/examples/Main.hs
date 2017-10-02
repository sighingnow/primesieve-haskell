{-# LANGUAGE OverloadedLists #-}

module Main where

import Foundation
import qualified Math.Linear.Matrix.Mutable as MM
import Math.Complex
import Math.Linear.Matrix
import qualified Math.Linear.Matrix as M
import Math.Linear.Linalg

app2 :: IO ()
app2 = do
    m <- MM.ones 3 4 :: IO (MM.IOMat Float)
    MM.shift 4 m
    MM.times 3 m
    -- putStrLn . show $ m

-- app3 :: IO ()
-- app3 = do
--     let a = M.transpose . M.transpose $ (M.ones 4 3 :: M.Mat Float)
--     print (M.map (+ 1) (M.map (+ 2) (M.ones 3 3 :: M.Mat Float)))
--     print (M.times 3 (M.times 4 a))

app4 :: IO ()
app4 = do
    let a = M.identity 4 4 :: M.Mat Float
    putStrLn . show $ a
    putStrLn . show $ (a `add` a)
    let d = M.dot a a
    putStrLn . show $ d
    putStrLn . show $ (M.dot a a)

app5 :: IO ()
app5 = do
    let x = ones 4 4 :: Mat (Complex Float)
        y = ones 4 4 :: Mat (Complex Float)
    putStrLn . show $ x `dot` y
    putStrLn . show $ det x
    putStrLn . show $ det (x `dot` y)
    let z = fromList' 3 3 [1.0, 1.0, 1.0, 3.0, 2.0, 3.0, 7.0, 8.0, 9.0] :: Mat Float -- (Complex Float)
    putStrLn . show $ (z `at` (1, 2))
    putStrLn . show $ det z
    putStrLn . show $ rank z
    putStrLn . show $ inverse z

app6 :: IO ()
app6 = do
    let x = fromList' 3 2 [ 1.54818221,  0.42103339, -0.18791299,  0.72404307, -0.99845554, -0.2456467 ] :: Mat Float
        y = transpose x
        (qx, rx) = qr x
        (qy, ry) = qr y
        (qx', rx') = qr' x
        (qy', ry') = qr' y
    putStrLn . show $ qx
    putStrLn . show $ rx
    putStrLn . show $ qy
    putStrLn . show $ ry
    putStrLn . show $ qx'
    putStrLn . show $ rx'
    putStrLn . show $ qy'
    putStrLn . show $ ry'

app7 :: IO ()
app7 = do
    let x = fromList' 2 2 [1, 1, 0, 1] :: Mat Float
        y = [1, 1] :: UArray Float
    putStrLn . show $ transform x y

app8 :: IO ()
app8 = do
    let x = fromList' 3 2 [ 1.54818221,  0.42103339, -0.18791299,  0.72404307, -0.99845554, -0.2456467 ] :: Mat Float
        (u, s, vt) = svd x
        (u', s', vt') = svd' x
    putStrLn . show $ u
    putStrLn . show $ s
    putStrLn . show $ vt
    putStrLn . show $ u'
    putStrLn . show $ s'
    putStrLn . show $ vt'
    putStrLn . show $ singular x

main :: IO ()
main = do
    app8
