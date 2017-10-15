{-# OPTIONS_GHC -fplugin Plugin.Intuition -fplugin-opt Plugin.Intuition:--nodebug #-}
{-# OPTIONS_GHC -fplugin Plugin.Intuition.KnownNat #-}
{-# OPTIONS_GHC -fprint-potential-instances #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLists #-}

module Main where

import Foundation
import qualified Math.Linear.Matrix.Mutable as MM
import Math.Complex
import Math.Linear.Vector
import Math.Linear.Matrix
import Math.Linear.Linalg

import Foundation
import Foundation.Primitive

import Plugin.Intuition.KnownNat

app2 :: IO ()
app2 = do
    m <- MM.ones Proxy Proxy :: IO (MM.IOMat Float 3 4)
    MM.shift' 4 m
    MM.times' 3 m
    -- putStrLn . show $ m

-- app3 :: IO ()
-- app3 = do
--     let a = M.transpose . M.transpose $ (M.ones 4 3 :: M.Mat Float)
--     print (M.map (+ 1) (M.map (+ 2) (M.ones 3 3 :: M.Mat Float)))
--     print (M.times 3 (M.times 4 a))

app4 :: IO ()
app4 = do
    let a = identity Proxy Proxy :: Mat Float 4 4
    putStrLn . show $ a
    putStrLn . show $ a `add` a
    let d = dot a a
    putStrLn . show $ d
    putStrLn . show $ a `dot` a

app5 :: IO ()
app5 = do
    let x = ones Proxy Proxy :: Mat (Complex Float) 4 4
        y = ones Proxy Proxy :: Mat (Complex Float) 4 4
    putStrLn . show $ x `dot` y
    putStrLn . show $ det x
    putStrLn . show $ det (x `dot` y)
    let z = fromList' 3 3 [1.0, 1.0, 1.0, 3.0, 2.0, 3.0, 7.0, 8.0, 9.0] :: Mat Float 3 3 -- (Complex Float)
    putStrLn . show $ (z `at` (Proxy :: Proxy 1, Proxy :: Proxy 2))
    putStrLn . show $ det z
    putStrLn . show $ rank z
    putStrLn . show $ inverse z

app6 :: IO ()
app6 = do
    let x = fromList' 3 2 [ 1.54818221,  0.42103339, -0.18791299,  0.72404307, -0.99845554, -0.2456467 ] :: Mat Float 3 2
        y = transpose x
        (qx, rx) = qr x
        -- (qy, ry) = qr y
        -- (qx', rx') = qr' x
        -- (qy', ry') = qr' y
    -- return ()
    putStrLn . show $ qx
    -- putStrLn . show $ rx
    -- putStrLn . show $ qy
    -- putStrLn . show $ ry
    -- putStrLn . show $ qx'
    -- putStrLn . show $ rx'
    -- putStrLn . show $ qy'
    -- putStrLn . show $ ry'

app7 :: IO ()
app7 = do
    let x = fromList' 2 2 [1, 1, 0, 1] :: Mat Float 2 2
        y = [1, 1] :: Vec Float 2
    putStrLn . show $ transform x y

-- app8 :: IO ()
-- app8 = do
--     let x = fromList' 3 2 [ 1.54818221,  0.42103339, -0.18791299,  0.72404307, -0.99845554, -0.2456467 ] :: Mat Float 3 2
--         (u, s, vt) = svd x
--         (u', s', vt') = svd' x
--     putStrLn . show $ u
--     putStrLn . show $ s
--     putStrLn . show $ vt
--     putStrLn . show $ u'
--     putStrLn . show $ s'
--     putStrLn . show $ vt'
--     putStrLn . show $ singular x

app9 :: IO ()
app9 = do
    -- let x = fromList' 2 2 [ 1.54818221,  0.42103339, -0.18791299,  0.72404307 ] :: Mat Float
    let x = fromList' 2 2 [ 4, 3, 7, 9 ] :: Mat Float 2 2
        y = fromList' 2 2 [ 4, 3, 7, 9 ] :: Mat (Complex Float) 2 2
        (lam, ut, v) = eigen x
    putStrLn . show $ lam
    -- putStrLn . show $ ut
    putStrLn . show $ ut
    putStrLn . show $ v
    putStrLn . show $ y `dot` v
    putStrLn . show $ eigenvals x

main :: IO ()
main = do
    app9
