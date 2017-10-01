-----------------------------------------------------------
-- |
-- module:                      Math.Complex
-- copyright:                   (c) 2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Coordinate `Data.Complex` with foundation's numeric and primitive type support.
--
{-# OPTIONS_GHC -Wno-orphans #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UnboxedTuples #-}

module Math.Complex
  ( Complex (..)
  , conjugate
  ) where

import Foundation
import Foundation.Class.Storable
import Foundation.Primitive

import Data.Complex (Complex(..), conjugate)
import GHC.Exts

instance Additive a => Additive (Complex a) where
    azero = azero :+ azero
    (a :+ b) + (c :+ d) = (a + c) :+ (b + d)

instance (Subtractive a, Difference a ~ a) => Subtractive (Complex a) where
    type Difference (Complex a) = (Complex a)
    (a :+ b) - (c :+ d) = (a - c) :+ (b - d)

instance (Additive a, Subtractive a, Difference a ~ a, Multiplicative a) => Multiplicative (Complex a) where
    midentity = midentity :+ azero
    (a :+ b) * (c :+ d) = (a * c - b * d) :+ (b * c + a * d)

instance (Additive a, Subtractive a, Difference a ~ a, Divisible a) => Divisible (Complex a) where
    (a :+ b) / (c :+ d) = ((a * c + b * d) / (c * c + d * d)) :+ ((b * c - a * d) / (c * c + d * d))

offsetComplex :: Offset (Complex a) -> (# Int#, Int# #)
offsetComplex !(Offset (I# i)) = (# n, n +# 1# #)
    where !n = uncheckedIShiftL# i 1#

{-# INLINE offsetComplex #-}

instance PrimType (Complex Float) where
    primSizeInBytes _ = primSizeInBytes (Proxy :: Proxy Float) + primSizeInBytes (Proxy :: Proxy Float)
    {-# INLINE primSizeInBytes #-}

    primShiftToBytes _ = 3 -- TODO may be wrong
    {-# INLINE primShiftToBytes #-}

    primBaUIndex ba n = F# (indexFloatArray# ba n1) :+ F# (indexFloatArray# ba n2)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primBaUIndex #-}

    primMbaURead mba n = primitive $ \s1 -> let !(# s2, r1 #) = readFloatArray# mba n1 s1
                                                !(# s3, r2 #) = readFloatArray# mba n2 s2
                                             in (# s3, F# r1 :+ F# r2 #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primMbaURead #-}

    primMbaUWrite mba n ((F# w1) :+ (F# w2)) = primitive $ \s1 -> let !s2 = writeFloatArray# mba n1 w1 s1
                                                                   in (# writeFloatArray# mba n2 w2 s2, () #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primMbaUWrite #-}

    primAddrIndex addr n = F# (indexFloatOffAddr# addr n1) :+ F# (indexFloatOffAddr# addr n2)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primAddrIndex #-}

    primAddrRead addr n = primitive $ \s1 -> let !(# s2, r1 #) = readFloatOffAddr# addr n1 s1
                                                 !(# s3, r2 #) = readFloatOffAddr# addr n2 s2
                                              in (# s3, F# r1 :+ F# r2 #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primAddrRead #-}

    primAddrWrite addr n ((F# w1) :+ (F# w2)) = primitive $ \s1 -> let !s2 = writeFloatOffAddr# addr n1 w1 s1
                                                                    in (# writeFloatOffAddr# addr n2 w2 s2, () #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primAddrWrite #-}

instance PrimType (Complex Double) where
    primSizeInBytes _ = primSizeInBytes (Proxy :: Proxy Double) + primSizeInBytes (Proxy :: Proxy Double)
    {-# INLINE primSizeInBytes #-}

    primShiftToBytes _ = 5 -- TODO may be wrong
    {-# INLINE primShiftToBytes #-}

    primBaUIndex ba n = D# (indexDoubleArray# ba n1) :+ D# (indexDoubleArray# ba n2)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primBaUIndex #-}

    primMbaURead mba n = primitive $ \s1 -> let !(# s2, r1 #) = readDoubleArray# mba n1 s1
                                                !(# s3, r2 #) = readDoubleArray# mba n2 s2
                                             in (# s3, D# r1 :+ D# r2 #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primMbaURead #-}

    primMbaUWrite mba n ((D# w1) :+ (D# w2)) = primitive $ \s1 -> let !s2 = writeDoubleArray# mba n1 w1 s1
                                                                   in (# writeDoubleArray# mba n2 w2 s2, () #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primMbaUWrite #-}

    primAddrIndex addr n = D# (indexDoubleOffAddr# addr n1) :+ D# (indexDoubleOffAddr# addr n2)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primAddrIndex #-}

    primAddrRead addr n = primitive $ \s1 -> let !(# s2, r1 #) = readDoubleOffAddr# addr n1 s1
                                                 !(# s3, r2 #) = readDoubleOffAddr# addr n2 s2
                                              in (# s3, D# r1 :+ D# r2 #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primAddrRead #-}

    primAddrWrite addr n ((D# w1) :+ (D# w2)) = primitive $ \s1 -> let !s2 = writeDoubleOffAddr# addr n1 w1 s1
                                                                    in (# writeDoubleOffAddr# addr n2 w2 s2, () #)
        where !(# n1, n2 #) = offsetComplex n
    {-# INLINE primAddrWrite #-}

instance Storable (Complex Float) where
    peek (Ptr addr) = primAddrRead addr (Offset 0)
    poke (Ptr addr) = primAddrWrite addr (Offset 0)

instance Storable (Complex Double) where
    peek (Ptr addr) = primAddrRead addr (Offset 0)
    poke (Ptr addr) = primAddrWrite addr (Offset 0)

instance StorableFixed (Complex Float) where
    size _ = size (Proxy :: Proxy Float) + size (Proxy :: Proxy Float)
    alignment _ = alignment (Proxy :: Proxy Float)

instance StorableFixed (Complex Double) where
    size _ = size (Proxy :: Proxy Double) + size (Proxy :: Proxy Double)
    alignment _ = alignment (Proxy :: Proxy Double)
