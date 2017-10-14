-----------------------------------------------------------
-- |
-- module:                      Math.Numeric.Random
-- copyright:                   (c) 2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Random data generator according to some special distribution constraints.
--
{-# LANGUAGE ScopedTypeVariables #-}

module Math.Numeric.Random where

import Foundation
import Foundation.Class.Storable (Ptr)
import Foundation.Collection
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Foreign
import Foundation.Primitive
import Foundation.Random

import GHC.TypeLits

import Prelude (fromIntegral)
import Foreign.C.String (castCharToCChar)
import Foreign.Ptr (nullPtr)
import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import Math.Linear.Vector
import qualified Math.Linear.Matrix.Mutable as Mutable
import Math.Linear.Matrix (Mat(..), unsafeWith, unsafeOp, unsafeUnaryOp, unsafeFactorizeOp, unsafeDecomposeOp)

-- | Random number generator powered by PCG.
randInt32 :: IO Int32
randInt32 = c_pcg32_random

foreign import ccall unsafe "pcg_basic.h pcg32_random" c_pcg32_random
    :: IO Int32

randInt32Arr :: forall n. KnownNat n => IO (Vec Int32 n)
randInt32Arr = do
    arr <- mutNew nlen
    withMutableVPtr arr $ \parr ->
        c_pcg32_random_array (integralDownsize nlen) parr
    unsafeFreeze arr
  where nlen = CountOf . fromIntegral . natVal $ (Proxy :: Proxy n)

foreign import ccall unsafe "pcg_random.h pcg32_random_array" c_pcg32_random_array
    :: Int32 -> Ptr Int32 -> IO ()

randFloat :: IO Float
randFloat = c_pcg32_random_float

foreign import ccall unsafe "pcg_random.h pcg32_random_float" c_pcg32_random_float
    :: IO Float

randFloatArr :: forall n. KnownNat n => IO (Vec Float n)
randFloatArr = do
    arr <- mutNew nlen
    withMutableVPtr arr $ \parr ->
        c_pcg32_random_float_array (integralDownsize nlen) parr
    unsafeFreeze arr
  where nlen = CountOf . fromIntegral . natVal $ (Proxy :: Proxy n)

foreign import ccall unsafe "pcg_random.h pcg32_random_float_array" c_pcg32_random_float_array
    :: Int32 -> Ptr Float -> IO ()

-- | Generate an array between given upper and lower bounds, like numpy's linspace.
linspace :: forall a n. (I.Elem a, KnownNat n)
    => Float -> Float -> Vec a n
linspace start end = unsafePerformIO $ do
    rv <- mutNew nlen
    withMutableVPtr rv $ \prv ->
        I.call $ I.linspace prv start end (integralDownsize nlen)
    unsafeFreeze rv
  where nlen = CountOf . fromIntegral . natVal $ (Proxy :: Proxy n)

-- | Random for standard normal distribution. The given @nlen@ must be an even number.
boxmuller :: forall n. KnownNat n => IO (Vec Float n)
boxmuller = do
    arr <- mutNew nlen
    withMutableVPtr arr $ \parr ->
        c_boxmuller (integralDownsize nlen) parr
    unsafeFreeze arr
  where nlen = CountOf . fromIntegral . natVal $ (Proxy :: Proxy n)

foreign import ccall unsafe "pcg_random.h boxmuller" c_boxmuller
    :: Int32 -> Ptr Float -> IO ()

