-----------------------------------------------------------
-- |
-- module:                      Math.Numeric.Random
-- copyright:                   (c) 2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Random data generator according to some special distribution constraints.
--

module Math.Numeric.Random where

import Foundation
import Foundation.Class.Storable (Ptr)
import Foundation.Collection
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Foreign
import Foundation.Primitive
import Foundation.Random

import Foreign.C.String (castCharToCChar)
import Foreign.Ptr (nullPtr)
import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import qualified Math.Linear.Matrix.Mutable as Mutable
import Math.Linear.Matrix (Mat(..), unsafeWith, unsafeOp, unsafeUnaryOp, unsafeFactorizeOp, unsafeDecomposeOp)

-- | Random number generator powered by PCG.
randInt32 :: IO Int32
randInt32 = c_pcg32_random

foreign import ccall unsafe "pcg_basic.h pcg32_random" c_pcg32_random
    :: IO Int32

randInt32Arr :: CountOf Int32 -> IO (UArray Int32)
randInt32Arr nlen = do
    arr <- mutNew nlen
    withMutablePtr arr $ \parr ->
        c_pcg32_random_array (integralDownsize nlen) parr
    unsafeFreeze arr

foreign import ccall unsafe "pcg_random.h pcg32_random_array" c_pcg32_random_array
    :: Int32 -> Ptr Int32 -> IO ()

randFloat :: IO Float
randFloat = c_pcg32_random_float

foreign import ccall unsafe "pcg_random.h pcg32_random_float" c_pcg32_random_float
    :: IO Float

randFloatArr :: CountOf Float -> IO (UArray Float)
randFloatArr nlen = do
    arr <- mutNew nlen
    withMutablePtr arr $ \parr ->
        c_pcg32_random_float_array (integralDownsize nlen) parr
    unsafeFreeze arr

foreign import ccall unsafe "pcg_random.h pcg32_random_float_array" c_pcg32_random_float_array
    :: Int32 -> Ptr Float -> IO ()

-- | Generate an array between given upper and lower bounds, like numpy's linspace.
linspace :: I.Elem a
    => Float -> Float -> Int32 -> UArray a
linspace start end size = unsafePerformIO $ do
    rv <- mutNew (integralCast size)
    withMutablePtr rv $ \prv ->
        I.call $ I.linspace prv start end size
    unsafeFreeze rv

-- | Random for normal distribution.


