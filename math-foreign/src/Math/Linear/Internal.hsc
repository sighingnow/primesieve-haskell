-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Internal
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Interactive with blas and lapack library via FFI.
--
{-# OPTIONS_HADDOCK hide #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Math.Linear.Internal where

-- import           Control.Monad    ( when )
-- import           Data.Complex     ( Complex(..) )
-- import           Foreign.Ptr
-- import           Foreign.Storable

import Foundation
import Foundation.Class.Storable
import Foundation.Foreign
import Foundation.Numerical
import Foundation.Primitive


import Control.Monad (when)
import GHC.Num (Num)
import qualified Foreign.Storable as Foreign.Storable

import Math.Complex

class (PrimType a, Num a, Storable a, Foreign.Storable.Storable a) => Elem a where
    -- | Encode type into integer and tranfer it to C library.
    typeid :: a -> Int

instance Elem Float where
    typeid _ = 1
    {-# INLINE typeid #-}

instance Elem Double where
    typeid _ = 2
    {-# INLINE typeid #-}

instance Elem (Complex Float) where
    typeid _ = 3
    {-# INLINE typeid #-}

instance Elem (Complex Double) where
    typeid _ = 4
    {-# INLINE typeid #-}

-- | Handle status code returned by pure C code.
call :: Int -> IO ()
call status = when (status /= 0) $ error ("ffi error, return " <> show status)

{-# INLINE call #-}

-- | Handle status code returned by impure C code.
call' :: IO Int -> IO ()
call' status = status >>= \err -> when (err /= 0) $ error ("ffi error, return " <> show err)

{-# INLINE call' #-}

#let ccall name, args = "foreign import ccall unsafe \"linear.h %s\" c_%s :: Int -> %s\n%s :: forall a . Elem a => %s\n%s = c_%s (typeid (undefined :: a))\n{-# INLINE %s #-}", #name, #name, args, #name, args, #name, #name, #name

#ccall identity,    "Ptr a -> Int -> Int -> Int"
#ccall random_,     "Ptr a -> Int -> Int -> IO Int"
#ccall diag,        "Ptr a -> Int -> Int -> Ptr a -> Int"
#ccall diagonal,    "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall sum,         "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall product,     "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall mean,        "Ptr a -> Ptr a -> Int -> Int -> Int"

#ccall lower,       "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall upper,       "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall transpose,   "Ptr a -> Ptr a -> Int -> Int -> Int"

#ccall shift,       "Ptr a -> Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall times,       "Ptr a -> Ptr a -> Ptr a -> Int -> Int -> Int"

#ccall add,         "Ptr a -> Int -> Int -> Int -> Ptr a -> Ptr a -> Int"
#ccall minus,       "Ptr a -> Int -> Int -> Int -> Ptr a -> Ptr a -> Int"
#ccall mult,        "Ptr a -> Int -> Int -> Int -> Ptr a -> Ptr a -> Int"
#ccall division,    "Ptr a -> Int -> Int -> Int -> Ptr a -> Ptr a -> Int"
#ccall dot,         "Ptr a -> Int -> Int -> Int -> Ptr a -> Ptr a -> Int"

#ccall inner,       "Ptr a -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Int"

#ccall det,         "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall trace,       "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall rank,        "Ptr Int -> Ptr a -> Int -> Int -> Int"
#ccall norm,        "Ptr a -> Ptr a -> Int -> Int -> Int"
#ccall inverse,     "Ptr a -> Ptr a -> Int -> Int -> Int"

#ccall lu,          "Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int"
#ccall qr,          "Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int"
#ccall svd,         "Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int"
#ccall jordan,      "Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int"
#ccall cholesky,    "Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int"
#ccall schur,       "Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int -> Int -> Ptr a -> Int"

#ccall transform,   "Ptr a -> Int -> Int -> Ptr a -> Ptr a -> Int"
