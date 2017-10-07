-----------------------------------------------------------
-- |
-- module:                      Math.Numeric.Special
-- copyright:                   (c) 2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Special numerical functions via FFI.
--

module Math.Numeric.Special where

import Foundation

{------------------------------------------------------------------------------
\[
  erf(x) = \frac{1}{\sqrt{\pi}} \int_{ -x}^{x} e^{ -t^2} \,dt
\]
------------------------------------------------------------------------------}

foreign import ccall unsafe "math.h erff" c_erff :: Float -> Float

erf :: Float -> Float
erf = c_erff

{-# INLINABLE erf #-}

foreign import ccall unsafe "math.h erf" c_erf :: Double -> Double

erf' :: Double -> Double
erf' = c_erf

{-# INLINABLE erf' #-}

{------------------------------------------------------------------------------
\[
  tgamma: tgamma(x) = \int_{0}^{\infty} t^{x-1} e^{ -t} \,dt
\]
------------------------------------------------------------------------------}

foreign import ccall unsafe "math.h tgammaf" c_tgammaf :: Float -> Float

tgamma :: Float -> Float
tgamma = c_tgammaf

{-# INLINABLE tgamma #-}

foreign import ccall unsafe "math.h tgamma" c_tgamma :: Double -> Double

tgamma' :: Double -> Double
tgamma' = c_tgamma

{-# INLINABLE tgamma' #-}
