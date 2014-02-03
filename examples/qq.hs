{-# LANGUAGE ViewPatterns, AllowAmbiguousTypes #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE QuasiQuotes #-}
module T1 where
import DimMat
import Data.Dimensions.SI
import Data.Dimensions.Poly
import Data.Dimensions.Unsafe


-- probably there's a nicer way to write  kg m?
x00 = 1 % Meter  .* 1 % kilo Gram
x10 = 1 % Second .* 1 % kilo Gram
x01 = 1 % Meter  .* 1 % Second
x11 = 1 % Second .* 1 % Second


{- |

You get the same thing back (except for different type defaulting)

>>> xs @@> (zero,zero)
1.0 m kg
>>> x00
1.0 m kg

>>> xs @@> (pos1,zero)
1.0 kg s
>>> x10
1.0 kg s

>>> xs @@> (pos1,pos1)
1.0 s^2
>>> x11
1 s^2

>>> xs @@> (zero,pos1)
1.0 m s
>>> x01
1 m s

>>> det xs
0.0 m kg s^2

-}
xs = [matD| x00,x01; x10,x11 |]

{- |

>>> y
1.0 m kg
>>> x00
1.0 m kg

-}
y = case xs of [matD| a,_;_,_ |] -> a

{- |

>>> z
6><6 1   kg^-1 s 1   kg^-1 s 1   kg^-1 s
m kg 1.0 1.0     0.0 0.0     0.0 0.0
kg s 1.0 1.0     0.0 0.0     0.0 0.0
m kg 0.0 0.0     1.0 1.0     0.0 0.0
kg s 0.0 0.0     1.0 1.0     0.0 0.0
m kg 0.0 0.0     0.0 0.0     1.0 1.0
kg s 0.0 0.0     0.0 0.0     1.0 1.0

>>> toHLists z
H[H[1.0 m kg, 1.0 m s, 0.0 m kg, 0.0 m s, 0.0 m kg, 0.0 m s], H[1.0 kg s, 1.0 s^2, 0.0 kg s, 0.0 s^2, 0.0 kg s, 0.0 s^2], H[0.0 m kg, 0.0 m s, 1.0 m kg, 1.0 m s, 0.0 m kg, 0.0 m s], H[0.0 kg s, 0.0 s^2, 1.0 kg s, 1.0 s^2, 0.0 kg s, 0.0 s^2], H[0.0 m kg, 0.0 m s, 0.0 m kg, 0.0 m s, 1.0 m kg, 1.0 m s], H[0.0 kg s, 0.0 s^2, 0.0 kg s, 0.0 s^2, 1.0 kg s, 1.0 s^2]]

-}
z = diagBlock $ hBuild xs xs xs
