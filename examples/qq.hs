{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE QuasiQuotes #-}
module T1 where
import DimMat
import qualified Prelude as P


x00 = (1::Double) *~ (metre * kilo gram)
x10 = 1 *~ (second * kilo gram)
x01 = 1 *~ (metre * second)
x11 = 1 *~ (second * second)


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

