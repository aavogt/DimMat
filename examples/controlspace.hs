{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE QuasiQuotes #-}
-- http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace
module T2 where
import DimMat
import qualified Prelude as P

a22 = (-0.1818 :: Double) *~ (second^neg1)
a23 =  2.6727 *~ (meter * second ^neg4)
a42 = (-0.4545) *~ (second / meter)
a43 = 31.1818 *~ (second ^ neg2)

-- these hints might be needed even after fixing matD?
d0 = _0 :: Dimensionless Double
d1 = _1 :: Dimensionless Double


{- |

>>> a0
2><2   1       m s^-3 
s^-1   -0.1818 2.6727 
m^-1 s -0.4545 31.1818

-}
a0 = [matD| a22, a23;
            a42, a43 |]

{- |

>>> a
4><4   1   1       m s^-3  1  
1      0.0 1.0     0.0     0.0
s^-1   0.0 -0.1818 2.6727  0.0
1      0.0 0.0     0.0     1.0
m^-1 s 0.0 -0.4545 31.1818 0.0

-}
a = [matD| d0,  _1,  _0, _0;
           _0, a22, a23, _0;
           d0,  _0,  _0, _1;
           _0, a42, a43, _0 |]


