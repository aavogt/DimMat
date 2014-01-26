{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
module T1 where
import DimMat
import Numeric.Units.Dimensional.TF.Prelude
-- import Numeric.LinearAlgebra.Algorithms
import qualified Prelude as P


{- Example from http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace -}

-- not really stated on that page, but if you go back a couple of pages in their derivation
-- you can see that the type of u is a 1x1 matrix whose sole element is a force

a22 = (-0.1818 :: Double) *~ (second^neg1)
a23 =  2.6727 *~ (meter * second^neg2)
a42 = (-0.4545) *~ (second^neg1 * meter^neg1)
a43 = 31.1818 *~ (second^neg2)

-- example state value
x = [matD| 1.0 *~ meter;
           0.2 *~ (meter / second);
           _0 :: Dimensionless Double;
           0.1 *~ (second^neg1) |]

-- example control input
u = [matD| (0 :: Double) *~ newton |]

--xDot = (a `multiply` x) `add` (b `multiply` u)
--y = (c `multiply` x) `add` (d `multiply` u)

{- | data type encoding units required by
http://en.wikibooks.org/wiki/Control_Systems/State-Space_Equations#State-Space_Equations

To refresh:

> dxs/dt = A xs + B us
> ys = C xs + D us

the units of d/dt are dtinv.
-}
data LiSystem dtinv (xs :: [*]) (ys :: [*]) (us :: [*]) e where
    LiSystem ::
        (a ~ [dtinv ': ai,DOne ': aj], -- 11 entry needs units of dtinv
         b ~ [b11   ': bi,DOne ': bj], -- matrices are in the canonical form
         c ~ [c11   ': ci,DOne ': cj],
         d ~ [d11   ': di,DOne ': dj],
         PPUnits a, PPUnits b, PPUnits c, PPUnits d,
         MultiplyCxt a xs dxs, -- Ax ~ dx/dt
         MultiplyCxt b us dxs, -- Bu ~ dx/dt
         MultiplyCxt c xs ys,  -- Cx ~ y
         MultiplyCxt d us ys   -- Du ~ y
         )
        => Unit dtinv e -- ^ tends to be (second^neg1)?
        -> DimMat a e -- ^ A
        -> DimMat b e -- ^ B
        -> DimMat c e -- ^ C
        -> DimMat d e -- ^ D
        -> LiSystem dtinv xs ys us e
-- need a show instance?

{- | usually we are time-invariant

 this type doesn't take a parameter for the dimension of the variable
 that you are integrating over, because it is pretty much always DTime,
 and it's probably confusing to generalize it?

 also it is called an LTI system (emphasis on the T) and not
 a linear-and-integration-variable-invariant system
-}
type LtiSystem = LiSystem (Div DOne DTime)

type ExampleSystem = LtiSystem
    [DLength, DVelocity, DOne, DFrequency]
    [DLength, DOne]
    '[DForce]
    Double

pendulum :: ExampleSystem
pendulum = LiSystem (second^neg1) a b c d
           where
             a = [matD| _0, _1, _0, _0;
                        _0, a22, a23, _0;
                        _0, _0, _0, _1;
                        _0, a42, a43, _0 |]
             b = [matD| (0 :: Double) *~ (second / kilo gram);
                        1.8182 *~ ((kilo gram)^neg1);
                        0.0 *~ (second / (meter * kilo gram));
                        4.5455 *~ (meter * kilo gram)^neg1 |]
             c = [matD| _1, (0 :: Double) *~ (second), 0 *~ meter, 0 *~ (meter *second);
                        _0, _0, _1, _0 |]
             d = [matD| (0 :: Double) *~ (meter / newton);
                        (0 :: Double) *~ (newton^neg1) |]
