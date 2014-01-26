{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
module T1 where
import DimMat
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.Units.Dimensional.TF
import qualified Prelude as P
import Text.PrettyPrint.ANSI.Leijen


{- Example from http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace -}

-- not really stated on that page, but if you go back a couple of pages in their derivation
-- you can see that the type of u is a 1x1 matrix whose sole element is a force

massOfCart = (0.5 :: Double) *~ (kilo gram)
massOfPendulum = (0.2 :: Double) *~ (kilo gram)
coefficientOfFrictionForCart = (0.1 :: Double) *~ (newton / (meter / second))
lengthToPendulumCenterOfMass = (0.3 :: Double) *~ meter
massMomentOfInertiaOfPendulum = (0.006 :: Double) *~ (kilo gram * meter^pos2)
g = (9.8 :: Double) *~ (meter / second^pos2)

p = massMomentOfInertiaOfPendulum*(massOfCart+massOfPendulum)+(massOfCart*massOfPendulum*lengthToPendulumCenterOfMass*lengthToPendulumCenterOfMass)

a22 = negate (massMomentOfInertiaOfPendulum+massOfPendulum * lengthToPendulumCenterOfMass * lengthToPendulumCenterOfMass) * coefficientOfFrictionForCart / p
a23 = (massOfPendulum * massOfPendulum * g * lengthToPendulumCenterOfMass * lengthToPendulumCenterOfMass) / p
a42 = negate (massOfPendulum * lengthToPendulumCenterOfMass * coefficientOfFrictionForCart) / p
a43 = massOfPendulum * g * lengthToPendulumCenterOfMass*(massOfCart + massOfPendulum)/p

b21 = (massMomentOfInertiaOfPendulum + (massOfPendulum * lengthToPendulumCenterOfMass * lengthToPendulumCenterOfMass)) / p
b41 = massOfPendulum * lengthToPendulumCenterOfMass / p

-- example state value
x = [matD| 1.0 *~ meter;
           0.2 *~ (meter / second);
           _0 :: Dimensionless Double;
           0.1 *~ (second^neg1) |]

dx = scale (_1 / (1 *~ second)) x

-- example control input
u = [matD| (0 :: Double) *~ newton |]

--xDot = (a `multiply` x) `add` (b `multiply` u)

y = [matD| 1 *~ meter; _0 :: Dimensionless Double  |]

-- :t isLTI (1 *~ second) x u y, given that x u and y have monomorphic
-- types properly infers constraints on the a,b,c,d arguments!
isLTI time x u y a b c d =
    (scale (_1 /time) x `add` multiply a x `add` multiply b u,
     y `add` multiply c x `add` multiply d u)

testIsLTI =
  (\ a b c d -> case isLTI (1 *~ second) x u y a b c d of
   _ -> do
    print $ vsep
        [text "A = " </> indent 0 (pretty (konst _0 `asTypeOf` a)),
         text "B = " </> indent 0 (pretty (konst _1 `asTypeOf` b)),
         text "C = " </> indent 0 (pretty (konst _2 `asTypeOf` c)),
         text "D = " </> indent 0 (pretty (konst _3 `asTypeOf` d))]
    ) undefined undefined undefined undefined


{- | data type encoding units required by
http://en.wikibooks.org/wiki/Control_Systems/State-Space_Equations#State-Space_Equations

To refresh:

> dxs/dt = A xs + B us
> ys = C xs + D us

the units of d/dt are dtinv.
-}

type family DivideVectors (to :: [*]) (from :: [*]) :: [[*]]
type instance DivideVectors ts fs = [ MapDiv (Head fs) ts, MapMul (Head fs) (MapRecip fs) ]

data LiSystem (iv :: *) (xs :: [*]) (ys :: [*]) (us :: [*]) e = LiSystem
                                                                {
                                                                  a :: DimMat (DivideVectors (MapDiv iv xs) xs) e,
                                                                  b :: DimMat (DivideVectors (MapDiv iv xs) us) e,
                                                                  c :: DimMat (DivideVectors ys xs) e,
                                                                  d :: DimMat (DivideVectors ys us) e
                                                                }
deriving instance
    (PPUnits (DivideVectors (MapDiv iv xs) xs),
     PPUnits (DivideVectors (MapDiv iv xs) us),
     PPUnits (DivideVectors ys xs),
     PPUnits (DivideVectors ys us),
     Show e)
    => Show (LiSystem iv xs ys us e)
-- need a show instance?

{- | usually we are time-invariant

 this type doesn't take a parameter for the dimension of the variable
 that you are integrating over, because it is pretty much always DTime,
 and it's probably confusing to generalize it?

 also it is called an LTI system (emphasis on the T) and not
 a linear-and-integration-variable-invariant system
-}
type LtiSystem = LiSystem DTime

type ExampleSystem = LtiSystem
    [DLength, DVelocity, DPlaneAngle, DAngularVelocity]
    [DLength, DPlaneAngle]
    '[DForce]
    Double

pendulum :: ExampleSystem
pendulum = LiSystem { a = a', b = b', c = c', d = d' }
         where
           a' = [matD| _0, _1, _0, _0;
                       _0, a22, a23, _0;
                       _0, _0, _0, _1;
                       _0, a42, a43, _0 |]
           b' = [matD| _0;
                       b21;
                       _0;
                       b41 |]
           c' = [matD| _1, _0, _0, _0;
                       _0, _0, _1, _0 |]
           d' = [matD| _0;
                       _0 |]

poles x = eigenvalues . a $ x

evaluate sys x u = let xDot = ((a sys) `multiply` x) `add` ((b sys) `multiply` u)
                       y = ((c sys) `multiply` x) `add` ((d sys) `multiply` u)
    in (xDot, y)
