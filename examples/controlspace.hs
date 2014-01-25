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

-- these hints are needed
d0 = _0 :: Dimensionless Double
dt0 = (0 :: Double) *~ (second^neg1)
dtm0 = (0 :: Double) *~ (second^neg1 * meter^neg1)

a = [matD| dt0, _1, _0, _0;
           _0, a22, a23, _0;
           dtm0, _0, _0, _1;
           _0, a42, a43, _0 |]

b = [matD| (0 :: Double) *~ (second / kilo gram);
           1.8182 *~ ((kilo gram)^neg1);
           0.0 *~ (second / (meter * kilo gram));
           4.5455 *~ (meter * kilo gram)^neg1 |]

c = [matD| _1, (0 :: Double) *~ (second), 0 *~ meter, 0 *~ (meter *second);
           _0, _0, _1, _0 |]

d = [matD| (0 :: Double) *~ (meter / newton);
           (0 :: Double) *~ (newton^neg1) |]

-- example state value
x = [matD| 1.0 *~ meter;
           0.2 *~ (meter / second);
           d0;
           0.1 *~ (second^neg1) |]

-- example control input
u = [matD| (0 :: Double) *~ newton |]

xDot = (a `multiply` x) `add` (b `multiply` u)
y = (c `multiply` x) `add` (d `multiply` u)

-- returns the shape of a matrix that, when right-multiplied by a column vector whose dimensions are from, produces a column vector whose dimensions are to
type family DivideVectors (to :: [*]) (from :: [*]) :: [[*]]
-- TODO: implement this type family
-- does it need to be converted to constraint form for inference reasons?

-- draft at creating a type to wrap all this (and simplify creating the matrices...)
-- need to make type functions to compute these types from the xs, ys, us that determine them
-- this type doesn't take a parameter for the dimension of the variable that you are integrating over, because it is pretty much always DTime, and it's probably confusing to generalize it?
-- also it is called an LTI system (emphasis on the T) and not a linear-and-integration-variable-invariant system
data ContinuousLtiSystem (xs :: [*]) (ys :: [*]) (us :: [*]) t = LtiSystem 
                                                                 {
                                                                   a' :: DimMat (DivideVectors (MapDiv DTime xs) xs) t,
                                                                   b' :: DimMat (DivideVectors (MapDiv DTime xs) us) t,
                                                                   c' :: DimMat (DivideVectors ys xs) t,
                                                                   d' :: DimMat (DivideVectors ys us) t
                                                                 }
deriving instance
    (PPUnits (DivideVectors ys us),
     PPUnits (DivideVectors ys xs),
     PPUnits (DivideVectors (MapDiv DTime xs) xs),
     PPUnits (DivideVectors (MapDiv DTime xs) us),
     Show t)
    => Show (ContinuousLtiSystem xs ys us t) 

type ExampleSystem = ContinuousLtiSystem '[DLength, DVelocity, DOne, DFrequency] '[DLength, DOne] '[DForce] Double

pendulum :: ExampleSystem
pendulum = undefined -- if we put definitions of a, b, c, d here where their types are already determined, we wouldn't need to litter them with so many freaking annotations
