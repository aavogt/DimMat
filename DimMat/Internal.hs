{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE RankNTypes #-}
-- similar to Doug McLean's work
module DimMat.Internal (
   -- * hmatrix
   (@@>),
   multiply,
   trans,
   inv,
   det,
   -- * AD
   diff,

   -- * actually internal
   toDM,
   DimMatFromTuple,
   DimMat(..),
  ) where
import Foreign.Storable (Storable)
import GHC.Exts (Constraint)
import qualified Numeric.AD as AD
import qualified Numeric.AD.Types as AD
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.Units.Dimensional.TF
import qualified Prelude as P
import qualified Numeric.NumType.TF as N

import qualified Numeric.LinearAlgebra as H
import qualified Data.Packed.Matrix as H
import qualified Numeric.LinearAlgebra.LAPACK as H

{- |
>>> let ke velocity = velocity*velocity*(1*~kilo gram)
>>> diff ke (3 *~ (metre/second))
6.0 m^-1 kg^-1 s

-}
diff :: (Num a) =>
        (forall s. AD.Mode s => Dimensional v x (AD.AD s a)
                             -> Dimensional v y (AD.AD s a))
        -> Dimensional v x a -> Dimensional v (Div x y) a
diff f z = Dimensional $ AD.diff (unD . f . Dimensional) (unD z)


{- ???
grad :: (Num a, Storable (AD.AD t a)) =>
        (forall s. AD.Mode s => DimVec r (AD.AD s a)
                             -> Quantity c (AD.AD s a))
     -> DimVec r a
     -> DimMat '[c] r a
grad f x = DimMat (H.fromLists [AD.grad (f . DimVec . H.fromList) (H.toList x)])
-}

unD (Dimensional a) = a                                                                               


-- | convention is that `head colUnits` should be 'DOne'. That will ensure
-- that we only have one type for a given matrix
newtype DimMat (rowUnits :: [*]) (colUnits :: [*]) a = DimMat (H.Matrix a)
newtype DimVec (units :: [*]) a = DimVec (H.Vector a)

type family MapDiv (a :: *) (xs :: [*]) :: [*]
type instance MapDiv a (x ': xs) = Div x a ': MapDiv a xs
type instance MapDiv a '[] = '[]

type family MapFst (a :: *) :: [*]
type instance MapFst ((a,_t) , as) = a ': MapFst as
type instance MapFst () = '[]

type family Fst (a :: *) :: *
type instance Fst (a,b) = a
type family Snd (a :: *) :: *
type instance Snd (a,b) = b

type family FromPairs (a :: *) :: [*]
type instance FromPairs (a,b) = a ': FromPairs b
type instance FromPairs () = '[]

type family UnDQuantity (a :: [*]) :: [*]
type instance UnDQuantity (x ': xs) = UnDQuantity1 x ': UnDQuantity xs
type instance UnDQuantity '[] = '[]

type family UnDQuantity1 (a :: *) :: *
type instance UnDQuantity1 (Dimensional DQuantity x t) = x

-- at some point the 't' is accessible?

{- | given @ijs :: [[Quantity a]]@ (except the : and [] constructors are
actually (,) and (), ie. a HList), calculate a @DimMat rowUnits colUnits@,
where the outer product of rowUnits and colUnits gives the units at each
index in the ijs.  The first element of colUnits is DOne.
-}
type family DimMatFromTuple ijs :: * -> *
type instance DimMatFromTuple ijs =
        DimMat (UnDQuantity (MapFst ijs))
               (DOne ': MapDiv (UnDQuantity1 (Fst (Fst ijs)))
               (UnDQuantity (FromPairs (Snd (Fst ijs)))))

toDM :: ijs -> DimMatFromTuple ijs t
toDM = error "toDM"

type family Inner (a :: [k]) (b :: [k]) :: k
type instance Inner (a ': as) (b ': bs) = Mul a b
type instance Inner '[] '[] = '[]

type family InnerCxt (a :: [k]) (b :: [k]) :: Constraint
type instance InnerCxt (a ': b ': c) (x ': y ': z) = (Mul a x ~ Mul b y, InnerCxt (b ': c) (y ': z))
type instance InnerCxt '[a] '[b] = ()
type instance InnerCxt '[] '[] = ()
-- and anything else is an error

type family MapMul (a :: k) (xs :: [k]) :: [k]
type instance MapMul a '[] = '[]
type instance MapMul a (x ': xs) = Mul a x ': MapMul a xs

type family Product (a :: [k]) :: k
type instance Product (a ': as) = Mul a (Product as)
type instance Product '[] = DOne

type family MapRecip (a :: [k]) :: [k]
type instance MapRecip (a ': as) = Div DOne a ': MapRecip as
type instance MapRecip '[] = '[]

type family At (a :: [k]) n :: k
type instance At (a ': as) N.Z = a
type instance At (a ': as) (N.S n) = At as n

-- | put the DimMat into a canonical form, which has a DOne as the first
-- elemen of the colUnits (2nd) list. Instances for kinds @*@ and @* -> *@
type family Canon (a :: k) :: k
type instance Canon (DimMat r c) = DimMat (MapMul (Head c) r) (MapDiv (Head c) c)
type instance Canon (DimMat r c x) = DimMat (MapMul (Head c) r) (MapDiv (Head c) c) x

type family Head (a :: [k]) :: k
type instance Head (a ': as) = a

(@@>) :: (N.NumType i, N.NumType j, Storable a) => DimMat ri rj a
    -> (i, j)
    -> Quantity ( (ri `At` i) `Mul` (rj `At` j) ) a
DimMat m @@> (i,j) = Dimensional (m H.@@> (N.toNum i,N.toNum j))

multiply :: H.Product a => DimMat ri ci a -> DimMat rj cj a
    -> DimMat (MapMul (Inner ci rj) ri) cj a
multiply (DimMat a) (DimMat b) = DimMat (H.multiply a b)

type TransCxt ri ci ri' ci' =
         ( Canon (DimMat ci' ri') ~ DimMat ci' ri',
           Canon (DimMat ci ri) ~ DimMat ri' ci')
trans :: TransCxt ri ci ri' ci' => DimMat ri ci a -> DimMat ri' ci' a
trans (DimMat a) = DimMat (H.trans a)

type InvCxt ri ci ri' ci' =
        ( MapRecip ci ~ ri',
          MapRecip ri ~ ci',
          MapRecip ri' ~ ci,
          MapRecip ci' ~ ri)

inv :: (H.Field a, InvCxt ri ci ri' ci') => DimMat ri ci a -> DimMat ri' ci' a
inv (DimMat a) = DimMat (H.inv a)

det :: H.Field a => DimMat ri ci a
        -> Dimensional v (Product ri `Mul` Product ci) a
det (DimMat a) = Dimensional (H.det a)

-- more needed!
