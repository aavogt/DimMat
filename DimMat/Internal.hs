{-# LANGUAGE OverlappingInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
-- similar to Doug McLean's work
module DimMat.Internal (
   -- * hmatrix
   (@@>),
   multiply,
   trans,
   inv,
   det,
   expm,
   -- ** H.Container 
   add,
   sub,
   scale, scaleRecip,
   equal,
   cmap,
   rank,
   mul,
   divide,
   arctan2,
   scalar,
   konst,
   conj,
   addConstant,

   -- *** dimension
   cols, rows,
   colsNT, rowsNT,
   -- * AD
   diff,

   -- * actually internal
   toDM,
   DimMatFromTuple,
   DimMat(..),
  ) where
import GHC.Exts (Constraint)
import qualified Numeric.AD as AD
import qualified Numeric.AD.Types as AD
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.Units.Dimensional.TF
import qualified Prelude as P
import qualified Numeric.NumType.TF as N

import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.LAPACK as H

import Data.Proxy
import Text.PrettyPrint.ANSI.Leijen
import Data.List (transpose)

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
    where unD (Dimensional a) = a


{- ???
grad :: (Num a, Storable (AD.AD t a)) =>
        (forall s. AD.Mode s => DimVec r (AD.AD s a)
                             -> Quantity c (AD.AD s a))
     -> DimVec r a
     -> DimMat '[c] r a
grad f x = DimMat (H.fromLists [AD.grad (f . DimVec . H.fromList) (H.toList x)])
-}


-- | convention is that `head colUnits` should be 'DOne'. That will ensure
-- that we only have one type for a given matrix. Use a GADT for type-signature
-- convenience: going back to a newtype and adding constraints on all the functions
-- for DimMat might happen.
data DimMat (sh :: [[*]]) a where
     DimMat :: (H.Container H.Matrix a, H.Field a)
        => H.Matrix a -> DimMat [ri, DOne ': ci] a
     DimVec :: H.Vector a -> DimMat '[sh] a

-- very crude
instance (Show a, PPUnits sh) => Show (DimMat sh a) where
    showsPrec _ (DimMat m) = case ppUnits (Proxy :: Proxy sh) of
        [rs,cs] -> 
            displayS $
            renderPretty 0.1 80 $ vcat $
            map (hsep . onHead dullgreen) $
            transpose $ map (onHead dullgreen . map string . pad) $
            zipWith (\a b -> a:b)
                ((show (H.rows m) ++ "><"++ show (H.cols m)) : cs) $
            transpose $
            zipWith (\r e -> r : map show e) rs (H.toLists m)
        where
            onHead f (x:xs) = f x : xs
            onHead _ [] = []
            pad :: [String] -> [String]
            pad [] = []
            pad xs = let
                w = maximum (map length xs)
                in map (\x -> take w $ x ++ replicate w ' ') xs


class PPUnits (sh :: [[*]]) where
    ppUnits :: Proxy sh -> [[String]]
instance (PPUnits' x, PPUnits xs) => PPUnits (x ': xs) where
    ppUnits _ = ppUnits' (Proxy :: Proxy x) : ppUnits (Proxy :: Proxy xs)
instance PPUnits '[] where
    ppUnits _ = []

class PPUnits' (sh :: [*]) where
    ppUnits' :: Proxy sh -> [String]
instance (PPUnits' xs) => PPUnits' (DOne ': xs) where
    ppUnits' _ = "1" : ppUnits' (error "ppUnits'" :: Proxy xs)
instance (Show x, PPUnits' xs) => PPUnits' (x ': xs) where
    ppUnits' _ = show (error "ppUnits'" :: x) : ppUnits' (error "ppUnits'" :: Proxy xs)
instance PPUnits' '[] where
    ppUnits' _ = []
-- | put the DimMat into a canonical form, which has a DOne as the first
-- element of the colUnits (2nd) list. Instances for kinds @*@ and @* -> *@
-- 
-- should not be needed, since the GADT does not allow making non-canonical
-- DimMats
type family Canon (a :: k) :: k
type instance Canon (DimMat [r,c]) = DimMat [MapMul (Head c) r, MapDiv (Head c) c]
type instance Canon (DimMat [r,c] x) = DimMat [MapMul (Head c) r, MapDiv (Head c) c] x

-- | @\a xs -> map (map (const a)) xs@
type family MapMapConst (a::k) (xs :: [[l]]) :: [[k]]
type instance MapMapConst a (x ': xs) = MapConst a x ': MapMapConst a xs
type instance MapMapConst a '[] = '[]

-- | @\a xs -> map (const a) xs@
type family MapConst (a :: k) (xs :: [l]) :: [k]
type instance MapConst a (x ': xs) = a ': MapConst a xs
type instance MapConst a '[] = '[]

-- | \a xs -> map (/a) xs
type family MapDiv (a :: *) (xs :: [*]) :: [*]
type instance MapDiv a (x ': xs) = Div x a ': MapDiv a xs
type instance MapDiv a '[] = '[]

-- | map fst
type family MapFst (a :: *) :: [*]
type instance MapFst ((a,_t) , as) = a ': MapFst as
type instance MapFst () = '[]

type family Fst (a :: *) :: *
type instance Fst (a,b) = a
type family Snd (a :: *) :: *
type instance Snd (a,b) = b

-- | convert from (a,(b,(c,(d,())))) to '[a,b,c,d]
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
        DimMat [UnDQuantity (MapFst ijs),
               DOne ': MapDiv (UnDQuantity1 (Fst (Fst ijs)))
               (UnDQuantity (FromPairs (Snd (Fst ijs))))]

-- | just for types
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

-- | xs*ys = zs: knowing two (one that is not zero) will give the third
type family ZipWithMul (xs :: [k]) (ys :: [k]) (zs :: [k]) :: Constraint
type instance ZipWithMul (x ': xs) (y ': ys) (z ': zs)
        = (MultEq x y z, ZipWithMul xs ys zs)
type instance ZipWithMul '[] '[] '[] = ()

-- | zipWith (zipWith Mul) xs ys ~ zs
type family ZipWithZipWithMul (xs :: [[k]]) (ys :: [[k]]) (zs :: [[k]]) :: Constraint
type instance ZipWithZipWithMul (x ': xs) (y ': ys) (z ': zs)
    = (ZipWithMul x y z, ZipWithZipWithMul xs ys zs)
type instance ZipWithZipWithMul '[] '[] '[] = ()

-- | 'product'
type family Product (a :: [k]) :: k
type instance Product (a ': as) = Mul a (Product as)
type instance Product '[] = DOne

-- | @map recip@
type family MapRecip (a :: [k]) :: [k]
type instance MapRecip (a ': as) = Div DOne a ': MapRecip as
type instance MapRecip '[] = '[]

-- | 'length'
type family Len (a :: [*]) :: *
type instance Len '[] = N.Z
type instance Len (a ': as) = N.S (Len as)

-- | @At a n@ is the type-level version of @a !! n@
type family At (a :: [k]) n :: k
type instance At (a ': as) N.Z = a
type instance At (a ': as) (N.S n) = At as n

-- | @AtEq a n b m c@ calculates @(At a n `Mult` At b m) ~ c@,
-- but can also 
type family AtEq (a :: [*]) n (b :: [*]) m (c :: *) :: Constraint
type instance AtEq (a ': as) N.Z (b ': bs) N.Z c = (MultEq a b c)
type instance AtEq (a ': as) (N.S n) bs m c = AtEq as n bs m c
type instance AtEq as N.Z (b ': bs) (N.S m) c = AtEq as N.Z bs m c

-- | multiplication with information going in any direction
type MultEq a b c = (Mul a b ~ c, Div c a ~ b, Div c b ~ a)

type family Head (a :: [k]) :: k
type instance Head (a ': as) = a

type family Tail (a :: [k]) :: [k]
type instance Tail (a ': as) = as

(@@>) :: (N.NumType i, N.NumType j, AtEq ri i ci j ty)
    => DimMat [ri,ci] a
    -> (i, j)
    -> Quantity ty a
DimMat m @@> (i,j) = Dimensional (m H.@@> (N.toNum i,N.toNum j))

multiply :: (H.Product a,
            sh' ~ [MapMul (Inner ci rj) ri, cj])
    => DimMat [ri,ci] a -> DimMat [rj,cj] a
    -> DimMat sh' a
multiply (DimMat a) (DimMat b) = DimMat (H.multiply a b)

trans :: (one ~ DOne,
         sh  ~ [a11 ': ri, one ': ci],
         sh' ~ [a11 ': ci, one ': ri])
         => DimMat sh a -> DimMat sh' a
trans (DimMat a) = DimMat (H.trans a)

type family AreRecips (a :: [k]) (b :: [k]) :: Constraint
type instance AreRecips (a ': as) (b ': bs) =
        (a ~ Div DOne b, b ~ Div DOne a, AreRecips as bs)
type instance AreRecips '[] '[] = ()

type SameLength a b = (SameLength' a b, SameLength' b a)

-- | copied from HList
class SameLength' es1 es
instance (es2 ~ '[]) => SameLength' '[] es2
instance (SameLength' xs ys, es2 ~ (y ': ys)) => SameLength' (x ': xs) es2

-- | if any [k] in the list's length is known, then all other [k] lists in the outer list
-- will be forced to have the same length
type family SameLengths (a :: [[k]]) :: Constraint
type instance SameLengths (a ': b ': bs) = (SameLength a b, SameLengths (b ': bs))
type instance SameLengths '[b] = ()
type instance SameLengths '[] = ()


type family InvCxt (sh :: [[*]]) (sh' :: [[*]]) :: Constraint
       
type instance InvCxt
    [a11 ': ri, dOne ': ci]
    [a11' ': ri', dOne' ': ci'] =
        (SameLengths [ri,ci,ri',ci'], AreRecips ri ci', AreRecips ci ri',
        dOne ~ DOne, dOne' ~ DOne, a11 ~ Div DOne a11', a11' ~ Div DOne a11)

inv :: (InvCxt sh sh', sh' ~ [ri2 ': _1 , DOne ': ci2]) => DimMat sh a -> DimMat sh' a
inv (DimMat a) = DimMat (H.inv a)

det :: (SameLengths [ri,ci]) => DimMat [ri,ci] a
        -> Quantity (Product ri `Mul` Product ci) a
det (DimMat a) = Dimensional (H.det a)

-- | multiplying by a scalar, like 'H.scalar' from hmatrix's 'H.Container' class
type ScaleCxt time ri ri' =
    (MapMul time ri ~ ri',
     MapDiv time ri' ~ ri,
     Head ri' `Div` Head ri ~ time)

{- | 'H.expm'

@y t = expm (scale t a) `multiply` y0@ solves the DE @y' = Ay@ where y0 is the
value of y at time 0

-}
expm :: (MapRecip ci ~ ri, MapRecip ri ~ ci)
    => DimMat [ri,ci] a
    -> DimMat [ri,ci] a
expm (DimMat a) = DimMat (H.expm a)

{- | 'H.scale'

-}
scale :: (ScaleCxt e ri ri')
    => Quantity e a -> DimMat [ri,ci] a -> DimMat [ri',ci] a
scale (Dimensional t) (DimMat a) = DimMat (H.scale t a)

{- | 'H.scaleRecip'
-}
scaleRecip :: (ScaleCxt e' ri ri', Div DOne e ~ e', Div DOne e' ~ e)
    => Quantity e a -> DimMat [ri,ci] a -> DimMat [ri',ci] a
scaleRecip (Dimensional t) (DimMat a) = DimMat (H.scale t a)

liftH2 :: (m ~ DimMat [ri,ci] a,
          h ~ H.Matrix a) => (h -> h -> h) -> m -> m -> m
liftH2 f (DimMat a) (DimMat b) = DimMat (f a b)

add a b = liftH2 H.add a b
sub a b = liftH2 H.sub a b

mul :: (ZipWithZipWithMul sh sh' sh'',
       sh'' ~ [_1 ': _2, DOne ': _3]) => DimMat sh a -> DimMat sh' a -> DimMat sh'' a
mul (DimMat a) (DimMat b) = DimMat (H.mul a b)

divide :: (ZipWithZipWithMul sh' sh'' sh,
          sh'' ~ [_1 ': _2,DOne ': _3]) => DimMat sh a -> DimMat sh' a -> DimMat sh'' a
divide (DimMat a) (DimMat b) = DimMat (H.divide a b)

arctan2 :: (sh' ~ MapMapConst DOne sh) => DimMat sh a -> DimMat sh a -> DimMat sh' a
arctan2 (DimMat a) (DimMat b) = DimMat (H.arctan2 a b)

equal :: (m ~ DimMat [ri,ci] a) => m -> m -> Bool
equal (DimMat a) (DimMat b) = H.equal a b

-- how to nicely defunctionalize? See Fun' Fun in HList?
cmap :: (ScaleCxt deltaE ri ri')
    => proxy deltaE
    -> (forall l m t i th n j e. 
        (N.NumType l, N.NumType m, N.NumType t, N.NumType i,
         N.NumType n, N.NumType j,
         e ~ Dim l m t i th n j) =>
           Quantity e a -> Quantity (Mul e deltaE) b)
    -> DimMat [ri,ci] a -> DimMat [ri',ci] b
cmap _ f m = error "cmap not implemented"
    {- H.mapMatrixWithIndex
     will be useful. 
    -}

rank (DimMat a) = H.rank a
rows (DimMat a) = H.rows a
cols (DimMat a) = H.cols a

-- use proxy?
rowsNT :: DimMat [ri,ci] a -> Len ri
rowsNT = error "rowsNT"

colsNT :: DimMat [ri,ci] a -> Len ci
colsNT = error "colsNT"

scalar :: (H.Field a,
          sh ~ ['[u], '[DOne]]) => Quantity u a -> DimMat sh a
scalar (Dimensional a) = DimMat (H.scalar a)

konst :: forall sh u us ones a _1. (H.Field a,
                              N.NumTypeI (Len ones),
                              N.NumTypeI (Len us),
                              ones ~ (DOne ': _1))
    => Quantity u a -> DimMat [us, ones] a
konst (Dimensional a) = DimMat (H.konst a
    (N.toNum (undefined :: Len us),
     N.toNum (undefined :: Len ones)))

type family CanAddConst (a :: k) (m :: [[k]]) :: Constraint
type instance CanAddConst a [as, ones] = (AllEq a as, AllEq DOne ones)

type family AllEq (a :: k) (xs :: [k]) :: Constraint
type instance AllEq a (x ': xs) = (a ~ x, AllEq a xs)
type instance AllEq a '[] = ()

addConstant :: (H.Field a, CanAddConst u sh)
    => Quantity u a
    -> DimMat sh a
    -> DimMat sh a
addConstant (Dimensional a) (DimMat b) = DimMat (H.addConstant a b)

conj :: DimMat sh a -> DimMat sh a
conj (DimMat a) = DimMat (H.conj a)

{- TODO:
  H.build ::
    H.IndexOf c
    -> hmatrix-0.15.2.0:Numeric.ContainerBoot.ArgOf c e -> c e
  H.atIndex :: c e -> H.IndexOf c -> e
  H.minIndex :: c e -> H.IndexOf c
  H.maxIndex :: c e -> H.IndexOf c
  H.minElement :: c e -> e
  H.maxElement :: c e -> e
  H.sumElements :: c e -> e
  H.prodElements :: c e -> e
  H.step :: H.RealElement e => c e -> c e
  H.cond :: H.RealElement e => c e -> c e -> c e -> c e -> c e -> c e
  H.find :: (e -> Bool) -> c e -> [H.IndexOf c]
  H.assoc :: H.IndexOf c -> e -> [(H.IndexOf c, e)] -> c e
  H.accum :: c e -> (e -> e -> e) -> [(H.IndexOf c, e)] -> c 

  (Conjugate transpose for complex matrices?)
  Pseudoinverse and related decompositions by LAPACK (SVD, QR etc.)
  Friendly syntax for introduction (more matD)
  Friendly syntax for elimination (matD as pattern)
  A pretty-printer that includes the types of each entry
  A clean way to get the n x n dimensionless identity matrix

  check that all types that could could be inferred are
-}
