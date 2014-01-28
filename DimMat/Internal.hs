{-# LANGUAGE CPP #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}
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
{- | an incomplete extension of dimensional-tf to work with hmatrix and ad
Haddocks are organized to follow hmatrix

note: for subscripting, use the 'HNat' while the units still use 'N.NumType'

TODO:

*  Friendly syntax for introduction (more matD)

*  Friendly syntax for elimination (matD as pattern)

*  A pretty-printer that multiplies out the dimensions at each place?
   The current show instance makes you mentally multiply out row and column units
   (which helps you to see the big picture, but may be more work in other cases)

*  check that all types that could could be inferred are

*  default columns/rows to dimensionless?

*  missing operations (check export list comments)
-}
module DimMat.Internal (
   -- * Data.Packed.Vector
   (@>),
   -- * Data.Packed.Matrix
   -- ** dimension
   cols, rows,
   colsNT, rowsNT,
   hasRows, hasCols,
   -- (><),
   trans,
   -- reshape, flatten, fromLists, toLists, buildMatrix,
   ToHLists(toHLists), toHList,
   FromHLists(fromHLists), fromHList,
   (@@>),

   -- asRow, asColumn, fromRows, toRows, fromColumns, toColumns
   -- fromBlocks
#if MIN_VERSION_hmatrix(0,15,0)
   diagBlock,
#endif
   -- toBlocks, toBlocksEvery, repmat, flipud, fliprl
   -- subMatrix, takeRows, dropRows, takeColumns, dropColumns,
   -- extractRows, diagRect, takeDiag, mapMatrix,
   -- mapMatrixWithIndexM, mapMatrixWithIndexM_, liftMatrix,
   -- liftMatrix2, liftMatrix2Auto, fromArray2D,

   ident, -- where to put this?
   -- * Numeric.Container
   -- constant, linspace,
   diag,
   ctrans,
   -- ** Container class
   scalar,
   conj,
   scale, scaleRecip,
   addConstant,
   add,
   sub,
   mul,
   divide,
   equal,
   arctan2,
   hconcat,
   vconcat,
   cmap,
   konst,
   zeroes,
   -- build, atIndex, minIndex, maxIndex, minElement, maxElement,
   -- sumElements, prodElements, step, cond, find, assoc, accum,
   -- Convert
   -- ** Product class
   multiply,
   -- dot, absSum, norm1, norm2, normInf,
   -- optimiseMult, mXm, mXv, vXm, (<.>),
   -- (<>), (<\>), outer, kronecker,
   -- ** Random numbers
   -- ** Element conversion
   -- ** Input/Output
   -- ** Experimental

   -- * Numeric.LinearAlgebra.Algorithms
   -- | incomplete wrapper for "Numeric.LinearAlgebra.Algorithms"

   -- ** Linear Systems
   -- linearSolve, luSolve, cholSolve, linearSolveLS, linearSolveSVD,
   inv,
   PInv(pinv), 
   pinvTol,
   det,
   -- invlndet,
   rank,
   -- rcond,
   -- ** Matrix factorizations

   -- *** Singular value decomposition
   -- *** Eigensystems
   -- $eigs
   wrapEig, wrapEigOnly,
   EigCxt,
   -- **** eigenvalues and eigenvectors
   eig,
   eigC,
   eigH,
   eigH',
   eigR,
   eigS,
   eigS',
   eigSH,
   eigSH',

   -- **** eigenvalues
   eigOnlyC,
   eigOnlyH,
   eigOnlyR,
   eigOnlyS,
   eigenvalues,
   eigenvaluesSH,
   eigenvaluesSH',

   -- *** QR
   -- *** Cholesky
   -- *** Hessenberg
   -- *** Schur
   -- *** LU 

   -- ** Matrix functions
   -- sqrtm, matFunc
   expm,

   -- ** Nullspace
   -- ** Norms
   -- ** Misc
   -- ** Util 

   -- * Automatic Differentiation
   -- $ad
   diff,

   -- * actually internal
   toDM,
   DimMatFromTuple,
   DimMat(..),
   AtEq, MapMul, Inner, SameLengths, Product, MapRecip,
   ZipWithZipWithMul, MapMapConst, CanAddConst, PPUnits,
   PPUnits', Head, MapDiv, AreRecips, ZipWithMul, PairsToList,
   DiagBlock, MapConst, SameLength', AppendShOf,
   MultiplyCxt, MapMultEq, Trans,
   Tail,MapMultEq', AppendEq, MultEq, Append, AppendEq',
   DropPrefix,
   RmDimensional(RmDimensional),
   AreRecips',ToHList,UnDimMat,AllEq,
   HListFromList,AddDimensional,ToHListRows',
   ToHListRow,AddQty,
   HMapOutWith(..),
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
import qualified Numeric.LinearAlgebra.LAPACK as H

import Text.PrettyPrint.ANSI.Leijen
import Data.List (transpose)

import Data.HList.CommonMain hiding (MapFst)

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


{- $ad

TODO: gradients, hessians, etc.

Types for derivative towers can see hlist's @HList\/Data\/HList\/broken\/Lazy.hs@

Complications include the fact that AD.grad needs a traversable,
but hmatrix stuff is not traversable (due needing Storable). In ipopt-hs
I got around this problem by copying data. Perhaps that is the solution?

> BROKEN
> grad :: (Num a, AreRecips i iinv, H.Element a, Storable a,
>           MapMultEq o iinv r) =>
>         (forall s. (AD.Mode s, H.Container H.Vector (AD.AD s a),
>                     Storable (AD.AD s a), H.Field (AD.AD s a))
>                 => DimMat '[i] (AD.AD s a)
>                 -> Quantity o (AD.AD s a))
>      -> DimMat '[i] a
>      -> DimMat '[r] a
> grad f (DimVec x) = DimMat (H.fromLists [AD.grad (unQty . f . DimVec . H.fromList) (H.toList x)])
>     where unQty (Dimensional a) = a
-}


{- | Matrix with statically checked units (and dimensions). This wraps up
HMatrix. The `sh` type parameter here contains @[row,column]@ units.
The outer product of @row@ and @column@ lists gives a matrix of units of
the same size as the data contained within.

some pain happens to ensure that sh always follows the convention
that the first element of the column units is 'DOne'. One remaining
potential error is that you can have @ri ~ '[]@, which would make for
a 0-row matrix.
-}
data DimMat (sh :: [[*]]) a where
     DimMat :: (H.Container H.Matrix a, H.Field a)
        => H.Matrix a -> DimMat [ri, DOne ': ci] a
     DimVec :: (H.Container H.Vector a, H.Field a)
        => H.Vector a -> DimMat '[sh] a

-- very crude
instance (Show a, PPUnits sh) => Pretty (DimMat sh a) where
    pretty (DimVec v) = case ppUnits (proxy :: Proxy sh) of
        [rs] -> vcat
             [ dullgreen (string label) <+> string (show e)
                | (e,label) <- H.toList v `zip` pad rs ]
    pretty (DimMat m) = case ppUnits (proxy :: Proxy sh) of
        [rs,cs] -> 
            vcat $
            map (hsep . onHead dullgreen) $
            transpose $ map (onHead dullgreen . map string . pad) $
            zipWith (\a b -> a:b)
                ((show (H.rows m) ++ "><"++ show (H.cols m)) : cs) $
            transpose $
            zipWith (\r e -> r : map show e) rs (H.toLists m)
        where
            onHead f (x:xs) = f x : xs
            onHead _ [] = []

instance Pretty (DimMat sh a) => Show (DimMat sh a) where
    showsPrec p x = showsPrec p (pretty x)

pad :: [String] -> [String]
pad [] = []
pad xs = let
    w = maximum (map length xs)
    in map (\x -> take w $ x ++ replicate w ' ') xs


class PPUnits (sh :: [[*]]) where
    ppUnits :: Proxy sh -> [[String]]
instance (PPUnits' x, PPUnits xs) => PPUnits (x ': xs) where
    ppUnits _ = ppUnits' (proxy :: Proxy x) : ppUnits (proxy :: Proxy xs)
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

-- | @\\a xs -> map (map (const a)) xs@
type family MapMapConst (a::k) (xs :: [[l]]) :: [[k]]
type instance MapMapConst a (x ': xs) = MapConst a x ': MapMapConst a xs
type instance MapMapConst a '[] = '[]

-- | @\\a xs -> map (const a) xs@
type family MapConst (a :: k) (xs :: [l]) :: [k]
type instance MapConst a (x ': xs) = a ': MapConst a xs
type instance MapConst a '[] = '[]

-- | @\\a xs -> map (/a) xs@
type family MapDiv (a :: *) (xs :: [*]) :: [*]
type instance MapDiv a (x ': xs) = Div x a ': MapDiv a xs
type instance MapDiv a '[] = '[]

-- | @map fst@
type family MapFst (a :: *) :: [*]
type instance MapFst ((a,_t) , as) = a ': MapFst as
type instance MapFst () = '[]

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
actually (,) and (), ie. a HList that doesn't use the HList constructors),
calculate a @DimMat rowUnits colUnits@, where the outer product of rowUnits
and colUnits gives the units at each index in the ijs.  The first element
of colUnits is DOne.
-}
type family DimMatFromTuple ijs :: * -> *
type instance DimMatFromTuple ijs =
        DimMat [UnDQuantity (MapFst ijs),
               DOne ': MapDiv (UnDQuantity1 (Fst (Fst ijs)))
               (UnDQuantity (FromPairs (Snd (Fst ijs))))]


-- | just for types produced by the matD quasiquote
toDM :: ijs -> DimMatFromTuple ijs t
toDM = error "toDM"

-- | @Inner a b = 'H.dot' a b@
type family Inner (a :: [k]) (b :: [k]) :: k
type instance Inner (a ': as) (b ': bs) = Mul a b
type instance Inner '[] '[] = '[]

-- | @InnerCxt t a b = t ~ 'H.dot' a b@
type family InnerCxt (t :: k) (a :: [k]) (b :: [k]) :: Constraint
type instance InnerCxt t (a ': c) (x ': z) = (MultEq a x t, InnerCxt t c z)
type instance InnerCxt t '[] '[] = ()

-- | MapMul a xs = map (a*) xs
type family MapMul (a :: k) (xs :: [k]) :: [k]
type instance MapMul a '[] = '[]
type instance MapMul a (x ': xs) = Mul a x ': MapMul a xs

-- | MapMulEq a xs ys = map (a*) xs ~ ys
--
-- should work in any direction (given xs and ys, solve a),
-- unlike 'MapMul' that goes forwards only
type MapMultEq a xs ys = (SameLengths [xs,ys], MapMultEq' a xs ys)

-- | helper for 'MapMultEq'
type family MapMultEq' (a :: k) (xs :: [k]) (ys :: [k]) :: Constraint
type instance MapMultEq' a (x ': xs) (y ': ys) = (MultEq a x y, MapMultEq' a xs ys)
type instance MapMultEq' a '[] '[] = ()

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

-- | @AtEq a n b m c@ calculates @(At a n \`Mult\` At b m) ~ c@,
-- but also can infer part of the `a` if the `b` and `c` are known
type family AtEq (a :: [*]) (n :: HNat) (b :: [*]) (m :: HNat) (c :: *) :: Constraint
type instance AtEq (a ': as) HZero (b ': bs) HZero c = (MultEq a b c)
type instance AtEq (a ': as) (HSucc n) bs m c = AtEq as n bs m c
type instance AtEq as HZero (b ': bs) (HSucc m) c = AtEq as HZero bs m c

-- | multiplication with information going in any direction
type MultEq a b c = (Mul a b ~ c, Div c a ~ b, Div c b ~ a)

type family Head (a :: [k]) :: k
type instance Head (a ': as) = a

type family Tail (a :: [k]) :: [k]
type instance Tail (a ': as) = as

-- | Data.Packed.Vector.'H.@>'
(@>) :: (HNat2Integral i)
    => DimMat '[units] a
    -> Proxy i
    -> Quantity (HLookupByHNat i units) a
DimVec v @> i = Dimensional (v H.@> hNat2Integral i)

-- | Data.Packed.Matrix.'H.@@>'
(@@>) :: (HNat2Integral i, HNat2Integral j, AtEq ri i ci j ty)
    => DimMat [ri,ci] a
    -> (Proxy i, Proxy j)
    -> Quantity ty a
DimMat m @@> (i,j) = Dimensional (m H.@@> (hNat2Integral i,hNat2Integral j))

type family MultiplyCxt (sh1 :: [[*]]) (sh2 :: [*]) (sh3 :: [*]) :: Constraint
type instance MultiplyCxt [r11 ': r,ci] rj ri' =
    ( InnerCxt (Inner ci rj) ci rj,
      MapMultEq (Inner ci rj) (r11 ': r) ri',
      SameLengths [ci,rj],
      SameLengths [r11 ': r, ri'])

{- | does H.'H.mXm' and H.'H.mXv'.

vXm and vXv (called dot) might be supported in the future too
-}
multiply :: (H.Product a,
             sh ~ [ _1 ': __1 ,DOne ': __2 ],
             MultiplyCxt sh rj ri')
    => DimMat sh a -> DimMat (rj ': cj) a
    -> DimMat (ri' ': cj) a
multiply (DimMat a) (DimMat b) = DimMat (H.multiply a b)
multiply (DimMat a) (DimVec b) = DimVec (H.mXv a b)

-- | @(Trans sh1 sh2) =>@ asserts that the two matrices are transposes
type family Trans (sh1 :: [[*]]) (sh2 :: [[*]]) :: Constraint
type instance Trans [a11 ': ri, one ': ci]
                    [a11' ': ci', one' ': ri'] =
    (MapMultEq a11 ri' ri, MapMultEq a11 ci ci',
     SameLengths [ci, ci'],
     SameLengths [ri, ri'],
     one ~ DOne, one' ~ DOne,
     a11 ~ a11')

trans :: (Trans sh sh',
          -- need to assert we have valid shapes
          -- (and at least one row)
          sh' ~ [_1 ': __1 , DOne ': _2],
          sh ~ [_3 ': __3, DOne ': _4])
    => DimMat sh a -> DimMat sh' a
trans (DimMat a) = DimMat (H.trans a)

type AreRecips a b = (SameLength' a b, AreRecips' a b)

type family AreRecips' (a :: k) (b :: k) :: Constraint
type instance AreRecips' (a ': as) (b ': bs) = (AreRecips' a b, AreRecips' as bs)
type instance AreRecips' '[] '[] = ()
type instance AreRecips' a b = (a ~ Div DOne b, b ~ Div DOne a)

type SameLength' a b = (SameLength a b, SameLength b a)


-- | if any [k] in the list's length is known, then all other [k] lists in the outer list
-- will be forced to have the same length
type family SameLengths (a :: [[k]]) :: Constraint
type instance SameLengths (a ': b ': bs) = (SameLength' a b, SameLengths (b ': bs))
type instance SameLengths '[b] = ()
type instance SameLengths '[] = ()


inv :: (PInv sh sh', sh' ~ [_1 ': r , DOne ': c], SameLengths [r,c])
    => DimMat sh a -> DimMat sh' a
inv (DimMat a) = DimMat (H.inv a)

{- | type for a pseudoinverse (and inverse):

The single instance comes from looking at inverses from a 2x2 matrix (let's call A):

> a b
> c d

and the inverse * determinant of the original

>  d  -b
> -c   a

In the product A * A^-1 the diagonal is dimensionless ('DOne').

That happens if the row and column type-level unit lists are reciprocals of
eachother ('AreRecips'), so the constraint on the instance of PInv encodes
this exactly (plus some constraints requiring that sh and sh' are at least
1x1)
-}
class PInv (sh :: [[*]]) (sh' :: [[*]]) where
        pinv :: DimMat sh a -> DimMat sh' a
       
instance 
  (sh ~  [_1 ': _2, DOne ': _3],
   sh' ~ [r', c'], c' ~ (DOne ': _4),
   MultiplyCxt sh r' c'inv,
   AreRecips c'inv c') =>
   PInv sh sh' where
    pinv (DimMat a) = DimMat (H.pinv a)

pinvTol :: (PInv sh sh',
#if MIN_VERSION_hmatrix(0,15,0)
-- on hmatrix 13, the pinvTol function has type Double -> Matrix Double -> MatrixDouble, later they generalized to Field t => Double -> Matrix t -> Matrix t
            e ~ Double,
#endif
           sh' ~ [ri2 ': _1 , DOne ': ci2]) => Double -> DimMat sh a -> DimMat sh' a
pinvTol tol (DimMat a) = DimMat (H.pinvTol tol a)

det :: (SameLengths [ri,ci]) => DimMat [ri,ci] a
        -> Quantity (Product ri `Mul` Product ci) a
det (DimMat a) = Dimensional (H.det a)

{- | Numeric.LinearAlgebra.Algorithms.'H.expm'

@y t = expm (scale t a) \`multiply\` y0@ solves the DE @y' = Ay@ where y0 is the
value of y at time 0

-}
expm :: (AreRecips ri ci)
    => DimMat [ri,ci] a
    -> DimMat [ri,ci] a
expm (DimMat a) = DimMat (H.expm a)

{- | Numeric.Container.'H.scale'

-}
scale :: (MapMultEq e ri ri')
    => Quantity e a -> DimMat [ri,ci] a -> DimMat [ri',ci] a
scale (Dimensional t) (DimMat a) = DimMat (H.scale t a)

{- | Numeric.Container.'H.scaleRecip'
-}
scaleRecip :: (MapMultEq e' ri ri', AreRecips e e')
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

{- | @cmap f m@ gives a matrix @m'@

@f@ is applied to 

-}
class CMap f sh sh' e e' where
    cmap :: f -> DimMat sh e -> DimMat sh' e'
instance 
    (ToHLists sh e xs,
     FromHLists sh' e' xs',
     SameLengths [xs,xs'],
     HMapAux (HMap f) xs xs') =>
    CMap f sh sh' e e' where
  cmap f m = fromHLists (HMap f `hMap` (toHLists m :: HList xs) :: HList xs')
    -- maybe there's a way to implement in terms of the real cmap

{- | the slightly involved type here exists because
ci1 and ci2 both start with DOne, but ci2's contribution
to ci3 does not need a DOne at the start. Another way to
read the constraints here is:

> map (*rem) (a11 : ri) = b11 : bi
> ci3 = ci1 ++ map (*rem) ci2

The same idea happens with vconcat.
-}
hconcat :: (MultEq rem a11 b11,
            MapMultEq rem ri1 ri2,
            MapMultEq rem ci2 ci2',
            AppendEq ci1 ci2' ci3) =>
    DimMat [a11 ': ri1,ci1] a -> DimMat [b11 ': ri2, ci2] a -> DimMat [a11 ': ri1, ci3] a
hconcat (DimMat a) (DimMat b) = DimMat (H.fromBlocks [[a, b]])

vconcat :: (AppendEq ri1 ri2 ri3) =>
    DimMat [ri1,ci1] a -> DimMat [ri2,ci1] a -> DimMat [ri3, ci1] a
vconcat (DimMat a) (DimMat b) = DimMat (H.fromBlocks [[a],[b]])

rank (DimMat a) = H.rank a
rows (DimMat a) = H.rows a
cols (DimMat a) = H.cols a

-- | H.'H.rows' except type-level
rowsNT :: DimMat (ri ': ci) a -> Proxy (HLength ri)
rowsNT _ = proxy

-- | H.'H.cols' except type-level
colsNT :: DimMat [ri,ci] a -> Proxy (HLength ci)
colsNT _ = proxy

-- | (m `hasRows` n) constrains the matrix/vector @m@ to have @n@ rows
hasRows :: (HReplicate n DOne, SameLengths [HReplicateR n DOne, ri], -- forwards
            HLength ri ~ n -- backwards
    ) => DimMat (ri ': _1) a -> Proxy (n :: HNat) -> DimMat (ri ': _1) a
hasRows x _ = x

-- | (m `hasRows` n) constrains the matrix/vector @m@ to have @n@ rows
hasCols :: (HReplicate n DOne, SameLengths [HReplicateR n DOne, ci], -- forwards
            HLength ci ~ n -- backwards
    ) => DimMat [ri, ci] a -> Proxy (n :: HNat) -> DimMat [ri,ci] a
hasCols x _ = x

scalar :: (H.Field a,
          sh ~ ['[u], '[DOne]]) => Quantity u a -> DimMat sh a
scalar (Dimensional a) = DimMat (H.scalar a)

{- | Numeric.Container.'H.konst', but the size is determined by the type.

>>> let n = hSucc (hSucc hZero) -- 2
>>> konst ((1::Double) *~ second) `hasRows` n `hasCols` n
2><2 1   1  
s    1.0 1.0
s    1.0 1.0

-}
konst :: forall u us ones a _1.
    (H.Field a,
     HNat2Integral (HLength ones),
     HNat2Integral (HLength us),
     ones ~ (DOne ': _1),
     AllEq DOne _1,
     AllEq u us)
    => Quantity u a -> DimMat [us, ones] a
konst (Dimensional a) = DimMat (H.konst a
    (hNat2Integral (proxy :: Proxy (HLength us)),
     hNat2Integral (proxy :: Proxy (HLength ones))))

-- | identity matrix. The size is determined by the type.
ident :: forall ones a _1.
    (H.Field a, HNat2Integral (HLength ones), ones ~ (DOne ': _1)) =>
    DimMat [ones, ones] a
ident = DimMat (H.ident (hNat2Integral (proxy :: Proxy (HLength ones))))

-- | zero matrix. The size and dimension is determined by the type.
zeroes :: forall r c a _1 _2 _3. (H.Field a,
                              HNat2Integral (HLength r),
                              HNat2Integral (HLength c),
                              c ~ (DOne ': _1),
                              r ~ (_2 ': _3))
    => DimMat [r, c] a
zeroes = DimMat (H.konst 0
        (hNat2Integral (proxy :: Proxy (HLength r)),
         hNat2Integral (proxy :: Proxy (HLength c))))

type family CanAddConst (a :: k) (m :: [[k]]) :: Constraint
type instance CanAddConst a [as, ones] = (AllEq a as, AllEq DOne ones)
type instance CanAddConst a '[as] = (AllEq a as)

type family AllEq (a :: k) (xs :: [k]) :: Constraint
type instance AllEq a (x ': xs) = (a ~ x, AllEq a xs)
type instance AllEq a '[] = ()

addConstant :: (H.Field a, CanAddConst u sh)
    => Quantity u a
    -> DimMat sh a
    -> DimMat sh a
addConstant (Dimensional a) (DimMat b) = DimMat (H.addConstant a b)
addConstant (Dimensional a) (DimVec b) = DimVec (H.addConstant a b)

conj :: DimMat sh a -> DimMat sh a
conj (DimMat a) = DimMat (H.conj a)
conj (DimVec a) = DimVec (H.conj a)

-- | conjugate transpose
ctrans x = conj . trans $ x

diag :: (MapConst DOne v ~ c,
        c ~ (DOne ': _1)
        ) => DimMat '[v] t -> DimMat '[v,c] t
diag (DimVec a) = DimMat (H.diag a)

#if MIN_VERSION_hmatrix(0,15,0)
-- | 'H.blockDiag'. The blocks should be provided as:
--
-- @blockDiag $ 'hBuild' m1 m2 m3@
--
-- only available if hmatrix >= 0.15
diagBlock :: (db ~ DimMat [ri, DOne ': ci] e,
              HMapOut UnDimMat t (H.Matrix e),
              Num e, H.Field e, DiagBlock t db)  => HList t -> db
diagBlock pairs = DimMat (H.diagBlock (hMapOut UnDimMat pairs))
#endif

toHList :: forall e e1 ri result.  (ToHList e e1 ri result)
    => DimMat '[ri] e -> HList result 
toHList (DimVec v) = case hListFromList (H.toList v) :: HList e1 of
    e1 -> hMap AddDimensional e1

fromHList :: forall e ri list.
    (H.Field e,
     HMapOut RmDimensional list e, ToHListRow ri e list)
    => HList list -> DimMat '[ri] e
fromHList xs = DimVec (H.fromList (hMapOut RmDimensional xs))

data RmDimensional = RmDimensional
instance (x ~ Quantity d y) => ApplyAB RmDimensional x y where
        applyAB _ (Dimensional a) = a

class H.Field e => FromHLists sh e xs where
    fromHLists :: HList xs -> DimMat sh e

instance 
        (ToHListRows' ri ci e result,
         HMapOut (HMapOutWith RmDimensional) result [e],
         SameLengths [ri, result],
         (HList resultHead ': _2) ~ result,
         SameLengths [ci, resultHead],
         ci ~ (DOne ': _1), H.Field e,
         sh ~ [ri,ci]) =>
    FromHLists sh e result where
    fromHLists xs = DimMat (H.fromLists (hMapOut (HMapOutWith RmDimensional) xs))

newtype HMapOutWith f = HMapOutWith f
instance (HMapOut f l e, es ~ [e], HList l ~ hl) => ApplyAB (HMapOutWith f) hl es where
    applyAB (HMapOutWith f) = hMapOut f

class ToHLists sh e xs where
    toHLists :: DimMat sh e -> HList xs
instance (ToHListsCxt e e1 e2 ri ci xs) => ToHLists [ri,ci] e xs where
  toHLists (DimMat m) = case hListFromList (map hListFromList (H.toLists m) :: [HList e1]) :: HList e2 of
    e2 -> hMap (HMap AddDimensional) e2

type ToHList e e1 ri result =
    (HListFromList e e1,
     SameLengths [e1,result,ri],
     HMapAux AddDimensional e1 result,
     ToHListRow ri e result)

type family ToHListRow (a :: [*]) e (b :: [*]) :: Constraint
type instance ToHListRow (a ': as) e (b ': bs) = (Quantity a e ~ b, ToHListRow as e bs)

-- | performance (compile-time) is pretty bad
type ToHListsCxt e e1 e2 ri ci result =
    (HListFromList e e1,
     HListFromList (HList e1) e2,
     SameLengths [ci,e1], SameLengths [e2,result,ri],
     HMapAux (HMap AddDimensional) e2 result,
     ToHListRows' ri ci e result)

class HListFromList e e' where
        hListFromList :: [e] -> HList e'
instance HListFromList e '[] where
        hListFromList _ = HNil
instance (e ~ e', HListFromList e es) => HListFromList e (e' ': es) where
        hListFromList (e : es) = e `HCons` hListFromList es 

class ToHListRows' (ri :: [*]) (ci :: [*]) (e :: *) (rows :: [*])
instance ToHListRows' '[] ci e '[]
instance (ToHListRows' ri ci e rows,
          MapMultEq r ci ci',
          HMapCxt (AddQty e) (HList ci') hListRow ci' row')
  => ToHListRows' (r ': ri) ci e (hListRow ': rows)

data AddQty e
instance (qty ~ Quantity d e) => ApplyAB (AddQty e) d qty

data AddDimensional = AddDimensional
instance (Quantity t x ~ y) => ApplyAB AddDimensional x y where
        applyAB _ x = Dimensional x

data UnDimMat = UnDimMat
instance (DimMat sh a ~ x, H.Matrix a ~ y) => ApplyAB UnDimMat x y where
        applyAB _ (DimMat x) = x

class DiagBlock (bs :: [*]) t
instance (DiagBlock as as', AppendShOf a as' ~ t) => DiagBlock (a ': as) t 
instance (a ~ a') => DiagBlock '[a] a'

type family Append (a :: [k]) (b :: [k]) :: [k]
type instance Append (a ': as) b = a ': Append as b
type instance Append '[] b = b

type family AppendEq' (a :: [k]) (b :: [k]) (ab :: [k]) :: Constraint
type instance AppendEq' (a ': as) b (a' ': abs) = (a ~ a', AppendEq' as b abs)
type instance AppendEq' '[] b abs = (b ~ abs)

-- | a bit overkill?
--  @a ++ b = ab@
type AppendEq a b ab =
   (ab ~ Append a b,
    AppendEq' a b ab,
    SameLengths [DropPrefix a ab, b],
    SameLengths [DropPrefix b ab, a])

type family DropPrefix (a :: [k]) (ab :: [k2]) :: [k2]
type instance DropPrefix (a ': as) (a' ': abs) = DropPrefix as abs
type instance DropPrefix '[] bs = bs

-- | rework to follow AppendEq?
type family AppendDims (sh :: [[*]]) (sh' :: [[*]]) :: [[*]]
type instance AppendDims [a,b] [c,d] = [Append a c, Append b d]
type instance AppendDims '[a] '[b] = '[Append a b]

type family AppendShOf (a :: *) (b :: *) :: *
type instance AppendShOf (DimMat sh t) (DimMat sh' t) = DimMat (AppendDims sh sh') t

class PairsToList a t where
        pairsToList :: a -> [H.Matrix t]
instance PairsToList () t where
        pairsToList _ = []
instance (PairsToList b t, t' ~ t) => PairsToList (DimMat sh t',b) t where
        pairsToList (DimMat a,b) = a : pairsToList b

type family EigCxt (sh :: [[*]]) (eigenValue  :: [*]) (eigenVector  :: k) :: Constraint

type instance EigCxt [r,c] eigval [cinv,d] =
  ( MapConst DOne r ~ d,
    SameLengths [r,c,cinv,d,eigval],
    AreRecips c cinv,
    ZipWithMul r c eigval)

-- | when no eigenvectors are needed
type instance EigCxt [r,c] eigval '() =
  ( SameLengths [r,c,eigval], ZipWithMul r c eigval)

{- $eigs

The Hmatrix eig factors A into P and D where A = P D inv(P) and D is diagonal.

The units for eigenvalues can be figured out:

>               _____
>      -1       |  c
> P D P  = A =  |r
>               |

>       _______
>       |   d
> P   = |c
>       |

>       _______
>       |   -1
>       |  c
>  -1   |   
> P   = | -1
>       |d

So we can see that the dimension labeled `d-1` in P inverse is actually the
same `c` in `A`. The actual units of `d` don't seem to matter because the
`inv(d)` un-does any units that the `d` adds. So `d` can be all DOne. But
another choice, such as 1/c would be more appropriate, since then you can
expm your eigenvectors (not that that seems to be something people do)?

To get the row-units of A to match up, sometimes `D` will have units. 
The equation ends up as D/c = r

Please ignore the type signatures on 'eig' 'eigC' etc. instead look at the type of
'wrapEig' 'wrapEigOnly' together with the hmatrix documentation (linked).

Perhaps the convenience definitions `eig m = wrapEig H.eig m` should be in
another module.
-}

-- | 'wrapEig' H.'H.eig'
eig m = wrapEig H.eig m
-- | 'wrapEig' H.'H.eigC'
eigC m = wrapEig H.eigC m
-- | 'wrapEig' H.'H.eigH'
eigH m = wrapEig H.eigH m
-- | 'wrapEig' H.'H.eigH''
eigH' m = wrapEig H.eigH' m
-- | 'wrapEig' H.'H.eigR'
eigR m = wrapEig H.eigR m
-- | 'wrapEig' H.'H.eigS'
eigS m = wrapEig H.eigS m
-- | 'wrapEig' H.'H.eigS''
eigS' m = wrapEig H.eigS' m
-- | 'wrapEig' H.'H.eigSH'
eigSH m = wrapEig H.eigSH m
-- | 'wrapEig' H.'H.eigSH''
eigSH' m = wrapEig H.eigSH' m

-- | 'wrapEigOnly' H.'H.eigOnlyC'
eigOnlyC m = wrapEigOnly H.eigOnlyC m
-- | 'wrapEigOnly' H.'H.eigOnlyH'
eigOnlyH m = wrapEigOnly H.eigOnlyH m
-- | 'wrapEigOnly' H.'H.eigOnlyR'
eigOnlyR m = wrapEigOnly H.eigOnlyR m
-- | 'wrapEigOnly' H.'H.eigOnlyS'
eigOnlyS m = wrapEigOnly H.eigOnlyS m
-- | 'wrapEigOnly' H.'H.eigenvalues'
eigenvalues m = wrapEigOnly H.eigenvalues m
-- | 'wrapEigOnly' H.'H.eigenvaluesSH'
eigenvaluesSH m = wrapEigOnly H.eigenvaluesSH m
-- | 'wrapEigOnly' H.'H.eigenvaluesSH''
eigenvaluesSH' m = wrapEigOnly H.eigenvaluesSH' m

wrapEig :: (ones ~ (DOne ': _1), EigCxt [r,c] eigVal [cinv,ones],
    H.Field y, H.Field z)
    => (H.Matrix x -> (H.Vector y, H.Matrix z)) ->
    DimMat [r,c] x ->
    (DimMat '[eigVal] y, DimMat [cinv,ones] z)
wrapEig hmatrixFun (DimMat a) = case hmatrixFun a of
    (e,v) -> (DimVec e, DimMat v)

wrapEigOnly :: (EigCxt [r,c] eigVal '(), H.Field y)
    => (H.Matrix x -> H.Vector y) ->
    DimMat [r,c] x -> DimMat '[eigVal] y
wrapEigOnly hmatrixFun (DimMat a) = case hmatrixFun a of
    (e) -> DimVec e

