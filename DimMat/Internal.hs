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

  Friendly syntax for introduction (more matD)

  Friendly syntax for elimination (matD as pattern)

  A pretty-printer that multiplies out the dimensions at each place?
  The current show instance

  A clean way to get the n x n dimensionless identity matrix

  check that all types that could could be inferred are

  default columns/rows to dimensionless?

  missing operations (check export list comments)
-}
module DimMat.Internal (
   -- * Data.Packed.Vector
   -- * Data.Packed.Matrix
   -- ** dimension
   cols, rows,
   colsNT, rowsNT,
   -- (><),
   trans,
   -- reshape, flatten, fromLists, toLists, buildMatrix,
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
   cmap,
   konst,
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
   -- pinv, pinvtol,
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
   AtEq, MapMul, Inner, InvCxt, SameLengths, Product, MapRecip, ScaleCxt,
   ZipWithZipWithMul, MapMapConst, Len, CanAddConst, PPUnits,
   PPUnits', Head, MapDiv, AreRecips, ZipWithMul, PairsToList,
   DiagBlock, MapConst, SameLength', AppendShOf,
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
instance (Show a, PPUnits sh) => Show (DimMat sh a) where
    showsPrec _ (DimVec v) = case ppUnits (Proxy :: Proxy sh) of
        [rs] ->
            displayS $
            renderPretty 0.1 80 $ vcat
             [ dullgreen (string label) <+> string (show e)
                | (e,label) <- H.toList v `zip` pad rs ]
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

-- | just for types produced by the matD quasiquote
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

type family MapMultEq (a :: k) (xs :: [k]) (ys :: [k]) :: Constraint
type instance MapMultEq a (x ': xs) (y ': ys) = (MultEq a x y, MapMultEq a xs ys)
type instance MapMultEq a '[] '[] = ()

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

-- | 'length', encoded with 'N.S' and 'N.Z'
type family Len (a :: [*]) :: *
type instance Len '[] = N.Z
type instance Len (a ': as) = N.S (Len as)

-- | @At a n@ is the type-level version of @a !! n@
type family At (a :: [k]) n :: k
type instance At (a ': as) N.Z = a
type instance At (a ': as) (N.S n) = At as n

-- | @AtEq a n b m c@ calculates @(At a n \`Mult\` At b m) ~ c@,
-- but also can infer part of the `a` if the `b` and `c` are known
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

-- | Data.Packed.Matrix.'H.@@>'
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
scale :: (ScaleCxt e ri ri')
    => Quantity e a -> DimMat [ri,ci] a -> DimMat [ri',ci] a
scale (Dimensional t) (DimMat a) = DimMat (H.scale t a)

{- | Numeric.Container.'H.scaleRecip'
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
-- ideally any transformation of units (say replace all metre by seconds)
-- that still fits should work. This will involve re-doing the lookups
-- done in 'matD'
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

rowsNT :: DimMat [ri,ci] a -> Proxy (Len ri)
rowsNT _ = Proxy

colsNT :: DimMat [ri,ci] a -> Proxy (Len ci)
colsNT _ = Proxy

scalar :: (H.Field a,
          sh ~ ['[u], '[DOne]]) => Quantity u a -> DimMat sh a
scalar (Dimensional a) = DimMat (H.scalar a)

konst :: forall sh u us ones a _1. (H.Field a,
                              N.NumType (Len ones),
                              N.NumType (Len us),
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
conj (DimVec a) = DimVec (H.conj a)

-- | conjugate transpose
ctrans :: (one ~ DOne,
         sh  ~ [a11 ': ri, one ': ci],
         sh' ~ [a11 ': ci, one ': ri])
         => DimMat sh a -> DimMat sh' a
ctrans = conj . trans

diag :: (MapConst DOne v ~ c,
        c ~ (DOne ': _1)
        ) => DimMat '[v] t -> DimMat '[v,c] t
diag (DimVec a) = DimMat (H.diag a)

#if MIN_VERSION_hmatrix(0,15,0)
-- | 'H.blockDiag'. The blocks should be provided as:
--
-- > blockDiag (m1, (m2, (m3, ())))
--
-- XXX should we bring in HList for this?
--
-- only available if hmatrix >= 0.15
diagBlock :: (PairsToList t e, db ~ DimMat [ri, DOne ': ci] e,
              Num e, H.Field e, DiagBlock t db)  => t -> db
diagBlock pairs = DimMat (H.diagBlock (pairsToList pairs))
#endif

class DiagBlock (bs :: *) t
instance (DiagBlock as as', AppendShOf a as' ~ t) => DiagBlock (a,as) t 
instance (a ~ a') => DiagBlock (a,()) a'

type family Append (a :: [k]) (b :: [k]) :: [k]
type instance Append (a ': as) b = a ': Append as b
type instance Append '[] b = b

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

