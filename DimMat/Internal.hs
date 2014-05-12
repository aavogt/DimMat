{-# OPTIONS_GHC -fno-warn-dodgy-exports #-}
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
   -- broken
   -- ToHLists(toHLists), toHList, FromHLists(fromHLists), fromHList,
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
   recipMat,
   addConstant,
   add,
   sub,
   mulMat, mulMat,
   divideMat, divideVec,
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
   DimMat.Internal.dot,
   -- absSum, norm1, norm2, normInf,
   -- norm1, normInf,
   pnorm,
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
   {-
   wrapEig, wrapEigOnly,
   EigV, EigE,
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
   -}

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
#ifdef WITH_Ad
   diff,
#endif
   -- * actually internal
   module DimMat.Internal,
  ) where
import Foreign.Storable (Storable)      
import GHC.Exts (Constraint)
import Prelude


import Data.Type.Equality (type (==))
import qualified Prelude as P
import Prelude (Double)
import Numeric.Units.Dimensional hiding (Mul, Div)
import Numeric.Units.Dimensional.Prelude hiding (Mul, Div)
import qualified Numeric.Units.Dimensional as D

import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.LAPACK as H

import Text.PrettyPrint.ANSI.Leijen
import Data.List (transpose)

import Data.HList.CommonMain hiding (MapFst)


-- | a version of Numeric.Units.Dimensional.'D.Mul' which
-- requires the arguments to include the 'D.Dim' type constructor
class (D.Mul a b c) => Mul a b c

instance (D.Mul a b c,
       a ~ Dim l m t i th n j,
       b ~ Dim l' m' t' i' th' n' j',
       c ~ Dim l'' m'' t'' i'' th'' n'' j'') =>
       Mul a b c

-- | a version of Numeric.Units.Dimensional.'D.Div' which
-- requires the arguments to include the 'D.Dim' type constructor
class (D.Div a b c) => Div a b c

instance (D.Div a b c,
       a ~ Dim l m t i th n j,
       b ~ Dim l' m' t' i' th' n' j',
       c ~ Dim l'' m'' t'' i'' th'' n'' j'') =>
       Div a b c

#ifdef WITH_AD
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
#endif


{- $ad

TODO: gradients, hessians, etc.

Types for derivative towers can see hlist's @HList\/Data\/HList\/broken\/Lazy.hs@,
but laziness doesn't really make much sense if the @take@ that is eventually used
to get a finite list for printing etc.

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


{- | Generalization of 'Dimensional' to matrices and vectors. Units
in each coordinate are known at compile-time. This wraps up HMatrix.
-}
data D (sh :: ( *, [[ * ]])) e where
  DMat :: (H.Container H.Matrix e, H.Field e)
    => H.Matrix e -> D '(r1,[r, c]) e
      -- ^ the units at coordinate i,j are
      -- @(r1 ': r)_i (DOne ': c)_j@
  DVec :: (H.Container H.Vector e, H.Field e)
    => H.Vector e -> D '(r1, '[r]) e
      -- ^ the units at coordinate i are
      -- @(r1 ': r)_i@
  DScal :: (H.Field e) => e -> D '(r1,'[]) e
      -- ^ the same as Dimensional

-- very crude
instance (Show a, PPUnits sh) => Pretty (D sh a) where
    pretty (DVec v) = case ppUnits (Proxy :: Proxy sh) of
        [rs] -> vcat
             [ dullgreen (string label) <+> string (show e)
                | (e,label) <- H.toList v `zip` pad rs ]
    pretty (DMat m) = case ppUnits (Proxy :: Proxy sh) of
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

instance Pretty (D sh a) => Show (D sh a) where
    showsPrec p x = showsPrec p (pretty x)

pad :: [String] -> [String]
pad [] = []
pad xs = let
    w = maximum (map length xs)
    in map (\x -> take w $ x ++ replicate w ' ') xs


class PPUnits (sh :: k) where
    ppUnits :: Proxy sh -> [[String]]

instance forall (r1 :: *) (r::[*]) (c :: [*]) l m t i th n j.
      (PPUnits [r1 ': r, c], Show (Quantity r1 Int), PPUnits' c, PPUnits' r,
       r1 ~ Dim l m t i th n j) => PPUnits '(r1, [r,c]) where
    ppUnits _ = ppUnits (Proxy :: Proxy [r1 ': r, DOne ': c])

instance (PPUnits' x, PPUnits xs) => PPUnits (x ': xs) where
    ppUnits _ = ppUnits' (Proxy :: Proxy x) : ppUnits (Proxy :: Proxy xs)
instance PPUnits '[] where
    ppUnits _ = []

class PPUnits' (sh :: [ * ]) where
    ppUnits' :: Proxy sh -> [String]
instance (PPUnits' xs) => PPUnits' (DOne ': xs) where
    ppUnits' _ = "1" : ppUnits' (Proxy :: Proxy xs)
instance (ShowDimSpec x, PPUnits' xs) => PPUnits' (x ': xs) where
    ppUnits' _ = showDimSpec (Proxy :: Proxy x) : ppUnits' (Proxy :: Proxy xs)
instance PPUnits' '[] where
    ppUnits' _ = []

class ShowDimSpec a where
    showDimSpec :: Proxy a -> String

instance (Show (Quantity a Int), Dim l m t i th n j ~ a) => ShowDimSpec a where
    showDimSpec _ = case drop 2 $ show (1 *~ (Dimensional 1 :: Unit a Int)) of
          "" -> "1"
          x -> x


-- | @a*b = c@ when any are lists
class MultEq (a :: k1) (b :: k2) (c :: k3)

-- instance (Zip3 MultEq aas bbs ccs) => MultEq aas bbs ccs
instance Zip3 Mul as bs cs => MultEq as bs cs
instance Zip1 Mul as b  c  => MultEq as b  c
instance Zip1 Mul bs a  c  => MultEq a  bs c
instance (Zip1 Mul cs aInv b,
          Mul a aInv DOne) => MultEq a  b  cs

instance (SameLength as bs,
        Zip2 Mul as bs c)  => MultEq as bs c
instance (SameLength as cs,
          Zip2 Div as cs bInv,
          Mul b bInv DOne) => MultEq as b  cs
instance (Zip2 Div bs cs aInv,
          SameLength bs cs,
          Mul a aInv DOne) => MultEq a  bs cs
instance Mul a b c => MultEq a b c


class (SameLength a b, SameLength b c) =>
    Zip3
      (op :: k -> k -> k -> Constraint)
      (a :: [k])
      (b :: [k])
      (c :: [k])

instance (SameLength aas bbs,
          SameLength ccs bbs,
          op a b c,
          (a ': as) ~ aas,
          (b ': bs) ~ bbs,
          (c ': cs) ~ ccs,
          Zip3 op as bs cs) => Zip3 op aas bbs ccs
instance Zip3 op '[] '[] '[]


class (SameLength a b) =>
    Zip2
      (op :: k -> k -> k -> Constraint)
      (a :: [k])
      (b :: [k])
      (c ::  k)

instance (SameLength aas bbs,
      op a b c,
      (a ': as) ~ aas,
      (b ': bs) ~ bbs,
      Zip2 op as bs c) => Zip2 op aas bbs c

instance Zip2 op '[] '[] c

class Zip1
      (op :: k -> k -> k -> Constraint)
      (a :: [k])
      (b ::  k)
      (c ::  k)

instance ((a ': as) ~ aas,
    op a b c,
    Zip1 op as b c) => Zip1 op aas b c

instance Zip1 op '[] b c


type family HeadOf (x :: k) (xs :: [k]) :: Constraint
type instance HeadOf x (x' ': u1) = (x ~ x')
type instance HeadOf x '[] = ("DimMat.HeadOf" ~ "given empty list")

{- | given @ijs :: [[Quantity a]]@ (except the : and [] constructors are
actually (,) and (), ie. a HList that doesn't use the HList constructors),
calculate a @DimMat rowUnits colUnits@, where the outer product of rowUnits
and colUnits gives the units at each index in the ijs.  The first element
of colUnits is DOne.
-}


class (SameLength a ab) => Outer a b ab
instance Outer '[] b '[]
instance (SameLength aas ccs,
          (a ': as) ~ aas,
          (c ': cs) ~ ccs,
          MultEq a b c,
          Outer as b cs) 
  => Outer aas b ccs

class DimMatFromTuple ijs r1 r c e


type family TupleToHListU (a :: *) :: [*]
type instance TupleToHListU (a, b) = () ': TupleToHListU b
type instance TupleToHListU () = '[]

type family TuplesToHListU (a :: *) :: [[*]]
type instance TuplesToHListU (a, b) = TupleToHListU a ': TuplesToHListU b 
type instance TuplesToHListU () = '[]

instance (Outer (r1 ': r) (DOne ': c) ijs',
      DMFromTuple1 e ijs ijs',
      SameLength (TuplesToHListU ijs) ijs') => DimMatFromTuple ijs r1 r c e

-- | helper for 'DimMatFromTuple'
type family DMFromTuple1 e b (b' :: [[*]]) :: Constraint
type family DMFromTuple2 e b (b' :: [*]) :: Constraint
type family DMFromTuple3 e b b' :: Constraint

type instance DMFromTuple3 e (Quantity b e') b' = (e ~ e', b ~ b')
type instance DMFromTuple1 e (x, xs) (x' ': xs') = (TupleToHListU x `SameLength` x',
                                                    DMFromTuple2 e x x', DMFromTuple1 e xs xs')
type instance DMFromTuple1 e () xs = (xs ~ '[])
type instance DMFromTuple2 e (x, xs) (x' ': xs') = (DMFromTuple3 e x x', DMFromTuple2 e xs xs')
type instance DMFromTuple2 e () xs = (xs ~ '[])

-- | just for types produced by the matD quasiquote
toDM :: DimMatFromTuple ijs r1 r c e => ijs -> Proxy (D '(r1, [r, c]) e)
toDM _ = Proxy


-- | @InnerCxt t a b = t ~ 'H.dot' a b@
type family InnerCxt (t :: k) (a :: [k]) (b :: [k]) :: Constraint
type instance InnerCxt t (a ': as) (b ': bs) = (MultEq a b t, InnerCxt t as bs)
type instance InnerCxt t '[] '[] = ()

class (SameLength a b, InnerCxt c a b) => Inner (a :: [*]) (b :: [*]) (c :: *)

instance (SameLength aas bbs, InnerCxt c aas bbs) => Inner aas bbs c

-- | @ProdEq a b@ is @product a ~ b@
class ProdEq a b
instance (ProdEq as b', Mul a b' b) => ProdEq (a ': as) b
instance (dOne ~ DOne) => ProdEq '[] dOne

-- | @RecipEq a aInv@ is @a*aInv ~ DOne@ (or a list of DOne)
class RecipEq (a :: k) (aInv :: k)
instance (MultEq as aInvs DOne) => RecipEq as aInvs


-- | @AtEq a n b m c@ calculates @(At a n \`Mult\` At b m) ~ c@,
-- but also can infer part of the `a` if the `b` and `c` are known
type family AtEq2 (a :: [k]) (n :: HNat) (b :: [k]) (m :: HNat) (c :: k) :: Constraint
type instance AtEq2  (a ': as) HZero (b ': bs) HZero c = (MultEq a b c)
type instance AtEq2  (a ': as) (HSucc n) bs m c = AtEq2 as n bs m c
type instance AtEq2  as HZero (b ': bs) (HSucc m) c = AtEq2 as HZero bs m c

type family AtEq (a :: [k]) (n :: HNat) (b :: k) :: Constraint
type instance AtEq (a ': as) HZero b = (a ~ b)
type instance AtEq (a ': as) (HSucc n) b = AtEq as n b

-- | Data.Packed.Vector.'H.@>'
(@>) :: (HNat2Integral i, AtEq r i ri, MultEq r1 ri u)
    => D '(r1,'[r]) a
    -> Proxy i
    -> Quantity u a
DVec v @> i = Dimensional (v H.@> hNat2Integral i)

-- | Data.Packed.Matrix.'H.@@>'
(@@>) :: (HNat2Integral i, HNat2Integral j, AtEq2 (r1 ': r) i (DOne ': c) j ty)
    => D '(r1, [r,c]) a
    -> (Proxy i, Proxy j)
    -> Quantity ty a
DMat m @@> (i,j) = Dimensional (m H.@@> (hNat2Integral i,hNat2Integral j))

pnorm :: (AllEq r1 r, AllEq DOne c)
         => H.NormType -> D '(r1, [r, c]) a -> Quantity r1 (H.RealOf a)
pnorm normType (DMat a) = Dimensional (H.pnorm normType a)

{- | @AllEq a xs@ is like @all (all (a==)) xs@, @all (a ==) xs@, @a == xs@:
whichever amount of [ ] is peeled off before making the comparison (with ~)
-}
class AllEq (a :: k1) (xs :: k2)

instance (a ~ x, AllEq a xs) => AllEq a (x ': xs)
instance AllEq a '[]
instance AllEq '[] xs
instance (AllEq a xs, AllEq as xs) => AllEq (a ': as) xs



{- | @c = a `dot` b@ is one of:

> c_ij = sum_j a_ij b_jk
> c_k  = sum_j a_j  b_jk
> c_i  = sum_j a_ij b_j
> c    = sum_j a_j  b_j

-}
class Dot a b c where
    dot :: H.Element e => D a e -> D b e -> D c e

instance
    ( MultEq (a ': ra) b (c ': rc),
      MultEq ca rb b,
      shA ~ '(a,[ra, ca]),
      shB ~ '(b, rb ': cb),
      shC ~ '(c, rc ': cb))
    => Dot shA shB shC where
    dot (DMat a) (DMat b) = DMat (H.multiply a b)
    dot (DMat a) (DVec b) = DVec (H.mXv a b)
    {-
    dot (DVec a) (DMat b) = DVec (H.vXm a b)
    dot (DVec a) (DVec b) = DScal (H.dot a b)
    -}


class Trans a b where
    trans :: D a e -> D b e

instance (a ~ '(r1, [x,y]),
          b ~ '(r1, [y,x]))
    => Trans a b where
    trans (DMat a) = DMat (H.trans a)

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
class PInv a b where
  pinv :: D a e -> D b e

type family LengthSndTwo (a :: k) :: Constraint
type instance LengthSndTwo '(a, as) = SameLength as '[(), ()]


type AreRecips a b = MultEq a b DOne
instance (MultEq ra cb a,
          MultEq ca rb b,
          AreRecips a b,
          bSh ~ '(b, [rb, cb]),
          aSh ~ '(a, [ra, ca]))
    => PInv aSh bSh where
  pinv (DMat a) = DMat (H.pinv a)

inv :: (PInv a b,
        b ~ '(t1, [t2, t3]))
  => D a e -> D b e
inv (DMat a) = DMat (H.inv a)

       

pinvTol :: (PInv a b,
#if !MIN_VERSION_hmatrix(0,15,0)
-- on hmatrix 13, the pinvTol function has type Double -> Matrix Double -> MatrixDouble, later they generalized to Field t => Double -> Matrix t -> Matrix t
            a ~ Double,
#endif
           b ~ '(t1, [t2, t3]) )
  => Double -> D a e -> D b e
pinvTol tol (DMat a) = DMat (H.pinvTol tol a)


class Det a b where
    det :: D a e -> Quantity b e

instance (SameLength r c,
       ProdEq (r1 ': r) (pr :: *),
       ProdEq c (pc :: *),
       MultEq pr pc b,
       a ~ '(r1, [r, c])) =>
    Det a b where
  det (DMat a) = Dimensional (H.det a)

{- | Numeric.LinearAlgebra.Algorithms.'H.expm'

@y t = expm (scale t a) \`multiply\` y0@ solves the DE @y' = Ay@ where y0 is the
value of y at time 0

-}
expm :: (AreRecips r c)
    => D '(r1, [r, c]) a
    -> D '(r1, [r, c])  a
expm (DMat a) = DMat (H.expm a)

{- | Numeric.Container.'H.scale'

-}
class Scale a b c where
  scale :: Quantity a e -> D b e -> D c e

fromQty :: H.Field e => Quantity a e -> D '(a, '[]) e
fromQty (Dimensional a) = DScal a

toQty :: D '(a, '[]) e -> Quantity a e
toQty (DScal a) = Dimensional a

instance (MultEq a (r1 ': r) (r1' ': r'),
      b ~ '(r1, r ': rs),
      c ~ '(r1', r' ': rs)) =>
  Scale a b c where
  scale (Dimensional t) (DMat a) = DMat (H.scale t a)
  scale (Dimensional t) (DVec a) = DVec (H.scale t a)

{- | Numeric.Container.'H.scaleRecip'
-}
class ScaleRecip a b c where
  scaleRecip :: D '(a, '[]) e -> D b e -> D c e

class ScaleRecip1 (bool :: Bool) a b c where
  scaleRecip1 :: Proxy bool -> D '(a, '[]) e -> D b e -> D c e

instance
  (ScaleRecipCxt r1 r1' r r' rs rs' a b c
   , rs' ~ '[ t1 ]
   ) => ScaleRecip1 True a b c where
  scaleRecip1 _ (DScal t) (DMat a) = DMat (H.scaleRecip t a)

instance
  (ScaleRecipCxt r1 r1' r r' rs rs' a b c,
   rs' ~ '[]) =>
  ScaleRecip1 False a b c where
  scaleRecip1 _ (DScal t) (DVec a) = DVec (H.scaleRecip t a)

instance (bool1 ~ (HLength bs == HSucc (HSucc HZero)),
          bool2 ~ (HLength cs == HSucc (HSucc HZero)),
          bool1 ~ bool2,
    ScaleRecip1 bool1 a '(b, bs) '(c, cs)) => ScaleRecip a '(b, bs) '(c, cs) where
  scaleRecip = scaleRecip1 (Proxy :: Proxy bool1)


type ScaleRecipCxt (r1 :: *) (r1' :: *) r r' rs rs' (a :: *) b c =
  (MultEq a (r1' ': r') (r1 ': r) ,
   MultEq rs rs' DOne,
   b ~ '(r1, r ': rs),
   c ~ '(r1', r' ': rs'))


-- | a shortcut for @scaleRecip (DScal 1)@
recipMat :: forall b c e. (H.Field e, ScaleRecip DOne b c) => D b e -> D c e
recipMat m = scaleRecip (DScal 1 :: D '(DOne, '[]) e) m


liftH2 :: 
  ( forall h f. (H.Container f e, h ~ f e) => h -> h -> h) ->
    D a e -> D a e -> D a e
liftH2 f (DMat a) (DMat b) = DMat (f a b)
liftH2 f (DVec a) (DVec b) = DVec (f a b)

add a b = liftH2 H.add a b
sub a b = liftH2 H.sub a b

mulMat :: ( MultEq as bs cs, MultEq a b c, cs ~ '[t1 , t2] )
  => D '(a,as) e -> D '(b,bs) e -> D '(c,cs) e
mulMat (DMat a) (DMat b) = DMat (H.mul a b)

mulVec :: ( MultEq as bs cs, MultEq a b c, cs ~ '[t1] )
  => D '(a,as) e -> D '(b,bs) e -> D '(c,cs) e
mulVec (DVec a) (DVec b) = DVec (H.mul a b)

divideMat :: ( MultEq as cs bs, MultEq a c b, cs ~ '[t1 , t2] )
  => D '(a,as) e -> D '(b,bs) e -> D '(c,cs) e
divideMat (DMat a) (DMat b) = DMat (H.divide a b)

divideVec :: ( MultEq as cs bs, MultEq a c b, cs ~ '[t1] )
  => D '(a,as) e -> D '(b,bs) e -> D '(c,cs) e
divideVec (DVec a) (DVec b) = DVec (H.divide a b)

arctan2 :: (bs ~ MapMapConst DOne as) => D '(a,as) e -> D '(a,as) e -> D '(b,bs) e
arctan2 (DMat a) (DMat b) = DMat (H.arctan2 a b)
arctan2 (DVec a) (DVec b) = DVec (H.arctan2 a b)

equal :: D a e -> D a e -> Bool
equal (DMat a) (DMat b) = H.equal a b
equal (DVec a) (DVec b) = H.equal a b

{- | @cmap f m@ gives a matrix @m'@

@f@ is applied to 

-}
class CMap f a b e e' where
    cmap :: f -> D a e -> D b e'
    {-
instance 
    (ToHLists sh e xs,
     FromHLists sh' e' xs',
     SameLengths [xs,xs'],
     HMapAux (HMap f) xs xs') =>
    CMap f sh sh' e e' where
  cmap f m = fromHLists (HMap f `hMap` (toHLists m :: HList xs) :: HList xs')
    -- maybe there's a way to implement in terms of the real cmap
    -- -}

type family AppendEq' (a :: [k]) (b :: [k]) (ab :: [k]) :: Constraint
type instance AppendEq' (a ': as) b (a' ': abs) = (a ~ a', AppendEq' as b abs)
type instance AppendEq' '[] b abs = (b ~ abs)

-- | a bit overkill?
--  @a ++ b = ab@
type AppendEq a b ab =
   (ab ~ HAppendR a b,
    AppendEq' a b ab,
    SameLength (DropPrefix a ab) b,
    SameLength (DropPrefix b ab) a)


type instance HAppendR (x ': xs) ys = x ': HAppendR xs ys
type instance HAppendR '[] ys = ys


type family DropPrefix (a :: [k]) (ab :: [k]) :: [k]
type instance DropPrefix (a ': as) (a' ': abs) = DropPrefix as abs
type instance DropPrefix '[] bs = bs

{- | the slightly involved type here exists because
ci1 and ci2 both start with DOne, but ci2's contribution
to ci3 does not need a DOne at the start. Another way to
read the constraints here is:

> map (*rem) (a11 : ri) = b11 : bi
> ci3 = ci1 ++ map (*rem) ci2

The same idea happens with vconcat.
-}
hconcat ::
    ( MultEq (rem :: *) a b,
      MultEq rem ra rb,
      MultEq rem (DOne ': cb) cb',
      AppendEq ca cb' cc ) =>
    D '(a, [ra,ca]) e -> D '(b, [rb, cb]) e -> D '(a, [ra, cc]) e
hconcat (DMat a) (DMat b) = DMat (H.fromBlocks [[a, b]])

vconcat :: (AppendEq ra (b ': rb) rc) =>
    D '(a, '[ra,ca]) e -> D '(b, '[rb,ca]) e -> D '(a, '[rc,ca]) e
vconcat (DMat a) (DMat b) = DMat (H.fromBlocks [[a],[b]])

rank, rows, cols :: D t a -> Int 
rank (DMat a) = H.rank a
rows (DMat a) = H.rows a
cols (DMat a) = H.cols a

-- | H.'H.rows' except type-level
rowsNT :: D '(a, r ': c) e -> Proxy (HLength (a ': ri))
rowsNT _ = Proxy

-- | H.'H.cols' except type-level
colsNT :: D '(a, r ': c ': cs) e -> Proxy (HLength (DOne ': c))
colsNT _ = Proxy

-- | (m `hasRows` n) constrains the matrix/vector @m@ to have @n@ rows
hasRows :: (SameLength (HReplicateR n ()) r, -- forwards
            HLength r ~ n -- backwards
    ) => D '(a, ra ': ca) e -> Proxy (n :: HNat) -> D '(a, ra ': ca) e
hasRows x _ = x

-- | (m `hasRows` n) constrains the matrix/vector @m@ to have @n@ rows
hasCols :: (SameLength (HReplicateR n ()) ci, -- forwards
            HLength ci ~ n -- backwards
    ) => D '(a, ra ': ca ': rest) e -> Proxy (n :: HNat) -> D '(a, ra ': ca ': rest) e
hasCols x _ = x

-- | H.'H.scalar'
class (MapConst '[] as ~ as) => Scalar as where
    scalar :: D '(a, '[]) e -> D '(a, as) e

instance Scalar '[ '[] ] where
  scalar (DScal a) = DVec (H.scalar a)

instance Scalar '[ '[], '[] ] where
  scalar (DScal a) = DMat (H.scalar a)

{- | Numeric.Container.'H.konst', but the size is determined by the type.

>>> let n = hSucc (hSucc hZero) -- 2
>>> konst ((1::Double) *~ second) `hasRows` n `hasCols` n
2><2 1   1  
s    1.0 1.0
s    1.0 1.0

-}
konst :: forall e a ra ca.
    (H.Field e,
     HNat2Integral (HLength (a ': ra)),
     HNat2Integral (HLength (DOne ': ca)),
     AllEq DOne ca,
     AllEq a ra)
    => D '(a, '[]) e -> D '(a, '[ra, ca]) e
konst (DScal a) = DMat (H.konst a
    (hNat2Integral (Proxy :: Proxy (HLength (a ': ra))),
     hNat2Integral (Proxy :: Proxy (HLength (DOne ': ca)))))


-- | identity matrix. The size is determined by the type.
ident :: forall ones e.
    (H.Field e, HNat2Integral (HLength (DOne ': ones))) =>
    D '(DOne, [ones, ones]) e
ident = DMat (H.ident (hNat2Integral (Proxy :: Proxy (HLength (DOne ': ones)))))

-- | zero matrix. The size and dimension is determined by the type.
zeroes :: forall c a r e. (H.Field e,
                        HNat2Integral (HLength (a ': r)),
                        HNat2Integral (HLength (DOne ': c)))
    => D '(a, '[r, c]) e
zeroes = DMat (H.konst 0
        (hNat2Integral (Proxy :: Proxy (HLength (a ': r))),
         hNat2Integral (Proxy :: Proxy (HLength (DOne ': c)))))

type family CanAddConst (a :: k) (m :: [[k]]) :: Constraint
type instance CanAddConst a [as, ones] = (AllEq a as, AllEq '[] ones)
type instance CanAddConst a '[as] = (AllEq a as)

addConstant :: (H.Field e, CanAddConst a sh)
    => D '(a, '[]) e
    -> D '(a, sh) e
    -> D '(a, sh) e
addConstant (DScal a) (DMat b) = DMat (H.addConstant a b)
addConstant (DScal a) (DVec b) = DVec (H.addConstant a b)

conj :: D sh a -> D sh a
conj (DMat a) = DMat (H.conj a)
conj (DVec a) = DVec (H.conj a)

-- | conjugate transpose
ctrans x = conj . trans $ x

diag :: (MapConst DOne r ~ c, SameLength r c)
  => D '(a, '[r]) t -> D '(a, '[r,c]) t
diag (DVec a) = DMat (H.diag a)

#if MIN_VERSION_hmatrix(0,15,0)
-- | 'H.blockDiag'. The blocks should be provided as:
--
-- @blockDiag $ 'hBuild' m1 m2 m3@
--
-- only available if hmatrix >= 0.15
diagBlock :: (HMapOut UnDimMat (b ': bs) (H.Matrix e),
              Num e, H.Field e, AppendShOf b bs (D '(a, sh) e),
              sh ~ '[r,c])
  => HList (b ': bs)
  -> D '(a, sh) e
diagBlock pairs = DMat (H.diagBlock (hMapOut UnDimMat pairs))
#endif

data UnDimMat = UnDimMat
instance (D sh a ~ x, H.Matrix a ~ y) => ApplyAB UnDimMat x y where
        applyAB _ (DMat x) = x

class DiagBlock (bs :: [*]) t

-- | @AppendShOf a [b,c,d] aas@ makes aas have the type of a matrix that
-- has a,b,c,d along the diagonal
class AppendShOf a (as :: [*]) aas
instance 
 (e ~ f, f ~ g,
  AppendShOf (D ab e) ds z,
  AppendDims a b ab,
     
  -- constraints to force D in the type
  x ~ D a e,
  y ~ D b f,
  z ~ D c g) =>
  AppendShOf x (y ': ds) z 
instance (x ~ z) => AppendShOf x '[] z

type family AppendDims (a :: (*, [[*]])) (b :: (*, [[*]])) (c :: (*, [[*]])) :: Constraint
type instance AppendDims '(a, [ra,ca]) '(b, [rb,cb]) '(c, [rc,cc])
  = (c ~ a, AppendEq ra (b ': rb) rc, AppendEq ca cb cc)
-- how to handle vectors?
--type instance AppendDims '(a, '[ra]) '(b, '[rb]) = '(a, '[HAppendR ra (b ': rb)])

class ToHList sh e result where
    toHList :: D sh e -> HList result

-- | given a vector like @x = DimMat '[units] e@ this does something like
-- @[ (x \@> i) | i <- [1 .. n] ]@, if we had comprehensions for HLists
instance 
    (HListFromList e e1,
     SameLength result e1,
     HMapAux AddDimensional e1 result,
     ToHListRow (r ': rs) e result) =>
   ToHList '(r, '[rs]) e result where
  toHList (DVec v) = case hListFromList (H.toList v) :: HList e1 of
      e1 -> hMap AddDimensional e1

class HListFromList e e' where
        hListFromList :: [e] -> HList e'
instance HListFromList e '[] where
        hListFromList _ = HNil
instance (e ~ e', HListFromList e es) => HListFromList e (e' ': es) where
        hListFromList (e : es) = e `HCons` hListFromList es 

type family ToHListRow (a :: [*]) e (b :: [*]) :: Constraint
type instance ToHListRow (a ': as) e (b ': bs) = (Quantity a e ~ b, ToHListRow as e bs)

data AddDimensional = AddDimensional
instance (Quantity t x ~ y) => ApplyAB AddDimensional x y where
        applyAB _ x = Dimensional x

class FromHList list sh e where
  fromHList :: HList list -> D sh e

instance 
    (H.Field e,
     HMapOut RmDimensional list e,
     ToHListRow (r ': rs) e list) =>
  FromHList list '(r, '[rs]) e where
  fromHList xs = DVec (H.fromList (hMapOut RmDimensional xs))

data RmDimensional = RmDimensional
instance (x ~ Quantity d y) => ApplyAB RmDimensional x y where
        applyAB _ (Dimensional a) = a


class FromHLists lists sh e where
  fromHLists :: HList lists -> D sh e


-- | [[Dim e unit]] -> DimMat units e
instance 
  (ToHListRows' (r1 ': r) c e lists,
   HMapOut (HMapOutWith RmDimensional) lists [e],
   H.Field e) =>
  FromHLists lists '(r1, [r,c]) e where
    fromHLists xs = DMat (H.fromLists (hMapOut (HMapOutWith RmDimensional) xs))

newtype HMapOutWith f = HMapOutWith f
instance (HMapOut f l e, es ~ [e], HList l ~ hl) => ApplyAB (HMapOutWith f) hl es where
    applyAB (HMapOutWith f) = hMapOut f

class ToHListRows' (r :: [*]) (c :: [*]) (e :: *) (rows :: [*])
instance ToHListRows' '[] c e '[]

instance (ToHListRows' r c e rows,
          MultEq r c c',
          HMapCxt (AddQty e) (HList c') hListRow c' row')
  => ToHListRows' (r1 ': r) c e (hListRow ': rows)

data AddQty u
instance (qty ~ Quantity u e) => ApplyAB (AddQty u) e qty

class ToHLists sh e xs where
    toHLists :: D sh e -> HList xs

-- | DimMat units e -> [[Dim e unit]]
instance
    (HListFromList e e1,
     HListFromList (HList e1) e2,
     HMapAux (HMap AddDimensional) e2 xs,
     ToHListRows' ri ci e xs,
     SameLength e2 xs,
     (r1 ': r) ~ ri, (DOne ': c) ~ ci )
  => ToHLists '(r1, [r,c]) e xs where
  toHLists (DMat m) = case hListFromList (map hListFromList (H.toLists m) :: [HList e1]) :: HList e2 of
    e2 -> hMap (HMap AddDimensional) e2




{-




class PairsToList a t where
        pairsToList :: a -> [H.Matrix t]
instance PairsToList () t where
        pairsToList _ = []
instance (PairsToList b t, t' ~ t) => PairsToList (DimMat sh t',b) t where
        pairsToList (DimMat a,b) = a : pairsToList b

class EigV (sh :: [[ [DimSpec *] ]])
           (eigenValue  :: [[DimSpec *]])
           (eigenVector :: [[[DimSpec *]]])

instance
  ( SameLengths [r,c,r',c',rinv,cinv,eigval,erinv],
    -- ZipWithMul r c eigval,
    MapConst '[] r ~ eigval,
    PInv [r',c'] [rinv,cinv],
    -- AreRecips r' cinv,
    -- AreRecips c' rinv,
    cinv ~ c,
    c ~ ('[] ': _1),
    c' ~ ('[] ': _2),
    ZipWithMul eigval rinv erinv,
    MultiplyCxt [r',c'] erinv r,
    sh ~ [r,c],
    sh' ~ [r',c'])
    =>  EigV sh eigval sh'
-- | when no eigenvectors are needed
type family EigE (sh :: [[ [DimSpec *] ]]) (eigenValue  :: [ [DimSpec *] ]) :: Constraint
type instance EigE [r,c] eigval = ( SameLengths [r,c,eigval], ZipWithMul r c eigval)

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

{-
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
-}

wrapEig :: (c' ~ ('[] ': _1),
            EigV [r,c] eigVal [r',c'],
    H.Field y, H.Field z)
    => (H.Matrix x -> (H.Vector y, H.Matrix z)) ->
    DimMat [r,c] x ->
    (DimMat '[eigVal] y, DimMat [r',c'] z)
wrapEig hmatrixFun (DimMat a) = case hmatrixFun a of
    (e,v) -> (DimVec e, DimMat v)

wrapEigOnly :: (EigE [r,c] eigVal, H.Field y)
    => (H.Matrix x -> H.Vector y) ->
    DimMat [r,c] x -> DimMat '[eigVal] y
wrapEigOnly hmatrixFun (DimMat a) = case hmatrixFun a of
    (e) -> DimVec e

-}

-- | @\\a xs -> map (map (const a)) xs@
type family MapMapConst (a::k) (xs :: [[l]]) :: [[k]]
type instance MapMapConst a (x ': xs) = MapConst a x ': MapMapConst a xs
type instance MapMapConst a '[] = '[]

-- | @\\a xs -> map (const a) xs@
type family MapConst (a :: k) (xs :: [l]) :: [k]
type instance MapConst a (x ': xs) = a ': MapConst a xs
type instance MapConst a '[] = '[]

-- | @map fst@
{-
type family MapFst (a :: *) :: [*]
type instance MapFst ((a,_t) , as) = a ': MapFst as
type instance MapFst () = '[]
-}

-- | convert from (a,(b,(c,(d,())))) to '[a,b,c,d]
type family FromPairs (a :: *) :: [*]
type instance FromPairs (a,b) = a ': FromPairs b
type instance FromPairs () = '[]

{-
-- | @\\a xs -> map (/a) xs@
type family MapDiv (a :: k) (xs :: [k]) :: [k]
type instance MapDiv a (x ': xs) = (x @- a) ': MapDiv a xs
type instance MapDiv a '[] = '[]

type family UnDQuantity (a :: [*]) :: [ [*] ]
type instance UnDQuantity (x ': xs) = UnDQuantity1 x ': UnDQuantity xs
type instance UnDQuantity '[] = '[]

type family UnDQuantity1 (a :: *) :: [*] 
type instance UnDQuantity1 (Unit t x) = x

type family DimMatFromTuple ijs :: * -> *
type instance DimMatFromTuple ijs =
        DimMat [UnDQuantity (MapFst ijs),
               '[] ': MapDiv (UnDQuantity1 (Fst (Fst ijs)))
               (UnDQuantity (FromPairs (Snd (Fst ijs))))]
type family Append (a :: [k]) (b :: [k]) :: [k]
type instance Append (a ': as) b = a ': Append as b
type instance Append '[] b = b
-}
