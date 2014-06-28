{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TemplateHaskell #-}
-- | a variation on <http://hackage.haskell.org/package/hmatrix-syntax>,
-- which constructs a matrix which has dimensions stored.
module DimMat.QQ (matD,blockD) where

-- import Data.Packed.Syntax.Internal
import Language.Haskell.TH


import Language.Haskell.TH as TH
import Language.Haskell.TH.Quote as TH
import Language.Haskell.TH.Syntax as TH

import Data.HList.CommonMain (hZero,hSucc)

import DimMat.Internal
import Data.Packed.Vector(
  Vector,
  )
import Data.Packed.Matrix(
  Matrix,
  rows,
  cols,
  )
import Data.Packed.ST(
  runSTVector,
  newUndefinedVector,
  unsafeWriteVector,
  runSTMatrix,
  newUndefinedMatrix,
  unsafeWriteMatrix,
 )
import Data.Packed.Development(
  MatrixOrder(..),
  at',
  atM',
  )
import Data.List(
  transpose
  )

import Data.Packed.Syntax.Internal
import Data.Proxy
import Language.Haskell.TH.Quote

import Numeric.Units.Dimensional (Dimensional(..))



matD = QuasiQuoter {
  quoteExp = \s -> case matListExp s of
    Right (_, _, x) -> buildMatST (map (map return) x)
    Left msg -> fail msg,
  quotePat = \s -> case matListPat s of
    Left msg -> fail msg
    Right (ni,nj,x) -> do
        l <- newName "l"
        viewP (lamE [varP l]
            (nestedTupE [ [| $(varE l) @@> ($(mkNat i),$(mkNat j)) |]
                        | i <- [0 .. ni-1], j <- [0 .. nj-1]]))
            (nestedTupP (map return (concat x))),
  quoteDec = error "matD",
  quoteType = error "matD"
 }

blockD = QuasiQuoter {
  quoteExp = \s -> case matListExp s of
    Right (_,_,x) -> buildBlockST (map (map return) x)
    Left msg -> fail msg,
  quotePat = error "blockD",
  quoteDec = error "blockD",
  quoteType = error "blockD"
}

nestedTupE = foldr (\a b -> tupE [a,b]) [| () |]
nestedTupP = foldr (\a b -> tupP [a,b]) [p| () |]

toHListE = foldr (\a b -> [| ($a,$b) |]) [| () |]
varPs = foldr (\a b -> tupP [a,b]) wildP
matrixNames = zipWith (\i -> zipWith (\j _ -> vij i j) [0 .. ])  [0 .. ]

toHListP = varPs . map (varPs . map varP) . matrixNames
vij i j = mkName ("x"++show i++"_"++show j)

buildBlockST :: [[Q TH.Exp]] -> Q TH.Exp
buildBlockST es =
  do
    let rows = map (foldWithNamedFunction 'hconcat) es
    let fullHFirst = foldWithNamedFunction 'vconcat rows
    let cols = map (foldWithNamedFunction 'vconcat) (transpose es)
    let fullVFirst = foldWithNamedFunction 'hconcat cols
    [| $(fullHFirst) `asTypeOf` $(fullVFirst) |]

foldWithNamedFunction :: Name -> [Q TH.Exp] -> Q TH.Exp
foldWithNamedFunction f [] = undefined
foldWithNamedFunction f [e] =  e
foldWithNamedFunction f (e:es) = [| $(varE f) $(e) $(foldWithNamedFunction f es) |]

-- adapted from Data.Packed.Syntax
buildMatST :: [[Q TH.Exp]] -> Q TH.Exp
buildMatST es = let r = length es
                    c = length (head es)
  in do
 p <- newName "p"
 caseE (toHListE (map toHListE es))
  [match (p `asP` toHListP es) 
    (normalB $ do
      m <- newName "m"
      letE
        [valD (varP m)
          (do
            n <- newName "n"
            normalB
              [| (DMat $ runSTMatrix $(doE $
                    bindS (varP n) [| newUndefinedMatrix RowMajor r c |] :
                    [ noBindS [| unsafeWriteMatrix $(varE n) i j
                                    ((\(Dimensional a) -> a) $(varE (vij i j))) |]
                    | i <- [0 .. r-1], j <- [0 .. c-1]] ++
                    [ noBindS [| return $(varE n) |] ])
                 ) `asProxyTypeOf` $(varE 'toDM) $(varE p) |])
          [] ]
        (foldr
                  (\a b -> [| $b `const` $a |])
                  [| $(varE m) `asProxyTypeOf` $(varE 'toDM) $(varE p) |]
                  -- a check that lookups actually work
                  [ [| $(varE (vij i j)) `asTypeOf`
                          ($(varE m) @@> ($(mkNat i),$(mkNat j))) |] 
                      | i <- [0 .. r-1],
                        j <- [0 .. c-1] ]
          )
    )
    []
  ]

mkNat :: Int -> ExpQ
mkNat 0 = [| hZero |]
mkNat n = [| hSucc $(mkNat (n-1)) |]

promotedListT :: [Q Type] -> Q Type
promotedListT = foldr (\a b -> [t| $promotedConsT $a $b|]) promotedNilT
