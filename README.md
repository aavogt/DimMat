# Compile-time checked units for matrices
Two related generalizations of `dimensional` are explored in this package:

* wrap hmatrix operations to check physical units and dimensions at compile time.
* automatic differentiation (`ad` package) that knows about units

This comes at a cost, since type errors become more complicated, and type signatures
for functions that operate on the involved quantities are more difficult (see
`examples/controlspace.hs`). In return, you get a stronger assurance that the
program does the right thing.

Additionally, there is an emphasis on type inference going both ways:
in ordinary haskell this is pretty much provided. For example if we have,
`z = x + y`, then if due to other parts of the program the compiler knows
`x :: (Int,_1,_2)`, `y :: (_3,Int,_4)` and `z :: (_5,_6,Int)`, it can
conclude that (`+`) is needed with type `(Int,Int,Int) -> (Int,Int,Int) -> (Int,Int,Int)` [^note].
This means using things with kind Constraint (class/type/type family) to


[^note]: which isn't available by default, but you can check http://hackage.haskell.org/package/NumInstances or write it by hand

Haddocks are available at http://aavogt.github.io/haddock/DimMat

# Installation
Get ghc-7.6 or ghc-7.8 and cabal-install. Then:
```sh
cabal install cabal-meta cabal-src
git clone https://github.com/aavogt/DimMat
cd DimMat/
cabal-meta install
```

# Related Work
## Units and AD or linear algebra
* [dimensional-vectors](https://github.com/bjornbm/dimensional-vectors) is smaller, and uses a `[[a]]` representation of the data instead of `Data.Packed.Matrix a`
* https://github.com/dmcclean/dimensional-dk-linalg shares some code with DimMat, but 
* https://github.com/bjornbm/dimensional-experimental/blob/master/Numeric/Units/Dimensional/AD.hs has the same wrapper for AD

## Matrix dimensions statically checked
These packages provide operations where the typechecker will prevent invalid operations, such as multiplying an m×n matrix with a p×q matrix when n /= p, at compile-time.

* http://hackage.haskell.org/package/vector-static
* http://hackage.haskell.org/package/linear
* http://hackage.haskell.org/package/Vec
* http://hackage.haskell.org/package/hmatrix-static
* [vectro](http://ofb.net/~frederik/vectro/) is a relatively early library
* [vector-space](http://hackage.haskell.org/package/vector-space) has limited operations on tuples (up to 4 elements)
* [tensor](http://hackage.haskell.org/package/tensor) has the number of indices at type level, but the range over which index varies is checked at runtime some of the time
* [storable-static-array](http://github.com/chowells79/storable-static-array) uses ghc-7.8 features

## units
* The [haskellwiki physical units](http://www.haskell.org/haskellwiki/Physical_units) mentions several methods to check units at compile time or at run-time
* https://github.com/haasn/units

## Records
[extensible records](http://www.haskell.org/haskellwiki/Extensible_record) have many needs in common with DimMat. Types in HMatrix like `fromBlocks :: [[Matrix t]] -> Matrix t` are (or will be) generalized to use [HList](http://hackage.haskell.org/package/HList) instead of ordinary lists (`[]`).
