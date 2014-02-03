# compile-time checked units for matrices
Two related generalizations of `dimensional-tf` are explored in this package:

* provide hmatrix operations that know about units too
* automatic differentiation (`ad` package) that knows about units

This comes at a cost, since type errors become more complicated, and type signatures
for functions that operate on the involved quantities are more difficult (see
`examples/controlspace.hs`). In return, you get a stronger assurance that your
program does the right thing.

Haddocks are available at http://aavogt.github.io/haddock/DimMat

## related work
* [dimensional-vectors](https://github.com/bjornbm/dimensional-vectors) is smaller, and uses a `[[a]]` representation of the data instead of `Data.Packed.Matrix a`
* http://www.haskell.org/haskellwiki/Physical_units

### matrix dimensions statically checked
These packages provide operations where the typechecker will prevent invalid operations, such as multiplying an m×n matrix with a p×q matrix when n /= p, at compile-time.

* http://hackage.haskell.org/package/vector-static
* http://hackage.haskell.org/package/linear
* http://hackage.haskell.org/package/Vec
* http://ofb.net/~frederik/vectro/
* [vector-space](http://hackage.haskell.org/package/vector-space) has limited operations on tuples (up to 4 elements)
* [tensor](http://hackage.haskell.org/package/tensor) has the number of indices at type level, but the range over which index varies is checked at runtime

### the record connection
[extensible records](http://www.haskell.org/haskellwiki/Extensible_record) have many needs in common with DimMat. Types in HMatrix like `fromBlocks :: [[Matrix t]] -> Matrix t` are (or will be) generalized to use [HList](http://hackage.haskell.org/package/HList) instead of ordinary lists (`[]`).
