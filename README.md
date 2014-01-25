# compile-time checked units for matrices
While `dimensional-tf` works well, it only works for scalars. Two related generalizations are explored in this package:

* automatic differentiation (`ad` package) that knows about units
* wrap operations provided by HMatrix so they know about units too

also included is a quasiquote (for expressions and patterns). See `examples/qq.hs`

Documentation at http://aavogt.github.io/haddock/DimMat

## related work
### matrix dimensions statically checked
http://hackage.haskell.org/package/vector-static

http://hackage.haskell.org/package/Vec

### other
HList
