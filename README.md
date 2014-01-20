# compile-time checked units for matrices
`dimensional-tf` works well, but it only works for scalars. Two related generalizations:

* automatic differentiation (`ad` package) that knows about units
* matrices/vectors with units

also included is a quasiquote (for expressions and patterns). See `examples/qq.hs`

## related work
### matrix dimensions statically checked
http://hackage.haskell.org/package/vector-static

http://hackage.haskell.org/package/Vec

### other
HList
