-- | see examples/
module DimMat (

   -- * Quasiquotes
   matD,
   blockD,
 

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
   ToHLists(toHLists), toHList, FromHLists(fromHLists), fromHList,
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
   mulMat, mulVec,
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
   Dot(..), 
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
   -- eigs
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
   -- ad
#ifdef WITH_Ad
   diff,
#endif


    -- * todo arrange
    DotS,
    MultEq,

    -- * to keep types looking ok
    D,
    module Data.HList.CommonMain,
    Complex, 

    -- ** "Numeric.NumType"
    Pos, Neg, Succ, Zero, Neg1,
    ) where
    
import Numeric.NumType
import DimMat.Internal
import DimMat.QQ
import Data.HList.CommonMain
import Data.Complex
