-- | see examples/
module DimMat (
    module DimMat.Internal,
    module DimMat.QQ,
    {-
    module Data.Dimensions,
    -}
    module Data.HList.CommonMain,
    -- * to keep types looking ok
    Complex, 
    -- Dim,
    ) where
    
import DimMat.Internal
import DimMat.QQ
{-
import Data.Dimensions
import Data.Dimensions.Unsafe (Dim)
-}
import Data.HList.CommonMain
import Data.Complex
