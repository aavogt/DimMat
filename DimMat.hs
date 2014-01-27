-- | see examples/
module DimMat (
    module DimMat.Internal,
    module DimMat.QQ,
    module Numeric.Units.Dimensional.TF.Prelude,
    module Data.HList.CommonMain,
    -- * to keep types looking ok
    Z,S,N,Zero,
    Complex,
    ) where
    
import DimMat.Internal
import DimMat.QQ
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.NumType.TF
import Data.HList.CommonMain
import Data.Complex
