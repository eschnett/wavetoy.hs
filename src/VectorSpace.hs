{-# LANGUAGE TypeFamilies #-}

module VectorSpace
    (VectorSpace(..))
        where

import Prelude hiding (zip, zipWith, zipWith3, unzip)
import Control.Exception.Base
-- import Data.Align
import Data.Array
import Data.Array.IArray (amap)
import Data.Key hiding ((!))



-- Modelled after http://hackage.haskell.org/package/vector-space-0.10.4

infixl 6 ^+^, ^-^
infixr 7 *^

class VectorSpace v where
    zeroV :: v
    (^+^) :: v -> v -> v
    negateV :: v -> v
    (^-^) :: v -> v -> v
    x ^-^ y = x ^+^ negateV y
    type Scalar v :: *
    (*^) :: Scalar v -> v -> v
    -- negateV x = (-1) *^ x

instance VectorSpace Double where
    zeroV = 0
    (^+^) = (+)
    negateV = negate
    type Scalar Double = Double
    (*^) = (*)

-- for fun
instance VectorSpace () where
    zeroV = ()
    () ^+^ () = ()
    negateV () = ()
    type Scalar () = ()
    () *^ () = ()

-- for fun
instance (VectorSpace v, VectorSpace w, Scalar v ~ Scalar w) =>
    VectorSpace (v, w) where
        zeroV = (zeroV, zeroV)
        (x,y) ^+^ (x',y') =  (x ^+^ x', y ^+^ y')
        -- (^+^) = zipWith (^+^)
        negateV = fmap negateV
        type Scalar (v, w) = Scalar v
        (*^) alpha = fmap (alpha *^)



instance Ix i => Zip (Array i) where
    zip a a' = listArray newbnds (zip (elems a) (elems a'))
        where newbnds = assert (bounds a == bounds a') (bounds a)

instance (VectorSpace e, Ix i, Num i) => VectorSpace (Array i e) where
        zeroV = listArray (1, 0) []
        xs ^+^ ys = listArray (bounds xs) [xs!i ^+^ ys!i | i <- indices xs]
        -- (^+^) = zipWith (^+^)
        negateV = amap negateV
        type Scalar (Array i e) = Scalar e
        (*^) alpha = amap (alpha *^)
