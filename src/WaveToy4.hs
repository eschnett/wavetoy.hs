{-# LANGUAGE DeriveFoldable, DeriveFunctor, TypeFamilies #-}

module WaveToy4
    (Cell(..), initialCell, energyCell, bndCell, rhsCell, errorCell,
     Grid(..), initialGrid, energyGrid, bndGrid, bndsGrid, rhsGrid, errorGrid,
     rk2, rk4)
    where

import Prelude hiding (zip, zipWith, zipWith3, unzip)
import Control.Exception.Base
import Data.Array
import Data.Array.IArray (amap)
import Data.Key hiding ((!))
import Data.AdditiveGroup
import Data.Fixed (mod')
import Data.VectorSpace



-- State vector

data Cell a = Cell { x, u, ρ, v :: a }
              deriving (Eq, Foldable, Functor, Read, Show)

instance Zip Cell where
    zip s s' = Cell (x s, x s') (u s, u s') (ρ s, ρ s') (v s, v s')

instance (AdditiveGroup a, Num a) => AdditiveGroup (Cell a) where
    zeroV = Cell zeroV zeroV zeroV zeroV
    (^+^) = zipWith (^+^)
    negateV = fmap negateV

instance (VectorSpace a, Num a) => VectorSpace (Cell a) where
    type Scalar (Cell a) = a
    (*^) alpha = fmap (alpha *)

standingWave :: Floating a => a -> a -> Cell a
standingWave t x = Cell x u ρ v
    where k = pi
          omega = sqrt(k^2)
          amp = 1 / k
          u = amp * sin(omega * t) * sin(k * x)
          ρ = amp * omega * cos(omega * t) * sin(k * x)
          v = amp * k * sin(omega * t) * cos(k * x)

planeWave :: Floating a => a -> a -> Cell a
planeWave t x = Cell x u ρ v
    where k = 2 * pi
          omega = sqrt(k^2)
          amp = 1 / (2 * pi)
          u = amp * sin(omega * t - k * x)
          ρ = amp * omega * cos(omega * t - k * x)
          v = - amp * k * cos(omega * t - k * x)

analyticCell :: Floating a => a -> a -> Cell a
analyticCell = standingWave

initialCell :: Floating a => a -> a -> Cell a
initialCell = analyticCell

energyCell :: Fractional a => Cell a -> a
energyCell s = 1/2 * ((ρ s)^2 + (v s)^2)

bndCell :: Real a => Bool -> Cell a -> Cell a -> Cell a
bndCell f s b = Cell x' (-u s) (-ρ s) (v s)
    where x' = 2 * x s - x b
-- bndCell f s b = Cell x' (u b) (ρ b) (v b)
--     where x' = x b + if not f then -1 else 1

rhsCell :: Fractional a => (Cell a, Cell a) -> Cell a -> Cell a
rhsCell (bm, bp) s = Cell 0 (ρ s) (ddx v) (ddx ρ)
    where ddx var = (var bp - var bm) / (x bp - x bm)

errorCell :: (AdditiveGroup a, Floating a) => a -> Cell a -> Cell a
errorCell t s = s ^-^ analyticCell t (x s)



-- Instance declarations for Array

instance Ix i => Zip (Array i) where
    zip a a' = listArray newbnds (zip (elems a) (elems a'))
        where newbnds = assert (bounds a == bounds a') (bounds a)

instance (AdditiveGroup e, Ix i, Num i) => AdditiveGroup (Array i e) where
    zeroV = listArray (1, 0) []
    (^+^) = zipWith (^+^)
    negateV = amap negateV

instance (VectorSpace e, Ix i, Num i) => VectorSpace (Array i e) where
    type Scalar (Array i e) = Scalar e
    (*^) alpha = amap (alpha *^)



-- Grid

newtype Grid a = Grid { getGrid :: Array Int a }
    deriving (Foldable, Functor, Read, Show)

instance Zip Grid where
    zip (Grid g) (Grid g') = Grid $ zip g g'

instance AdditiveGroup a => AdditiveGroup (Grid a) where
    zeroV = Grid $ listArray (1,0) []
    (^+^) = zipWith (^+^)
    negateV = fmap negateV

instance VectorSpace a => VectorSpace (Grid a) where
    type Scalar (Grid a) = Scalar a
    (*^) alpha = fmap (alpha *^)

initialGrid :: Floating a => Int -> a -> Grid (Cell a)
initialGrid n t = Grid $ listArray bnds elts
    where bnds = (0, n-1)
          elts = [initialCell t (x i) | i <- range bnds]
          x i = (fromIntegral i + 1/2) / fromIntegral n

energyGrid :: (AdditiveGroup a, Fractional a) => Grid (Cell a) -> a
energyGrid g = getSum $ foldMap (Sum . energyCell) g

bndGrid :: Real a => Bool -> Grid (Cell a) -> Cell a
bndGrid f (Grid g) = if not f
                     then bndCell f (g!lbnd) (g!(lbnd+1))
                     else bndCell f (g!ubnd) (g!(ubnd-1))
    where (lbnd, ubnd) = bounds g
-- bndGrid f (Grid g) = if not f
--                      then bndCell f (g!lbnd) (g!ubnd)
--                      else bndCell f (g!ubnd) (g!lbnd)
--     where (lbnd, ubnd) = bounds g

bndsGrid :: Real a => Grid (Cell a) -> (Cell a, Cell a)
bndsGrid g = (bndGrid False g, bndGrid True g)

rhsGrid :: Fractional a => (Cell a, Cell a) -> Grid (Cell a) -> Grid (Cell a)
rhsGrid bs (Grid g) = Grid $ listArray bnds elts
    where bnds = bounds g
          elts = if rangeSize bnds == 0
                 then []
                 else if rangeSize bnds == 1
                 then [rhs bm bp (g!lbnd)]
                 else ([rhs bm (bndp lbnd) (g!lbnd)] ++
                       [rhs (bndm i) (bndp i) (g!i) | i <- [lbnd+1..ubnd-1]] ++
                       [rhs (bndm ubnd) bp (g!ubnd)])
          (lbnd, ubnd) = bnds
          (bm, bp) = bs
          bndm i = g ! (i-1)
          bndp i = g ! (i+1)
          rhs bm bp elt = rhsCell (bm, bp) elt

errorGrid :: (AdditiveGroup a, Floating a) =>
             a -> Grid (Cell a) -> Grid (Cell a)
errorGrid t = fmap (errorCell t)



-- Integrator

rk2 :: (VectorSpace v, s ~ Scalar v, Fractional s) => s -> (v -> v) -> v -> v
rk2 h f y0 = let k1 = f y0
                 k2 = f (y0 ^+^ (h/2) *^ k1)
             in y0 ^+^ h *^ k2

rk4 :: (VectorSpace v, s ~ Scalar v, Fractional s) => s -> (v -> v) -> v -> v
rk4 h f y0 = let k1 = f y0
                 k2 = f (y0 ^+^ (h/2) *^ k1)
                 k3 = f (y0 ^+^ (h/2) *^ k2)
                 k4 = f (y0 ^+^ h *^ k3)
             in (y0 ^+^
                 (h/6) *^ k1 ^+^ (h/3) *^ k2 ^+^ (h/3) *^ k2 ^+^ (h/6) *^ k4)
