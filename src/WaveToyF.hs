{-# LANGUAGE DeriveFunctor, TypeFamilies #-}

module WaveToy4
    (State, -- StateBnd, StateBnds,
     initialState, bndState, bndsState, rhsState,
     Grid, -- StateArrayBnd, StateArrayBnds,
     initialArray, bndArray, bndsArray, rhsArray,
     rk2)
    where

import Prelude hiding (zip, zipWith, zipWith3, unzip)
import Control.Exception.Base
import Data.Array
import Data.Array.IArray (amap)
import Data.Key hiding ((!))
import Data.AdditiveGroup
import Data.VectorSpace



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



-- State vector

data State a = State { u :: a,
                       rho :: a,
                       v :: a }
               deriving (Functor, Read, Show)

instance Zip State where
    zip s s' = State (u s, u s') (rho s, rho s') (v s, v s')

instance (AdditiveGroup a, Num a) => AdditiveGroup (State a) where
    zeroV = State zeroV zeroV zeroV
    (^+^) = zipWith (^+^)
    negateV = fmap negateV

instance (VectorSpace a, Num a) => VectorSpace (State a) where
    type Scalar (State a) = a
    (*^) alpha = fmap (alpha *)

initialState :: Floating a => a -> State a
initialState x = State 0 (sin (2 * pi * x)) 0

bndState :: Num a => Bool -> State a -> State a
bndState f s = State (-u s) (-rho s) (v s)

bndsState :: Num a => State a -> (State a, State a)
bndsState s = (bndState False s, bndState True s)

rhsState :: Fractional a => (State a, State a) -> State a -> State a
rhsState bs s = State (rho s) (dx v bs s) (dx rho bs s)
    where dx var (bm, bp) s = (var bp - var bm) / 2



-- Grid

newtype Grid a = Grid { getGrid :: Array Int a }
    deriving (Functor, Read, Show)

instance Zip Grid where
    zip (Grid g) (Grid g') = Grid $ zip g g'

instance AdditiveGroup a => AdditiveGroup (Grid a) where
    zeroV = Grid $ listArray (1,0) []
    (^+^) = zipWith (^+^)
    negateV = fmap negateV

instance VectorSpace a => VectorSpace (Grid a) where
    type Scalar (Grid a) = Scalar a
    (*^) alpha = fmap (alpha *)

initialGrid :: Floating a => Int -> Grid (State a)
initialGrid n = Grid $ listArray bnds elts
    where bnds = (0, n-1)
          inds = [0..n-1]
          elts = [initialState (x i) | i <- inds]
          x :: Fractional a => Int -> a
          x i = (fromIntegral i + 1/2) / fromIntegral n

bndArray :: Num a => Bool -> Grid (State a) a -> State a
bndArray f (Grid g) = let (lbnd, ubnd) = bounds g
                          bnd = if f then ubnd else lbnd
                          elt = g ! bnd
                      in bndState f elt

bndsArray :: Num a => Grid (State a) -> (State a, State a)
bndsArray g = (bndArray False g, bndArray True g)

rhsArray :: Fractional a =>
            (State a, State a) -> Grid (State a) -> Grid (State a)
rhsArray bs (Grid g) = Grid $ listArray bnds elts
    where bnds = bounds g
          elts = [rhs bm (bndp lbnd) (g!lbnd)] ++
                 [rhs (bndm i) (bndp i) (g!i) | i <- [lbnd+1..ubnd-1]] ++
                 [rhs (bndm ubnd) bp (g!ubnd)]
          (lbnd, ubnd) = bnds
          (bm, bp) = bs
          bndm i = s ! (i-1)
          bndp i = s ! (i+1)
          rhs bm bp elt = rhsState (bm, bp) elt



rk2 :: (VectorSpace v, s ~ Scalar v, Fractional s) => s -> (v -> v) -> v -> v
rk2 dt rhs s0 = let r0 = rhs s0
                    s1 = s0 ^+^ (dt/2) *^ r0
                    r1 = rhs s1
                    s2 = s0 ^+^ dt *^ r1
                in s2



{-
data State a = State { u :: a,
                       rho :: a,
                       v :: a }
               deriving (Read, Show)

newtype StateBnd a = StateBnd (State a)
    deriving (Read, Show)
newtype StateBnds a = StateBnds (StateBnd a, StateBnd a)
    deriving (Read, Show)

initialState :: Floating a => a -> State a
initialState x = State 0 (sin (2 * pi * x)) 0

bndState :: Num a => Bool -> State a -> StateBnd a
bndState f s = StateBnd $ State (-u s) (-rho s) (v s)

bndsState :: Num a => State a -> StateBnds a
bndsState s = StateBnds (bndState False s, bndState True s)

rhsState :: Num a => StateBnds a -> State a -> State a
rhsState bs s = State (rho s) (dx v bs s) (dx rho bs s)
    where dx var (StateBnds (StateBnd bm, StateBnd bp)) s = var bp - var bm



newtype StateArray a = StateArray { getArray :: Array Int (State a) }
    deriving (Read, Show)

newtype StateArrayBnd a = StateArrayBnd (StateBnd a)
    deriving (Read, Show)
newtype StateArrayBnds a = StateArrayBnds (StateArrayBnd a, StateArrayBnd a)
    deriving (Read, Show)

-- deriving instance Read a => Read (StateArray a)
-- deriving instance Show a => Show (StateArray a)

initialArray :: Floating a => Int -> StateArray a
initialArray n = StateArray $ array bnds elts
    where bnds = (0, n-1)
          inds = [0..n-1]
          elts = [(i, initialState (x i)) | i <- inds]
          x :: Fractional a => Int -> a
          x i = (fromIntegral i + 1/2) / fromIntegral n

bndArray :: Num a => Bool -> StateArray a -> StateArrayBnd a
bndArray f s = let (lbnd, ubnd) = bounds $ getArray s
                   elt = getArray s ! if f then ubnd else lbnd
             in StateArrayBnd $ bndState f elt

bndsArray :: Num a => StateArray a -> StateArrayBnds a
bndsArray s = StateArrayBnds (bndArray False s, bndArray True s)

rhsArray :: Num a => StateArrayBnds a -> StateArray a -> StateArray a
rhsArray bs (StateArray s) = StateArray $ array bnds elts
    where bnds = bounds s
          elts = [(lbnd, rhs bm (bndp lbnd) (s!lbnd))] ++
                 [(i, rhs (bndm i) (bndp i) (s!i)) | i <- [lbnd+1..ubnd-1]] ++
                 [(ubnd, rhs (bndm ubnd) bp (s!ubnd))]
          (lbnd, ubnd) = bnds
          StateArrayBnds (StateArrayBnd bm, StateArrayBnd bp) = bs
          bndm i = StateBnd $ s ! (i-1)
          bndp i = StateBnd $ s ! (i+1)
          rhs bm bp elt = rhsState (StateBnds (bm, bp)) elt
-}
