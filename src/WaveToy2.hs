{-# LANGUAGE DeriveFunctor, TypeFamilies #-}

module WaveToy2
    (State, initialState, bndState, bndsState, rhsState,
     StateArray, initialArray, bndArray, bndsArray, rhsArray,
     rk2)
    where

import Prelude hiding (zip, zipWith, zipWith3, unzip)
import Data.Array
import Data.Key hiding ((!))

import VectorSpace



data State a = State { u :: a,
                       rho :: a,
                       v :: a }
               deriving (Functor, Read, Show)

instance Zip State where
    zip s s' = State (u s, u s') (rho s, rho s') (v s, v s')

instance (Num a, VectorSpace a) => VectorSpace (State a) where
    zeroV = State zeroV zeroV zeroV
    (^+^) = zipWith (^+^)
    negateV = fmap negateV
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



newtype StateArray a = StateArray { getArray :: Array Int (State a) }
    deriving (Functor, Read, Show)

instance Zip StateArray where
    zip (StateArray s) (StateArray s') = StateArray $ zipWith zip s s'

instance (Num a, VectorSpace a) => VectorSpace (StateArray a) where
    zeroV = StateArray $ listArray (1,0) []
    -- (^+^) = zipWith (zipWith (^+^))
    StateArray s ^+^ StateArray s' = StateArray $ zipWith (zipWith (^+^)) s s'
    negateV (StateArray s) = StateArray $ fmap (fmap negateV) s
    type Scalar (StateArray a) = a
    alpha *^ (StateArray s) = StateArray $ fmap (fmap (alpha *)) s

initialArray :: Floating a => Int -> StateArray a
initialArray n = StateArray $ listArray bnds elts
    where bnds = (0, n-1)
          inds = [0..n-1]
          elts = [initialState (x i) | i <- inds]
          x :: Fractional a => Int -> a
          x i = (fromIntegral i + 1/2) / fromIntegral n

bndArray :: Num a => Bool -> StateArray a -> State a
bndArray f s = let (lbnd, ubnd) = bounds $ getArray s
                   elt = getArray s ! if f then ubnd else lbnd
             in bndState f elt

bndsArray :: Num a => StateArray a -> (State a, State a)
bndsArray s = (bndArray False s, bndArray True s)

rhsArray :: Fractional a => (State a, State a) -> StateArray a -> StateArray a
rhsArray bs (StateArray s) = StateArray $ listArray bnds elts
    where bnds = bounds s
          elts = [rhs bm (bndp lbnd) (s!lbnd)] ++
                 [rhs (bndm i) (bndp i) (s!i) | i <- [lbnd+1..ubnd-1]] ++
                 [rhs (bndm ubnd) bp (s!ubnd)]
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
