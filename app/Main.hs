module Main where

import Prelude hiding (zip, zipWith, zipWith3, unzip)
import Control.Monad
import Control.Monad.Loops
import Data.Key hiding ((!))

import WaveToy4



data State = State { iter :: Int,
                     time :: Double,
                     grid :: Grid (Cell Double) }

rhs g = rhsGrid (bndsGrid g) g
next dt (State i t g) = State (i + 1) (t + dt) (rk2 dt rhs g)



-- Parameters
icfl = 2
rf = 1
nx = rf * 10
nt = 10 * icfl * nx
ninfo = icfl * nx
nfile = icfl * nx

dx = 1 / fromIntegral nx
dt = dx / fromIntegral icfl



output ninfo nfile dx (State i t g) = do
  let en = dx * energyGrid g
  when (ninfo > 0 && mod i ninfo == 0) $ do
    putStrLn $ ""
    putStrLn $ "iteration: " ++ show i
    putStrLn $ "time: " ++ show t
    -- putStrLn $ "grid: " ++ show g
    putStrLn $ "energy: " ++ show en
  when (nfile > 0 && mod i nfile == 0) $ do
    (if i == 0 then writeFile else appendFile) "wavetoy.out" (showGrid g)
    where showCell (Cell x u ρ v) =
              show x ++ "\t" ++ show u ++ "\t" ++ show ρ ++ "\t" ++ show v
          showLine (s, e) =
              show i ++ "\t" ++ show t ++ "\t" ++
              showCell s ++ "\t" ++
              showCell e ++ "\n"
          showGrid g =
              "# iter time   x u ρ v   ex eu eρ ev\n" ++
              foldMap showLine (zip g (errorGrid t g)) ++
              "\n"

drive = do
  -- Initial conditions
  let i = 0
  let t = 0
  let g = initialGrid nx t
  let s = State i t g
  output ninfo nfile dx s
  -- Evolution loop
  -- s <- iterateUntilM
  --      (\s -> iter s >= nt)
  --      (\s -> do let s' = step dt dx s
  --                output ninfo nfile dx s'
  --                return s')
  --      s
  s <- let done s = iter s >= nt
           step s = do let s' = next dt s
                       output ninfo nfile dx s'
                       return s'
       in iterateUntilM done step s
  return ()



main :: IO ()
main = do putStrLn "WaveToy in Haskell"
          putStrLn "Erik Schnetter <schnetter@gmail.com>"
          drive
          -- Done.
          putStrLn "\nDone."
