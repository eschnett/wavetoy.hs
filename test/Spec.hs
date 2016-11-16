import Test.HUnit

import Data.AdditiveGroup

import WaveToy4

testCellError = TestCase $ assertEqual "error," want have
    where want = zeroV::Cell Double
          have = errorCell 0 (initialCell 0 0)

tests = TestList [testCellError]

main :: IO ()
main = do counts <- runTestTT tests
          assert $ errors counts + failures counts == 0
          return ()
