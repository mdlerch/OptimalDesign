library(OptimalDesign)
context("Check efficiency calculations")

# Make a dataset
X1 <- genFactorial(n.levels = 2, n.factors = 2)
formula1 <- ~X1 + X2

# Large evaluation set
Xe <- genFactorial(n.levels = 100, n.factors = 2) / 99

test_that("2^2 factorial", {
              expect_equal(
                           getEff(formula1, design = X1, criteria = "D")$D,
                           100)
              expect_equal(
                           getEff(formula1, design = X1, criteria = "A")$A,
                           100)
              expect_more_than(
                           getEff(formula1, design = X1, evaluation = Xe, criteria = "I")$I,
                           1.28)
              expect_equal(
                           getEff(formula1, design = X1, evaluation = Xe, criteria = "G")$G,
                           50)
})
