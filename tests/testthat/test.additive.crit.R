library(OptimalDesign)
context("Two factor additive model")

# Make a dataset
X1 <- seq(-1, 1, .01)
X2 <- seq(-1, 1, .01)
X <- expand.grid(X1, X2)
names(X) <- c("X1", "X2")
formula <- ~X1 + X2

test_that("D > 99%", {
              expect_more_than(
                           getEff(formula,
                                  optimalDesign(formula, X, 4, "D", 1000000),
                                  criteria = "D")$D,
                           99)
})
test_that("A > 73", {
              expect_more_than(
                           getEff(formula,
                                  optimalDesign(formula, X, 4, "A", 1000000),
                                  criteria = "A")$A,
                           73)
})

test_that("I > 1", {
              expect_more_than(
                           getEff(formula,
                                  optimalDesign(formula, X, 4, "I", 1000000),
                                  evaluation = X,
                                  criteria = "I")$I,
                           1)
})

test_that("G > 49.99", {
              expect_more_than(
                           getEff(formula,
                                  optimalDesign(formula, X, 4, "D", 1000000),
                                  evaluation = X,
                                  criteria = "G")$G,
                           49.99)
})
