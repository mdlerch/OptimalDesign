library(OptimalDesign)
context("Two factor additive model")

# Make a dataset
X1 <- seq(-1, 1, .01)
X2 <- seq(-1, 1, .01)
X <- expand.grid(X1, X2)
names(X) <- c("X1", "X2")
formula <- ~X1 + X2

test_that("Find the 2^2 factorial with opt-montecarlo", {
    expect_more_than(
        getEff(formula,
               optimalDesign(formula, X, 4, "D", 1000000),
               criteria = "D")$D,
        99)
    expect_more_than(
        getEff(formula,
               optimalDesign(formula, X, 4, "A", 1000000),
               criteria = "A")$A,
    73)
    expect_less_than(
        getEff(formula,
              optimalDesign(formula, X, 4, "I", 1000000),
              evaluation = X,
              criteria = "I")$I,
        2)
    expect_less_than(abs(
        getEff(formula,
           optimalDesign(formula, X, 4, "G", 50000),
               evaluation = X,
               criteria = "G")$G
        - 100), 2)
})

test_that("Find the 2^2 factorial with opt-genetic", {
    expect_more_than(
        getEff(formula,
               geneticdesign(formula, X, 4, "D", iter = 100000, 10),
               criteria = "D")$D,
        99)
})

    # expect_more_than(
    #     getEff(formula,
    #            optimalDesign(formula, X, 4, "A", 1000000),
    #            criteria = "A")$A,
    # 73)
    # expect_less_than(abs(
    #     getEff(formula,
    #           optimalDesign(formula, X, 4, "I", 1000000),
    #           evaluation = X,
    #           criteria = "I")$I
    #     - 1.06), .1)
    # expect_less_than(abs(
    #     getEff(formula,
    #            optimalDesign(formula, X, 4, "D", 1000000),
    #            evaluation = X,
    #            criteria = "G")$G
    #     - 100), .1)
# })
