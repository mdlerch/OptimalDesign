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
    expect_less_than(abs(
        getEff(formula1, design = X1, evaluation = Xe, criteria = "I")$I
        - 0.64), .1)
    expect_equal(
        getEff(formula1, design = X1, evaluation = Xe, criteria = "G")$G,
        100)
})

Xccd <- data.frame(X1 = c(rep(-1, 5), rep(0, 8), rep(1, 5)),
                   X2 = c(-1, -1, 0, 1, 1, -1, rep(0, 6), 1, -1, -1, 0, 1, 1),
                   X3 = c(-1, 1, 0, -1, 1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 1, 0, -1, 1))
Xe <- genFactorial(n.levels = 100, n.factors = 3) / 99

fccd <- ~X1*X2 + X1*X3 + X2*X3 + I(X1^2) + I(X2^2) + I(X3^2)

test_that("K=3, 4 center point CCD second order model (Jobo pg 195)", {
    expect_less_than(abs(
        getEff(fccd, design = Xccd, criteria = "D")$D
        - 39.6635), .1)
    expect_less_than(abs(
        getEff(fccd, design = Xccd, criteria = "A")$A
        - 28.6826), .1)
    expect_less_than(abs(
        getEff(fccd, design = Xccd, evaluation = Xe, criteria = "I")$I
        - 0.5880), .1)
    expect_less_than(abs(
        getEff(fccd, design = Xccd, evaluation = Xe, criteria = "G")$G
        - 83.6451), .1)
})
