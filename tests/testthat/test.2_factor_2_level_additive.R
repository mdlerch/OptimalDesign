library(OptimalDesign)
context("Two level two factors additive model should be simple")

# candidate
x = seq(-1, 1, 1)
y = seq(-1, 1, 1)
Xsmall <- expand.grid(x, y)
names(Xsmall) <- c("x", "y")
x = seq(-1, 1, .1)
y = seq(-1, 1, .1)
Xlarge <- expand.grid(x, y)
names(Xlarge) <- c("x", "y")
X2squared <- data.frame(x = c(-1, 1, -1, 1), y = c(-1, -1, 1, 1))

# function to return the indices of original matrix, mat, that new matrix,
# matcher, matches
get_nums <- function(mat, matcher)
{
    matchrow <- function(r, mat)
    {
        intersect(which(r[1] == mat[ , 1]), which(r[2] == mat[ , 2]))
    }
    out <- apply(matcher, 1, function(x) matchrow(x, mat))
    sort(out)
}

test_that("2^2 small grid", {
              expect_equal(get_nums(Xsmall, X2squared),
                           sort(as.numeric(rownames(optimalDesign(~x+y, Xsmall, 4, criteria = "D")))))
})

test_that("2^2 large grid", {
              expect_equal(get_nums(Xlarge, X2squared),
                           sort(as.numeric(rownames( optimalDesign(~x+y, Xlarge, 4, criteria = "D")))))
})
