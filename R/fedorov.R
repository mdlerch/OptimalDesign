optimalDesign <- function(formula, candidate, n, iter = 100)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
#    set.seed(8675309)
    current <- initDesign(candidateX, n)

    # R indices start with 1
    M <- XprimeX(candidateX[current + 1, ])
    M_inv <- solve(M)

    current <- fedorovcpp(M_inv, candidateX, current, 1:nrow(candidate), iter)

    # R indices start with 1
    candidate[current + 1, ]
}

XprimeX <- function(X)
{
    t(X) %*% X
}
