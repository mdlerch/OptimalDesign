optimalDesign <- function(formula, candidate, n, iter = 100)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
#    set.seed(8675309)
    current <- initDesign(candidateX, n)
 
    M <- XprimeX(candidateX[current, ])
    M_inv <- solve(M)

    current <- fedorovcpp(M_inv, candidateX, current, 1:nrow(candidate), iter)

    candidate[current, ]
}

XprimeX <- function(X)
{
    t(X) %*% X
}