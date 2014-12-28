optimalDesign <- function(formula, candidate, n, iter = 100)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
    current <- initialDesign(candidateX, n)
    # initial <- initialDesignR(candidateX, n)

    M <- XprimeX(candidateX[current, ])
    M_inv <- solve(M)

    current <- fedorovcpp(M_inv, candidateX, current, 1:nrow(candidate), iter)

    candidate[current, ]
}

XprimeX <- function(X)
{
    t(X) %*% X
}

#
# Returns the index of points in candidateX set with the largest leverages. This prevents
# singular initial designs for the fedorov algorithm
#
# The leverages (actually the square roots of leverages) are computed by taking the l2 norm 
# (aka Frobenius norm) of the left singular row vectors obtained from the SVD of the 
# candidate set matrix
#

initialDesignR <- function(candidateX, designSize)
{
    l2norm <- function(v) norm(as.matrix(v), 'F')
    leverages <- apply(svd(candidateX)$u, 1, l2norm) # actually square roots of leverages

    order(leverages, decreasing=T)[1:designSize]
}
