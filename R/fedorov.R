optimalDesign <- function(formula, candidate, n, iter = 100)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
    current <- initDesign(candidateX, n)

    current <- fedorovcpp(candidateX, current, 1:nrow(candidate), iter)

    # R indices start with 1
    candidate[current + 1, ]
}
