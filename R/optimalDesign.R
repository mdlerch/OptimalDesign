optimalDesign <- function(formula, candidate, n, criteria = "D", iter = 10000)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
    current <- initDesign(candidateX, n)

    criteriums <- c("D", "A", "G", "IV")
    criteria <- match.arg(criteria, criteriums)

    if (criteria == "D")
    {
        current <- fedorovcpp(candidateX, current, 1:nrow(candidate), iter)
    } else {
        stop("Only D criteria currently supported")
    }

    # R indices start with 1
    candidate[current + 1, ]
}
