optimalDesign <- function(formula, candidate, n, criterion = "D", iter = 10000)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
    current <- initDesign(candidateX, n)

    criteria <- c("D", "A", "G", "IV")
    criterion <- match.arg(criterion, criteria)

    if (criterion == "D")
    {
        current <- fedorovcpp(candidateX, current, 1:nrow(candidate), 1, iter)
    }
    else if (criterion == "A")
    {
        current <- fedorovcpp(candidateX, current, 1:nrow(candidate), 2, iter)
    } else
    {
        stop("Only D and A criteria are currently supported")
    }

    # R indices start with 1
    candidate[current + 1, ]
}
