optimalDesign <- function(formula, candidate, n, criterion = "D", iter = 10000,
                          repeated = TRUE)
{

    if (repeated == FALSE && n > nrow(candidate))
    {
        stop("Requested design larger than candidate set and repeated trials not set")
    }

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
    current <- initDesign(candidateX, n, repeated)

    criteria <- c("D", "A", "G", "I")
    criterion <- match.arg(criterion, criteria)

    # initial set of legal swaps depends on repeated option
    if (repeated)
    {
        candidateidx <- (0:(nrow(candidate) - 1))
    }
    else
    {
        candidateidx <- NULL
        for (i in 0:(nrow(candidate) - 1))
        {
            if (!(i %in% current))
            {
                candidateidx <- c(candidateidx, i)
            }
        }
    }

    x=2
    
    if (criterion == "D")
    {
        current <- opt_montecarlocpp(candidateX, current, candidateidx, 1, iter, repeated)
    }
    else if (criterion == "A")
    {
        current <- opt_montecarlocpp(candidateX, current, candidateidx, 2, iter, repeated)
    }
    else if (criterion == "G")
    {
        current <- opt_montecarlocpp(candidateX, current, candidateidx, 3, iter, repeated)
    }
    else if (criterion == "I")
    {
        current <- fedorovcpp(candidateX, current, candidateidx, 4, iter, repeated)
    }
    else if (criterion == "I")
    {
        current <- fedorovcpp(candidateX, current, candidateidx, 4, iter, repeated)
    }
    else
    {
        stop("Only D, A, G, and I criteria are currently supported")
    }

    # R indices start with 1
    candidate[current + 1, ]
}
