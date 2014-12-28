optimalDesign <- function(formula, candidate, n, N = 100)
{

    # convert dataframe of design points to matrix of model points
    candidateX <- model.matrix(formula, candidate)

    # initial indices
    initial <- initialDesign(candidateX, n)
    # list of all indices
    all <- 1:nrow(candidateX)
    # initial <- initialDesignR(candidateX, n)

    M <- XprimeX(candidateX[initial, ])
    M_inv <- solve(M)

    current <- initial

    for (i in 1:N)
    {
        prop <- proposition(current, candidateX, all)
        if (moreEfficient(M_inv, prop$x, prop$y))
        {
            if ( (attempt <-
                  tryToSwitch(current, prop$new, prop$old, candidateX))$Y )
            {
                M_inv <- attempt$M_inv
                current <- attempt$current
            }
        }
    }
    candidateX[current, ]
}

# {{{ Small functions

DtoX <- function(D)
{
    cbind(rep(1, nrow(D)), D[ , 1], D[ , 2])
}

XprimeX <- function(X)
{
    t(X) %*% X
}

# }}} Small functions

moreEfficient <- function(M_inv, x, y)
{
    dxy <- x %*% M_inv %*% t(y)
    dx  <- x %*% M_inv %*% t(x)
    dy  <- y %*% M_inv %*% t(y)

    ((1 + dx) * (1 - dy) + dxy^2 - 1) > 0
}

proposition <- function(current, candidateX, all)
{
    old <- sample(current, 1)
    new <- sample(all[!(all %in% current)], 1)

    x <- DtoX(rbind(candidateX[new, ]))
    y <- DtoX(rbind(candidateX[old, ]))

    list(old = old, new = new, x = x, y = y)
}

tryToSwitch <- function(current, new, old, candidateX)
{
    current[current == old] <- new

    D <- candidateX[current, ]
    M <- XprimeX(DtoX(D))

    M_inv <- try(solve(M))

    if (class(M_inv) == "matrix")
    {
        return(list(current = current, M_inv = M_inv, Y = TRUE))
    } else {
        return(FALSE)
    }
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
