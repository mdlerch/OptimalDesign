optimalDesign <- function(candidate, initial, n = length(initial), N = 100)
{

    M <- XtoM(DtoX(candidate[initial, ]))
    M_inv <- solve(M)
    all <- 1:nrow(candidate)

    current <- initial

    for (i in 1:N)
    {
        prop <- proposition(current, candidate, all)
        if (moreEfficient(M_inv, prop$x, prop$y))
        {
            if ( (attempt <- 
                  tryToSwitch(current, prop$new, prop$old, candidate))$Y )
            {
                M_inv <- attempt$M_inv
                current <- attempt$current
            }
        }
    }
    candidate[current, ]
}

# {{{ Small functions

DtoX <- function(D)
{
    cbind(rep(1, nrow(D)), D[ , 1], D[ , 2])
}

XtoM <- function(X)
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

proposition <- function(current, candidate, all)
{
    old <- sample(current, 1)
    new <- sample(all[!(all %in% current)], 1)

    x <- DtoX(rbind(candidate[new, ]))
    y <- DtoX(rbind(candidate[old, ]))

    list(old = old, new = new, x = x, y = y)
}

tryToSwitch <- function(current, new, old, candidate)
{
    current[current == old] <- new

    D <- candidate[current, ]
    M <- XtoM(DtoX(D))

    M_inv <- try(solve(M))

    if (class(M_inv) == "matrix")
    {
        return(list(current = current, M_inv = M_inv, Y = TRUE))
    } else {
        return(FALSE)
    }
}
