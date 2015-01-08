#
# Design generation functions.
#


###########################################################################
##                           FACTORIAL DESIGN                            ##
###########################################################################


# user facing program (user doesn't need depth term)
genFactorial <- function(n_levels, n_factors)
{
    genFact_recursive(n_levels, n_factors, n_factors)
}

genFact_recursive <- function(n_levels, n_factors, depth)
{
    # center levels to zero and get scaling like in AlgDesign
    levels <- scale(1:n_levels - n_levels/2 -.5, scale=(n_levels%%2/2+.5))
    this <- rep(rep(levels, n_levels^(n_factors-depth)), rep(n_levels^(depth-1), n_levels^(n_factors-depth+1)))

    if (depth > 1)
        data.frame(unname(cbind(this, genFact_recursive(n_levels, n_factors, depth-1) )))
    else
        this
}

# Non-recursive method?
# genFactorial <- function(n_levels, n_factors)
# {
#     # get +/- end points
#     if (n_levels %% 2) # odd
#     {
#         E <- (n_levels - 1) / 2
#     }
#     else # even
#     {
#         E <- n_levels - 1
#     }
#     # get factorial points
#     X1 <- seq(-E, E, length.out = n_levels)
#     i <- 1
#     # repeat for number of variables
#     while (i < n_factors)
#     {
#         X1 <- cbind(X1, X1)
#         i <- i + 1
#         colnames(X1)[i] <- paste0("X", i)
#     }
#     expand.grid(X1)
# }


###########################################################################
##                       Central Composite Design                        ##
###########################################################################

genCCD <- function(n_factors, n_center, alpha)
{
    if (missing(alpha))
    {
        alpha <- sqrt(n_factors)
    }
    if (missing(n_center))
    {
        n_center <- 1
    }

    A <- as.matrix(genFactorial(2, n_factors))
    C <- matrix(0, nrow = n_center, ncol = n_factors)

    E <- CCDcol(n_factors, 1)

    if (n_factors > 1)
    {
        for (i in 2:n_factors)
        {
            E <- cbind(E, CCDcol(n_factors, i))
        }
    }
    E <- E * alpha

    as.data.frame(rbind(A, E, C))
}

CCDcol <- function(n_factors, col)
{
    column <- rep(0, n_factors * 2)
    column[2 * (col - 1) + c(1, 2)] <- c(1, -1)
    column
}


# genBBD <- function()
# {
# }

