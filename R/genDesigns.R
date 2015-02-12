#
# Design generation functions.
#


###########################################################################
##                           FACTORIAL DESIGN                            ##
###########################################################################


# Non-recursive method
#
# This is so simple I think we should trash the recursive one
# I thought it was cute but fuck it
genFactorial <- function(n.levels, n.factors)
{
    # construct levels
    if (n.levels %% 2) # odd
    {
        levels <- 1:n.levels - mean(1:n.levels)
    }
    else # even
    {
        levels <- (1:n.levels - mean(1:n.levels)) * 2
    }

    repeatedLevels <- data.frame(replicate(n.factors,levels))

    expand.grid(repeatedLevels)
}

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

    CCDcol <- function(n_factors, col)
    {
        column <- rep(0, n_factors * 2)
        column[2 * (col - 1) + c(1, 2)] <- c(1, -1)
        column
    }

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

# I think we'll need this before we get to box-behnken:
###########################################################################
##                       Fractional-Factorial Design                     ##
###########################################################################

# genFFD <- function()
# {
# }

###########################################################################
##                       Bahcks-Baihynkin Design                         ##
###########################################################################

# t = number of factors per block
# r = number of blocks in which a factor appears
# lambda = number of times that each pair of factors appear in the same block
#
# By definition n.factors must be greater than 2
# probably should have a restriction on n.blocks given n.factors
genBBD <- function(n.factors, n.blocks, n.centerPoints, t=2, r, lambda)
{
    
}

###########################################################################
##                       Plackett-Burman Design                          ##
###########################################################################

# I don't know how many designs we want to put in here.
# There's an infobox at the bottom of the wiki page on some of these designs
# that has a list of what I presume are important designs.

# genPBD <- function()
# {
# {
