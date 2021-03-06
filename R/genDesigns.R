#
# Design generation functions.
#


###########################################################################
##                           FACTORIAL DESIGN                            ##
###########################################################################


genFactorial <- function(n_factors, n_levels)
{
    # construct levels
    if (n_levels %% 2) # odd
    {
        levels <- 1:n_levels - mean(1:n_levels)
    }
    else # even
    {
        levels <- (1:n_levels - mean(1:n_levels)) * 2
    }

    repeatedLevels <- data.frame(replicate(n_factors,levels))

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

    A <- as.matrix(genFactorial(n_factors, 2))
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

###########################################################################
##                       Fractional-Factorial Design                     ##
###########################################################################

genLatin <- function(size, iterations)
{
    design = opt_genLatin(size, iterations)
    noquote(matrix(LETTERS[design + 1], size, size))
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
#genBBD <- function(n.factors, n.blocks, n.centerPoints, t=2, r, lambda)
#{
#    
#}

###########################################################################
##                       Plackett-Burman Design                          ##
###########################################################################

# I don't know how many designs we want to put in here.
# There's an infobox at the bottom of the wiki page on some of these designs
# that has a list of what I presume are important designs.

# genPBD <- function()
# {
# {
