#
# Collection of small utility functions
#

getEff <- function(formula, design, candidate)
{
    candidate <- model.matrix(formula, candidate/max(candidate))
    X <- model.matrix(formula, design/max(candidate))
    p <- ncol(X)
    N <- nrow(X)
    infoMat <- solve(t(X)%*%X)

    # evaluate efficiencies
    dEff <- 100 / (det(infoMat)^(1/p) * N)
    aEff <- 100*p / (sum(diag(infoMat)) * N)

    efficiencies <- list(dEff=dEff, aEff=aEff)

    # calculate iff candidate set provided
    # I don't know about these to be honest
    if(!missing(candidate))
    {
        spv <- apply(N * t(candidate) * (infoMat%*%t(candidate)), 2, sum)
        gEff <- 100 * sqrt(p/N) / max(spv^.5)
        iEff <- mean(spv^.5)

        efficiencies <- c(efficiencies, list(gEff=gEff, iEff=iEff))
    }

    efficiencies
}

#
# Generate factorial designs. It's recursive!
#

# user facing program (user doesn't need depth term)
genFactorial <- function(n_levels, n_terms)
{
    genFact_recursive(n_levels, n_terms, n_terms)
}

genFact_recursive <- function(n_levels, n_terms, depth)
{
    # center levels to zero and get scaling like in AlgDesign
    levels <- scale(1:n_levels - n_levels/2 -.5, scale=(n_levels%%2/2+.5))
    this <- rep(rep(levels, n_levels^(n_terms-depth)), rep(n_levels^(depth-1), n_levels^(n_terms-depth+1)))

    if (depth > 1)
        data.frame(unname(cbind(this, genFact_recursive(n_levels, n_terms, depth-1) )))
    else
        this
}

# Non-recursive method?
# genFactorial <- function(n_levels, n_terms)
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
#     while (i < n_terms)
#     {
#         X1 <- cbind(X1, X1)
#         i <- i + 1
#         colnames(X1)[i] <- paste0("X", i)
#     }
#     expand.grid(X1)
# }



genBBD <- function()
{

}

genCCD <- function()
{

}
