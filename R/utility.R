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
