#
# Collection of small utility functions
#

getEff <- function(formula, design, evaluation, criteria = c("D", "A", "I", "G"))
{
    X <- model.matrix(formula, design)
    p <- ncol(X)
    N <- nrow(X)
    infoMat <- t(X)%*%X

    # evaluate efficiencies
    D <- "D" %in% criteria
    A <- "A" %in% criteria
    I <- "I" %in% criteria
    G <- "G" %in% criteria
    if (D)
    {
        dEff <- 100 * det(infoMat)^(1/p) / N
    }
    if (A)
    {
        aEff <- 100*p / (sum(diag(M <- solve(infoMat))) * N)
    }
    if (I | G)
    {
        # calculate iff evaluation set provided
        # I don't know about these to be honest
        if(!missing(evaluation))
        {
            if (!exists("M"))
            {
                M <- solve(infoMat)
            }

            evaluation <- model.matrix(formula, evaluation)

            # Var(x \beta ) = x Var(\beta) x' = \sigma x (X'X)^{-1} x'
            # spv2 <- N * diag(evaluation %*% M %*% t(evaluation))
            # apply(N * t(evaluation) * (infoMat %*% t(evaluation)), 2, sum)
            spv <- N * apply(evaluation %*% M * evaluation, 1, sum)

            if (G)
            {
                gEff <- 100 * sqrt(p/N) / max(spv^.5)
            }
            if (I)
            {
                iEff <- mean(spv^.5)
            }
        }
        else
        {
            warning("Must give evaluation set for I and G efficiencies")
            if (I)
            {
                iEff <- NA
            }
            if (G)
            {
                gEff <- NA
            }
        }
    }

    efficiencies <- list()
    if (D)
    {
        efficiencies$D <- dEff
    }
    if (A)
    {
        efficiencies$A <- aEff
    }
    if (I)
    {
        efficiencies$I <- iEff
    }
    if (G)
    {
        efficiencies$G <- gEff
    }

    efficiencies
}
