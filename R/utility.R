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
            spv <- N * apply(evaluation %*% M * evaluation, 1, sum)

            apply(N * t(candidate) * (infoMat%*%t(candidate)), 2, sum)
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

#
# what to make a human readable spv equation. maybe a 'formula' too
#
# do we really need to run this sieve every time? and can it be quicker?
# john didn't write this
#
sieve <- function(n)
{
    if(n > 1e8) stop("n too large")
    primes <- rep(TRUE, n)
    primes[1] <- FALSE
    last.prime <- 2L
    fsqr <- floor(sqrt(n))
    while (last.prime <= fsqr)
    {
        primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
        sel <- which(primes[(last.prime+1):(fsqr+1)])
        if(any(sel)) last.prime <- last.prime + min(sel)
        else last.prime <- fsqr+1
    }
    which(primes)
}

spvEquation <- function(formula, design)
{
    design <- model.matrix(formula, design)
    infoMat <- solve(t(design)%*%design)

    # would like sieve to return specified number of primes. not primes leq
    primes <- c(1, sieve(100))[1:ncol(design)]
    primey <- outer(primes, primes)
    terms <- unique(c(primey))
    coef <- cbind( rep(0, length(terms)), terms)

    # an apply statement?
    for(i in 1:length(terms)) coef[i,1] <- sum(infoMat*(primey == terms[i]))

    # will later want to get terms right and all (like turn 4 into 2^2 to
    # indicate x^2)
    
    return(coef)
}
