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
# rename these functions as you see fit (or the variables for that matter)
#

genFact <- function(n_levels, n_terms, i=n_terms)
{
    # center levels to zero and get scaling like in AlgDesign
    levels <- scale(1:n_levels - n_levels/2 -.5, scale=(n_levels%%2/2+.5))
    this <- rep(rep(levels ,n_levels^(n_terms-i)), rep(n_levels^(i-1), n_levels^(n_terms-i+1)))

    if(i>1) data.frame(unname(cbind( this, genFact(n_levels,n_terms,i-1) )))
    else this
}

genBBD <- function()
{

}

genCCD <- function()
{

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
