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

#
# want to make a human readable spv equation. maybe a 'formula' too
#
# do we really need to run this sieve every time? and can it be quicker?
# john didn't write this. he prob ripped it off stack exchange
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

#
# Prime number decomposition function
#

primeDecomp <- function(x)
{
    div <- 2:x
    factors <- div[x %% div == 0]
    
    primeNums = c()

    # get the prime numbers you want
    while(length(factors) > 1 & !factors[length(factors)]%%1)
    {
        nextPrime = min(factors[!factors%%1])
        primeNums = c(primeNums, nextPrime)
        factors = (factors/nextPrime)[-which(factors==nextPrime)]		
    }

    if(length(factors)==1 & !factors[length(factors)]%%1) primeNums = c(primeNums, factors)

    # figure out those prime numbers' multiplicities
    # primes in first row, multiplicities is second
    as.data.frame(table(primeNums))
}

# This function calculates the SPV equation by a little trick keeping terms
# straight. Basically there is a tiny computer algebra system that represents
# variables and their products by prime numbers and their products
# respectively. It's the fund. theorem of arithmetic in action.
#
# right now this assumes there's an intercept term in in the model

spvEquation <- function(formula, design)
{
    # make model matrix and then info matirx
    modelMat <- model.matrix(formula, design)
    infoMat <- solve(t(modelMat) %*% modelMat)

    # would like sieve to return specified number of primes. not primes leq
    # prime numbers represent variables
    primes <- sieve(100)[1:ncol(design)]

    # get products of prime numbers corresponding to products of variables
    modelTerms = model.matrix(formula, data.frame(t(primes)))

    # outerproduct of vector of prime number products
    primeMat <- outer(as.numeric(modelTerms),as.numeric(modelTerms))

    # terms in spv equation as represented by products of primes
    terms <- unique(c(primeMat))

    # for each prime number product, p, found in the terms vector, sum up 
    # the elements in the info matrix that have indices equal to those in
    # the outer product matrix that have the prime number product p
    coef <- cbind( rep(0, length(terms)), terms)
    for(i in 1:length(terms)) coef[i,1] <- sum(infoMat*(primeMat == coef[i,2]))

    #########################################################
    # Now to make the spv function human readable
    #########################################################

    # compute prime number decompositions for prime number products
    coef.primed <- lapply(coef[,2], primeDecomp)

    # create data.frame for primes and the variable names to be able to index
    variable.names <- data.frame(primes, names(design))


    # initialize string to print out SPV
    spv.string <- as.character(coef[1])

    # perhaps should treat intercept differently in some of these data
    # structures since, in a way, it is a little different
    #
    # skips intercept term because added at initialization step above
    for(i in setdiff(which(coef[,1] != 0), 1))
  	{
        spv.string <- paste(spv.string, "+", coef[i])

        for(j in 1:nrow(coef.primed[[i]]))
        {
            prime <- coef.primed[[i]]$primeNums[j]
            mult <- coef.primed[[i]]$Freq[j]
            single.prime.string <- paste(variable.names[variable.names[,1]==prime, 2], "^", mult, sep="")

            spv.string <- paste(spv.string, "*", single.prime.string, sep="")
        }
    }

    spv.string
}
