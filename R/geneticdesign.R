geneticdesign <- function(formula, dataframe, n, criterion = "D", iter = 100000,
                          M = 10, evo = T)
{
    if (! (criterion %in% c("D", "A")))
    {
        stop("Only D criterion currently available")
    }

    # turn data frame into primes
    K_in <- length(names(dataframe))
    if (K_in > 10)
    {
        stop("Too many variables")
    }
    primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    newdata <- data.frame(ZZZJUNK = 1)
    for (i in 1:K_in)
    {
        newdata[[names(dataframe)[i]]] <- primes[i]
    }

    theformula <- model.matrix(formula, newdata)[-1]

    parents <- matrix(NA, nrow = n * K_in, ncol = M)
    for (i in 1:M)
    {
        parents[ , i] <- runif(n * K_in, -1, 1)
    }

    if (criterion == "D")
    {
        children <- opt_geneticrealcpp(parents, n, theformula, iter, (1:M) - 1, crit = 1, evo)
    } else if (criterion == "A")
    {
        children <- opt_geneticrealcpp(parents, n, theformula, iter, (1:M) - 1, crit = 2, evo)
    }

    return(children)



    # Cbest <- -5
    # for (i in 1:M)
    # {
    #     child <- data.frame(matrix(children[ , i], nrow = n, byrow = FALSE))
    #     names(child) <- names(dataframe)
    #     if (criterion == "D")
    #     {
    #         Cstar <- getEff(formula, child, criteria = "D")$D
    #     } else if (criterion == "A")
    #     {
    #         Cstar <- getEff(formula, child, criteria = "A")$A
    #     }

    #     if (Cstar > Cbest)
    #     {
    #         Cbest <- Cstar
    #         bestmodel <- child
    #     }
    # }
    # bestmodel
}
