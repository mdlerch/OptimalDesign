library(OptimalDesign)

M <- 10
n <- 4
K_in <- 2
parents <- matrix(NA, nrow = n * K_in, ncol = M)
for (i in 1:M)
{
    parents[ , i] <- runif(n * K_in, -1, 1)
}

formula <- ~X1 * X2
design <- data.frame(X1 = 2, X2 = 3)

out <- geneticdesign(formula = formula, dataframe = design, n = 4, crit = "D")
geneticdesign(formula = formula, dataframe = design, n = 4, crit = "A")


for (i in 1:M)
{
    child <- data.frame(matrix(children[ , i], nrow = n, byrow = FALSE))
    names(child) <- names(dataframe)
    if (criterion == "D")
    {
        Cstar <- getEff(formula, child, criteria = "D")$D
    } else if (criterion == "A")
    {
        Cstar <- getEff(formula, child, criteria = "A")$A
    }

    if (Cstar > Cbest)
    {
        Cbest <- Cstar
        bestmodel <- child
    }
}
bestmodel
