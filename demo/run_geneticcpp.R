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

out <- geneticdesign(formula = formula, dataframe = design, n = 4, crit = "D", evo = FALSE)

library(animation)
ani.options(interval=.1)
ani.options(autobrowse=FALSE)
ani.options(autoplay=FALSE)
saveGIF({
for (i in 1:100)
    {
        plot(0,0, type = 'n', xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
        points(out[1:4, 1, i], out[5:8, 1, i], pch = 19, col = c(1, 2, 3, 4), cex = 2)
    }
})


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
