library(OptimalDesign)

context <- "2^2 factorial with monte carlo"

# Make a dataset
X1 <- seq(-1, 1, .01)
X2 <- seq(-1, 1, .01)
X <- expand.grid(X1, X2)
names(X) <- c("X1", "X2")
formula <- ~X1 + X2

optD <- optimalDesign(formula, X, 4, "D", 1000000)
if (getEff(formula, optD, criteria = "D")$D < 99)
{
    stop(paste(context, "D"))
}

optA <- optimalDesign(formula, X, 4, "A", 1000000)
if (getEff(formula, optA, criteria = "A")$A < 73)
{
    stop(paste(context, "A"))
}

optI <- optimalDesign(formula, X, 4, "I", 1000000)
if (getEff(formula, optI, criteria = "I", evaluation = X)$I > 3)
{
    stop(paste(context, "I"))
}

optG <- optimalDesign(formula, X, 4, "G", 50000)
if (getEff(formula, optG, criteria = "G", evaluation = X)$G < 96)
{
    stop(paste(context, "G"))
}

context <- "2^2 factorial with genetic"

optD <- geneticdesign(formula, X, 4, "D", iter = 1000000, 10)
if (getEff(formula, optD, criteria = "D")$D < 99)
{
    stop(paste(context, "D"))
}
