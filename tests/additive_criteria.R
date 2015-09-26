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
    msg <- paste("Test: > 99. Observed:", Deff)
    stop(paste(context, "D", msg))
}

optA <- optimalDesign(formula, X, 4, "A", 1000000)
if (Aeff <- getEff(formula, optA, criteria = "A")$A < 73)
{
    msg <- paste("Test: > 73. Observed:", Aeff)
    stop(paste(context, "A", msg))
}

optI <- optimalDesign(formula, X, 4, "I", 1000000)
if (Ieff <- getEff(formula, optI, criteria = "I", evaluation = X)$I > 3)
{
    msg <- paste("Test: < 3. Observed:", Ieff)
    stop(paste(context, "I", msg))
}

optG <- optimalDesign(formula, X, 4, "G", 50000)
if (Geff <- getEff(formula, optG, criteria = "G", evaluation = X)$G < 96)
{
    msg <- paste("Test: > 96. Observed:", Ieff)
    stop(paste(context, "G", msg))
}

context <- "2^2 factorial with genetic"

optD <- geneticdesign(formula, X, 4, "D", iter = 1000000, 10)
if (Deff <- getEff(formula, optD, criteria = "D")$D < 99)
{
    msg <- paste("Test > 99. Observed:", Deff)
    stop(paste(context, "D", msg))
}
