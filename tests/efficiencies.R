library(OptimalDesign)

context <- "2^2 factorial efficiency calculations"

# Make a dataset
X1 <- genFactorial(n_factors = 2, n_levels = 2)
formula1 <- ~X1 + X2

# Large evaluation set
Xe <- genFactorial(n_factors = 2, n_levels = 100) / 99

efficiencies <- getEff(formula1, design = X1, evaluation = Xe)

if (!all.equal(efficiencies$D, 100))
{
    stop(paste(context, "D"))
}
if (!all.equal(efficiencies$A, 100))
{
    stop(paste(context, "A"))
}
if (efficiencies$I < .60 | efficiencies$I > .65)
{
    stop(paste(context, "I"))
}
if (!all.equal(efficiencies$G, 100))
{
    stop(paste(context, "G"))
}

# Jobo pg 195
context <- "K=3, 4 center points, CCD efficiencies"

Xccd <- data.frame(X1 = c(rep(-1, 5), rep(0, 8), rep(1, 5)),
                   X2 = c(-1, -1, 0, 1, 1, -1, rep(0, 6), 1, -1, -1, 0, 1, 1),
                   X3 = c(-1, 1, 0, -1, 1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 1, 0, -1, 1))
Xe <- genFactorial(n_factors = 3, n_levels = 100) / 99

fccd <- ~X1*X2 + X1*X3 + X2*X3 + I(X1^2) + I(X2^2) + I(X3^2)

efficiencies <- getEff(fccd, design = Xccd, evaluation = Xe)

if (efficiencies$D < 39.5 | efficiencies$D > 39.8 )
{
    stop(paste(context, "D"))
}
if (efficiencies$A < 28.5 | efficiencies$A > 28.8)
{
    stop(paste(context, "A"))
}
if (efficiencies$I < .55 | efficiencies$I > .58)
{
    stop(paste(context, "I"))
}
if (efficiencies$G < 83.5 | efficiencies$G > 83.8)
{
    stop(paste(context, "G"))
}
