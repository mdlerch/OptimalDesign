library(OptimalDesign)

X <- genFactorial(8, 2) / 7

optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "A", repeated = F)
optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "D", repeated = T)
optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "D", repeated = F)

