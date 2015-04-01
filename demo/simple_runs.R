library(OptimalDesign)

X <- genFactorial(8, 2) / 7
Xe <- genFactorial(50, 2) / 49

optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "A", repeated = F)
optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "D", repeated = T)
optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "D", repeated = F)

Ddesign <- optimalDesign(~X1+X2, candidate = X, n = 8, criterion = "D", repeated = T, iter = 100000)
Gdesign <- optimalDesign(~X1+X2, candidate = Xe, n = 8, criterion = "G", repeated = T, iter = 100000)
getEff(formula = ~X1 + X2, design = Ddesign, criteria = "G", evaluation = Xe)
getEff(formula = ~X1 + X2, design = Gdesign, criteria = "G", evaluation = Xe)
