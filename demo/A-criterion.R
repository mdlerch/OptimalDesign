library(OptimalDesign)

Xbig <- genFactorial(10, 2) / 9
Xsmall <- genFactorial(3, 2)


Xopt <- optimalDesign(~X1 + X2, Xbig, 4, criterion = "A")

getEff(~X1 + X2, design = Xsmall[1:4, ], criteria = "A")
getEff(~X1 + X2, design = Xopt, criteria = "A")

Xopt
Xsmall[1:4, ]

