require(AlgDesign)

X <- as.matrix(gen.factorial(levels = 3, nVars = 2))

optimalDesign(X, 1:4)


all <- 1:8
current <- 1:4

sample(all[!(all %in% current)], 1)


DD <- X[1:4, ]
XX <- DtoX(DD)
M <- t(XX) %*% XX

M_inv <- solve(M)
