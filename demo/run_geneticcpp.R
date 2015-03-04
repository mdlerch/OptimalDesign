library(RcppArmadillo)

Rcpp::sourceCpp("./src/genetic.cpp")

M <- 10
n <- 4
parents <- matrix(NA, nrow = n * 2, ncol = 10)
for (i in 1:M)
{
    parents[ , i] <- runif(n * 2, -1, 1)
}

children <- geneticcpp(parents, n, 1, (1:M) - 1)

cbind(children[ , 1], parents[ , 1])
