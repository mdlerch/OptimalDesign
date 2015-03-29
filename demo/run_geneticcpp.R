library(RcppArmadillo)

M <- 10
n <- 4
parents <- matrix(NA, nrow = n * 2, ncol = M)
for (i in 1:M)
{
    parents[ , i] <- runif(n * 2, -1, 1)
}

Rcpp::sourceCpp("./src/opt_genetic_real.cpp")

out <- opt_geneticrealcpp(parents, n, 10000, (1:M) - 1)

cbind(out[ , 1], parents[ , 1])
