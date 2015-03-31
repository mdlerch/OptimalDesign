library(RcppArmadillo)

M <- 10
n <- 9
K <- 4
parents <- matrix(NA, nrow = n * K, ncol = M)
for (i in 1:M)
{
    parents[ , i] <- runif(n * K, -1, 1)
}

Rcpp::sourceCpp("./src/opt_genetic_real.cpp")

formula <- ~X1 * X2
design <- data.frame(X1 = 2, X2 = 3)
X <- model.matrix(formula, design)
X <- as.numeric(X)

out <- opt_geneticrealcpp(parents, n, c(2, 3, 5, 7), 100000, (1:M) - 1)

cbind(out[ , 1], parents[ , 1])

