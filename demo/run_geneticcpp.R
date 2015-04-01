library(RcppArmadillo)

M <- 10
n <- 4
K_in <- 2
parents <- matrix(NA, nrow = n * K_in, ncol = M)
for (i in 1:M)
{
    parents[ , i] <- runif(n * K_in, -1, 1)
}

Rcpp::sourceCpp("./src/opt_genetic_real.cpp")

formula <- ~X1 * X2
design <- data.frame(X1 = 2, X2 = 3)
X <- model.matrix(formula, design)
X <- as.numeric(X)

out <- opt_geneticrealcpp(parents, n, c(2, 3, 6), 100000, (1:M) - 1)

cbind(out[ , 1], parents[ , 1])


library(OptimalDesign)


formula <- ~X1 * X2
design <- data.frame(X1 = 2, X2 = 3)

geneticdesign(formula = formula, dataframe = design, n = 4, criterion = "D")

