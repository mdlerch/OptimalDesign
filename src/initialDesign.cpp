#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec initDesign(const arma::mat& X)
{
    // get dimensions of candidate matrix
    int nr = X.n_rows;
    int nc = X.n_cols;

    // svd decomposition (dc -- divide and conquer is fastest and default)
    // use econ because we only need left side
    // I believe econ also, by default, gives nr by nc rather than nr by nr
    arma::mat U, V;
    arma::vec S;
    arma::svd_econ(U, S, V, X, "left");

    arma::vec l2norm(nr);

    // Leverages are the root L^2-norm of U.
    // We don't need actual leverages, just order so don't bother rooting.
    int i = 0;
    do {
        // method when U is nr by nr
        // l2norm(i) = arma::norm((U.row(i)).cols(0, nc - 1), 2);
        l2norm(i) = arma::norm(U.row(i), 2);
        i++;
    } while(i < nr);

    return l2norm;
}
