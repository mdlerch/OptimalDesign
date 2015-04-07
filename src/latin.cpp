#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

double get_delta_g(double, arma::mat, arma::mat);
double get_delta_d(arma::mat, arma::mat, arma::mat);
double get_delta_a(arma::mat, arma::mat, arma::mat);
double get_delta_i(arma::mat, arma::mat, arma::mat, arma::mat);
arma::vec delta_common(arma::mat, arma::mat, arma::mat);

// [[Rcpp::export]]
arma::uvec latin(int n, int iterations)
{

  // cells that contain zeros and ones
  arma::mat zeros;
  arma::mat ones;

  // index of Xc (complete candidate set) to insert
    int in;
    // index of current (to remove) out == current(out_c) - 1
    int out_c;
    // index of Xc (complete candidate set) to remove
    int out;

    // number of points in candidate set
    int N = candidateidx.n_rows;
    // number of design points to use
    int n = current.n_rows;

    // Swapping condition delta > 0 == good
    double delta;

    // Get (X'X)^{-1}
    arma::mat xpxinv;
    arma::mat X = Xc.rows(current);
    arma::inv(xpxinv, X.t() * X);
