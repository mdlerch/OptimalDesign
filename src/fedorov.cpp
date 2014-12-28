#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec fedorovcpp(const arma::mat& xpx, IntegerVector current,
                         IntegerVector complete)
{
    // NOTE: this allows multiple copies of the same design point

    int Ncomplete = complete.size();
    int Ncurrent = current.size();

    int in = 0;
    int out = 0;

    while (in == out)
    {
        in = rand() % Ncomplete;
        out = rand() % Ncurrent;
    }

    arma::rowvec row_in = xpx.row(in);
    arma::rowvec row_out = xpx.row(out);

    arma::vec dii = row_in * xpx * row_in.t();
    arma::vec doo = row_out * xpx * row_out.t();
    arma::vec dio = row_in * xpx * row_out.t();

    return ((1 + dii) * (1 - doo) + dio * dio - 1);

}
