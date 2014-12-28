#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::ivec fedorovcpp(const arma::mat& xpxinv, const arma::mat& X,
                         arma::ivec current, arma::ivec complete, int iter)
{
    // NOTE: this allows multiple copies of the same design point

    int i = 0;

    while (i < iter)
    {
        // 1. Propose a index to put _in_ and an index to take _out_

        // design point number
        int in = 0;
        // element number of current vector
        int out = 0;
        // design point number
        int out_c = 0;

        int N = complete.n_rows;
        int n = current.n_rows;

        // Exclude non swap
        while (in == out_c)
        {
            in = rand() % N;
            out = rand() % n;
            // Note current is indexed in R as starting with 1,
            // but starting with 0 in C++
            out_c = current(out) - 1;
        }

        // 2. Get the 3 fedorov values based on the in and out
        arma::mat row_in = X.row(in);
        arma::mat row_out = X.row(out_c);

        arma::mat diim = row_in * xpxinv * row_in.t();
        arma::mat doom = row_out * xpxinv * row_out.t();
        arma::mat diom = row_in * xpxinv * row_out.t();

        double dii = diim(0, 0);
        double doo = doom(0, 0);
        double dio = diom(0, 0);

        // 3. Calculate the change in det( (X'X)^{-1} )
        double delta = ((1 + dii) * (1 - doo) + dio * dio - 1);

        // 4. If delta > 0, accept, else revert (ie do nothing)
        // TODO: do we need to verify that changes will produce an invertible
        // matrix?

        if (delta > 0)
        {
            // in is index in c++ (begins at 0)
            // add 1 to put in R indexing system
            current(out) = in + 1;
        }

        i++;
    }

    return current;

}
