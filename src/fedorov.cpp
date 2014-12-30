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
    // index of X (complete candidate set) to insert
    int in;
    // index of X (complete candidate set) to remove
    int out;
    // index of current (to remove) out == current(out_c) - 1
    int out_c;

    // number of points in candidate set
    int N = complete.n_rows;
    // number of design points to use
    int n = current.n_rows;

    // swapping row vectors as 1 by C matrices
    arma::mat row_in(1, X.n_cols), row_out(1, X.n_cols);

    // Fedorov values as 1 by 1 matrices
    arma::mat diim(1, 1), doom(1, 1), diom(1, 1);
    // Fedorov values as doubles
    double dii, doo, dio;

    // Change in det( (X'X)^{-1})
    double delta;

    while (i < iter)
    {
        // 1. Propose a index to put _in_ and an index to take _out_
        // Don't pick two identical points
        do
        {
            in = rand() % N;
            out_c = rand() % n;
            out = current(out_c); // current now indexed from 0
        }
        while (in == out);

        // 2. Get the 3 Fedorov values based on the in and out
        // vectors to swap from design
        row_in = X.row(in);
        row_out = X.row(out);

        // Fedorov values as matrices
        diim = row_in * xpxinv * row_in.t();
        doom = row_out * xpxinv * row_out.t();
        diom = row_in * xpxinv * row_out.t();

        // Fedorov values as double
        dii = diim(0, 0);
        doo = doom(0, 0);
        dio = diom(0, 0);

        // 3. Calculate the change in det( (X'X)^{-1} )
        delta = ((1 + dii) * (1 - doo) + dio * dio - 1);

        // 4. If delta > 0, accept, else revert (ie do nothing)
        // TODO: do we need to verify that changes will produce an invertible
        // matrix?

        if (delta > 0)
        {
            // out_c is the index of current.  current(out_c) is an index of X
            // that is being removed in favor of the new "in" index of X.
            current(out_c) = in;
        }

        i++;
    }

    return current;

}
