#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec fedorovcpp(const arma::mat& Xc, arma::uvec current, int crit,
                      arma::ivec complete, int iter)
{
    // NOTE: this allows multiple copies of the same design point

    int i = 0;
    // index of Xc (complete candidate set) to insert
    int in;
    // index of Xc (complete candidate set) to remove
    int out;
    // index of current (to remove) out == current(out_c) - 1
    int out_c;

    // number of points in candidate set
    int N = complete.n_rows;
    // number of design points to use
    int n = current.n_rows;

    // swapping row vectors as 1 by C matrices
    arma::mat row_in(1, Xc.n_cols), row_out(1, Xc.n_cols);

    // Fedorov values as 1 by 1 matrices
    arma::mat diim(1, 1), doom(1, 1), diom(1, 1);
    // Fedorov values as doubles
    double dii, doo, dio;

    // Fedorov phi values
    arma::mat phiii(1, 1), phioo(1, 1), phiio(1, 1);

    // Change in det( (X'X)^{-1})
    double delta_d;
    double delta;

    // Get (X'X)^{-1}
    arma::mat xpxinv;
    arma::mat X = Xc.rows(current);
    arma::inv(xpxinv, X.t() * X);

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
        row_in = Xc.row(in);
        row_out = Xc.row(out);

        // Fedorov values as matrices
        diim = row_in * xpxinv * row_in.t();
        doom = row_out * xpxinv * row_out.t();
        diom = row_in * xpxinv * row_out.t();

        // Fedorov values as double
        dii = diim(0, 0);
        doo = doom(0, 0);
        dio = diom(0, 0);

        // 3. Calculate the change in det( (X'X)^{-1} )
        delta_d = ((1 + dii) * (1 - doo) + dio * dio - 1);

        if (crit == 1) // Criteria D
        {
            delta = delta_d;
        }
        else if (crit == 2) // Criteria A
        {
        }

        // 4. If delta_d > 0, accept, else revert (ie do nothing)
        if (delta_d > 0)
        {
            // out_c is the index of current.  current(out_c) is an index of Xc
            // that is being removed in favor of the new "in" index of Xc.
            current(out_c) = in;

            // if new design not invertible switch back
            X = Xc.rows(current);
            if (! arma::inv(xpxinv, X.t() * X))
            {
                current(out_c) = out;
            }
        }

        i++;
    }

    return current;

}
