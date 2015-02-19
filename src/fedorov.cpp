#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec fedorovcpp(const arma::mat& Xc, arma::uvec current,
                      arma::ivec complete, int crit, int iter)
{
    // NOTE: this allows multiple copies of the same design point
    // do we need the complete ivec as an argument? instead of complete.n_rows
    // on line 26 below, we could Xc.n_rows
    //
    // even if we can't get n_rows from Xc cause it's got the & address thing
    // or it's a constant, we could just pass a numeric instead of complete

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

    // Fedorov phi values as 1 by 1 matrices
    arma::mat phiiim(1, 1), phioom(1, 1), phiiom(1, 1), phioim(1, 1);
    // Fedorov phi values as doubles
    double phiii, phioo, phiio, phioi;

    // Change in det( (X'X)^{-1})
    double delta_d;
    double delta;

    // Get (X'X)^{-1}
    arma::mat xpxinv;
    arma::mat X = Xc.rows(current);
    arma::inv(xpxinv, X.t() * X);

    /**** SVD matrices for g-criterion ****/
    // objects for SVD of candidate set
    arma::mat U_can;
    arma::mat V_can;
    arma::mat D_can;
    arma::vec s_can;
    arma::svd_econ(U_can, s_can, V_can, Xc, "left");

    // objects for SVD of current design
    arma::mat U;
    arma::mat V;
    arma::vec s;
    arma::svd_econ(U, s, V, X, "right");
    s = 1/s;
    arma::mat Dinv = arma::diagmat(s);

    // special formula for calculating g-criterion for current design
    arma::mat svd_thing = U_can * V * Dinv;
    svd_thing = svd_thing % svd_thing;
    arma::vec leverages = arma::sum(svd_thing, 1);
    double g_crit = leverages.max();

    // objects for SVD of design with potential new point
    arma::uvec current_test;
    arma::mat X_test;
    arma::mat U_test;
    arma::mat V_test;
    arma::vec s_test;
    arma::mat Dinv_test;
    arma::mat svd_thing_test;
    arma::vec leverages_test;
    double g_crit_test;

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
            phiiim = row_in * xpxinv * xpxinv * row_in.t();
            phiiom = row_out * xpxinv * xpxinv * row_in.t();
            phioim = row_in * xpxinv * xpxinv * row_out.t();
            phioom = row_out * xpxinv * xpxinv * row_out.t();

            // convert phi values to double
            phiii = phiiim(0, 0);
            phiio = phiiom(0, 0);
            phioi = phioim(0, 0);
            phioo = phioom(0, 0);

            delta = ( (1 - dii) * phioo + dio * (phiio + phioi) - (1 + doo) * phiii ) / (1 + delta_d);
        }
        else if (crit == 3) // Criteria G
        {
            // calculate SVD for potential new design
            current_test = current;
            current_test(out_c) = in;
            X_test = Xc.rows(current_test);
            arma::svd_econ(U_test, s_test, V_test, X_test, "right");
            s_test = 1/s_test;
            Dinv_test = arma::diagmat(s_test);

            // special formula for g-criterion for potential new design
            svd_thing_test = U_can * V_test * Dinv_test;
            svd_thing_test = svd_thing_test % svd_thing_test;
            leverages_test = arma::sum(svd_thing_test, 1);
            g_crit_test = leverages_test.max();
            
            // if test g-criterion higher, delta > 0 indicates success
            delta = g_crit_test - g_crit;
        }

        // 4. If delta > 0, accept, else revert (ie do nothing)
        if (delta > 0)
        {   
            // out_c is the index of current.  current(out_c) is an index of Xc
            // that is being removed in favor of the new "in" index of Xc.
            current(out_c) = in;
            X = Xc.rows(current);

            if (crit == 3)
            {
                g_crit = g_crit_test;
            }

            // if new design not invertible, switch back
            if (! arma::inv(xpxinv, X.t() * X))
            {
                current(out_c) = out;

                if (crit == 3)
                {
                    g_crit = g_crit - delta;
                }
            }
        }

        i++;
    }

    return current;

}
