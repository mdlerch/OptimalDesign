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
arma::uvec opt_montecarlocpp(const arma::mat& Xc, arma::uvec current,
                             arma::uvec candidateidx, const arma::mat& Xe,
                             int crit, int iterations,
                             bool repeated)
{
    // Xc is the X matrix for all candidate points.
    // current is a vector of indexes to Xc to indicate which points are in
    //   the current design
    // candidateidx is a vector of indexes that are legal propositions for
    //   swapping with an element of current (all if repeated, all\current if
    //   not)
    // crit is the criteria of interest (1==D, 2==A, 3==G, 4==I)
    // iterations is the total number of iterations in the process
    // repeated is whether or not to allow candidate points to appear more than
    //   once in current.

    // index of candidateidx (legal set of indices) to insert
    int in_c;
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

    /**** Initial object(s) for i-criterion ****/
    arma::mat B = (Xe.t() * Xe) / N;

    /**** Initial object(s) for g-criterion ****/
    // objects for SVD of candidate set
    arma::mat U_can, V_can, D_can;
    arma::vec s_can;
    arma::svd_econ(U_can, s_can, V_can, Xe, "left");

    // objects for SVD of current design
    arma::mat U, V;
    arma::vec s;
    arma::svd_econ(U, s, V, X, "right");
    s = 1/s;
    arma::mat Dinv = arma::diagmat(s);

    // special formula for calculating g-criterion for current design
    arma::mat svd_thing = U_can * V * Dinv;
    svd_thing = svd_thing % svd_thing;
    arma::vec leverages = arma::sum(svd_thing, 1);
    double g_crit = leverages.max();

    int iter = 0;
    while (iter < iterations)
    {
        // 1. Propose an index to put _in_ and an index to take _out_
        // Don't pick two identical points
        do
        {
            in_c = rand() % N;
            in = candidateidx(in_c);
            out_c = rand() % n;
            out = current(out_c);
        }
        while (in == out);

        // 2. Calculate the change in the criterion
        if (crit == 1) // Criteria D
        {
            delta = get_delta_d(xpxinv, Xc.row(in), Xc.row(out));
        }
        else if (crit == 2) // Criteria A
        {
            delta = get_delta_a(xpxinv, Xc.row(in), Xc.row(out));
        }
        else if (crit == 3) // Criteria G
        {
            delta = get_delta_g(g_crit, Xc.rows(current), U_can);
        }
        else if (crit == 4) // Criteria I
        {
            delta = get_delta_i(xpxinv, Xc.row(in), Xc.row(out), B);
        }

        // 3. If delta > 0, accept, else revert (ie do nothing)
        if (delta > 0)
        {
            // out_c is the index of current.  current(out_c) is an index of Xc
            // that is being removed in favor of the new "in" index of Xc.
            current(out_c) = in;

            X = Xc.rows(current);

            // if not allowing repeated, update the legal candidates
            // this should be just swapping the one I am putting in with
            // the one I just pulled out
            if (!repeated)
            {
                candidateidx(in_c) = out;
            }

            // if new design not invertible, switch back
            // john found the 1e15 bound on stack overflow. might work
            if (cond(X.t() * X) > 1e15)
            {
                current(out_c) = out;
                if (!repeated)
                {
                    candidateidx(in_c) = in;
                }
            }
            else
            {
                xpxinv = inv(X.t() * X);

                if (crit == 3)
                {
                    g_crit = g_crit + delta;
                }
            }
        }

        iter++;
    }

    return current;
}
