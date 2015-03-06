#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;

arma::ivec primedecomp(int);
double get_delta_g(double, arma::mat, arma::mat);
double get_delta_d(arma::mat, arma::mat, arma::mat);
double get_delta_a(arma::mat, arma::mat, arma::mat);
arma::vec delta_common(arma::mat, arma::mat, arma::mat);
void grblend(arma::mat &, uint, arma::mat, arma::uvec, double);
void grcreep(arma::mat &, uint, double size, double);
void grmutat(arma::mat &, uint, double);

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat opt_geneticrealcpp(arma::mat parents, int n, int iterations, arma::uvec pidx)
{
    // parents is a matrix with each row being a vector of the design points
    //  with the first n being the first variable, the second n being the second
    //  varaible etc.
    // n is the design size
    // iterations is the number of iterations of the genetic algorithm
    // pidx is just 1:M where M is number of parents; easier to pass than write
    //  in c++!

    int M = parents.n_cols;
    arma::uvec second_parent;
    arma::cube xpxinv(1 + parents.n_rows / n, 1 + parents.n_rows / n, M);
    arma::mat X(n, 3), xpxi(3, 3);
    arma::vec Xv;
    arma::mat rowin(1, 3);

    arma::mat children = parents;

    int iter, child;

    int i;
    double delta;

    // hard code ~X1 + X2 for now.  Later use primedecomp to get actual design
    // hard code D criterion for now add other criteria later.

    // make x'x.inv

    int j;

    for (i = 0; i < M; ++i)
    {
        Xv = parents.col(i);
        for (j = 0; j < n; j++)
        {
            X(j, 0) = 1;
        }
        X.col(1) = Xv.subvec(0, n - 1);
        X.col(2) = Xv.subvec(n, 2 * n -1);
        arma::inv(xpxi, X.t() * X);
        xpxinv.slice(i) = xpxi;
    }


    iter = 0;
    while(iter < iterations)
    {
        second_parent = RcppArmadillo::sample(pidx, M, false);

        children = parents;


        for (child=0; child<M; ++child)
        {
            grblend(children, child, parents, second_parent, 0);
            grcreep(children, child, .3, 0);
            grmutat(children, child, .5);

            for (j = 0; j < n; ++j)
            {
                X(j, 0) = 1;
            }


            X.col(1) = parents.col(child).subvec(0, n - 1);
            X.col(2) = parents.col(child).subvec(n, 2 * n - 1);

            for (i = 0; i < n; ++i)
            {
                delta = 0;
                if (X(i, 1) != children(i, child) | X(i, 2) != children(i + n, child))
                {
                    rowin(0, 0) = 1;
                    rowin(0, 1) = children(i, child);
                    rowin(0, 2) = children(i + n, child);
                    delta += get_delta_d(xpxinv.slice(child), rowin, X.row(i));
                    X.row(i) = rowin;
                }
            }
            if (delta > 0)
            {

                if (arma::inv(xpxi, X.t() * X))
                {
                    parents.col(child) = children.col(child);
                    xpxinv.slice(child);
                }
            }
        }
        iter++;
    }

    return parents;
}

void grblend(arma::mat & children, uint child, arma::mat parents, arma::uvec parent2, double alpha)
{
    int i;

    for (i=0; i<children.n_rows; ++i)
    {
        if (R::runif(0, 1) < alpha)
        {
            children(i, child) = parents(i, parent2[i]);
        }
    }
}

void grcreep(arma::mat & children, uint child, double size, double alpha)
{
    int i;

    for (i=0; i<children.n_rows; ++i)
    {
        if (R::runif(0, 1) < alpha)
        {
            children(i, child) = children(i, child) + R::rnorm(0, size);
            if (children(i, child) > 1)
            {
                children(i, child) = 1;
            }
            else if (children(i, child) < -1)
            {
                children(i, child) = -1;
            }
        }
    }
}

void grmutat(arma::mat & children, uint child, double alpha)
{
    int i;

    for (i=0; i<children.n_rows; ++i)
    {
        if (R::runif(0, 1) < alpha)
        {
            children(i, child) = R::runif(-1, 1);
        }
    }
}

arma::ivec primedecomp(int num)
{
    arma::ivec primes(num);

    int i = 2;
    int idx = 0;

    while (i * i <= num)
    {
        if (num % i == 0)
        {
            primes(idx) = i;
            num = num / i;
            ++idx;
        }
        else
        {
            ++i;
        }
    }
    if (num > 1)
    {
        primes(idx) = num;
        ++idx;
    }

    primes.resize(idx);
    return primes;
}

arma::vec delta_common(arma::mat xpxinv, arma::mat row_in, arma::mat row_out)
{
    // Fedorov values as 1 by 1 matrices
    arma::mat diim(1, 1), doom(1, 1), diom(1, 1);
    // Fedorov values as doubles
    double dii, doo, dio;

    double delta;

    arma::vec common(4);


    diim = row_in * xpxinv * row_in.t();
    doom = row_out * xpxinv * row_out.t();
    diom = row_in * xpxinv * row_out.t();


    // Fedorov values as double
    dii = diim(0, 0);
    doo = doom(0, 0);
    dio = diom(0, 0);

    delta = ((1 + dii) * (1 - doo) + dio * dio - 1);


    common(0) = delta;
    common(1) = dii;
    common(2) = doo;
    common(3) = dio;

    return common;
}

double get_delta_d(arma::mat xpxinv, arma::mat row_in, arma::mat row_out)
{
    arma::vec common;

    common = delta_common(xpxinv, row_in, row_out);

    return common(0);
}

double get_delta_a(arma::mat xpxinv, arma::mat row_in, arma::mat row_out)
{
    // fedorov A phi values as matrices
    arma::mat phiiim, phioom, phiiom, phioim;
    // fedorov A phi values as double
    double phiii, phioo, phiio, phioi;

    double delta;

    arma::vec common;

    common = delta_common(xpxinv, row_in, row_out);

    double delta_d = common(0);
    double dii = common(1);
    double dio = common(2);
    double doo = common(3);

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

    return delta;
}

double get_delta_g(double g_crit_old, arma::mat X, arma::mat U_can)
{
    double delta;
    arma::mat Dinv, svd_thing;
    arma::mat U, V;

    arma::vec s, leverages;

    double g_crit;

    arma::svd_econ(U, s, V, X, "right");

    s = 1 / s;

    Dinv= arma::diagmat(s);

    svd_thing = U_can * V * Dinv;
    svd_thing = svd_thing % svd_thing;
    leverages = arma::sum(svd_thing, 1);
    g_crit = leverages.max();

    return g_crit - g_crit_old;
}


