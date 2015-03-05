#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>

using namespace Rcpp;



arma::ivec primedecomp(int);
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
    // pidx is just 1:M where M is number of parents (easier to pass than write
    //  in c++!

    int M = parents.n_cols;
    arma::uvec second_parent;
    arma::cube xpxinv(n, parents.n_rows / n, M);
    arma::mat xpxi;
    arma::mat X;
    arma::vec Xv;

    arma::mat children = parents;

    int iter, child;

    int i;

    // hard code ~X1 + X2 for now.  Later use primedecomp to get actual design
    // hard code D criterion.

    // make x'x.inv

    for (i = 0; i < M; ++i)
    {
        Xv = parents.col(i);
        X = join_cols(Xv.subvec(0, n - 1), Xv.subvec(n, 2 * n - 1));
        arma::inv(xpxi, X.t() * X);
        xpxinv.slice(i) = xpxi;
    }


    for (iter=0; iter<M; ++iter)
    {
        second_parent = RcppArmadillo::sample(pidx, M, false);

        children = parents;

        for (child=0; child<M; ++child)
        {
            grblend(children, child, parents, second_parent, 0);
            grcreep(children, child, .3, 0);
            grmutat(children, child, .5);
        }




    }

    return children;
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

    arma::vec common;

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
