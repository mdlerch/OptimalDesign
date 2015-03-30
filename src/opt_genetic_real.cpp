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
void grbound(arma::mat &, uint, double);
double adaptalpha(arma::vec, arma::uvec);

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat opt_geneticrealcpp(arma::mat parents, int n, arma::ivec formula,
                             int iterations, arma::uvec pidx)
{
    // parents is a matrix with each row being a vector of the design points
    //  with the first n being the first variable, the second n being the second
    //  varaible etc.
    // n is the design size
    // formula is a vector of prime numbers that allows us to build the X matrix
    //  based on just the inputs which are provided by parents.
    // iterations is the number of iterations of the genetic algorithm
    // pidx (parent index) is just 1:M where M is number of parents; easier to
    //  pass than write in c++!

    // number of parameters
    int K = formula.n_elem;

    // number of parents
    int M = parents.n_cols;

    // all parents have their own child
    arma::mat children = parents;

    // hard code ~X1 + X2 for now.  Later use primedecomp to get actual design
    arma::mat X(n, 3), xpxi(3, 3);
    arma::cube xpxinv(1 + parents.n_rows / n, 1 + parents.n_rows / n, M);

    // Change in criterion
    double delta;

    // hard code D criterion for now add other criteria later.
    for (int i = 0; i < M; ++i)
    {
        arma::vec Xv = parents.col(i);
        for (int j = 0; j < n; j++)
        {
            X(j, 0) = 1;
        }
        X.col(1) = Xv.subvec(0, n - 1);
        X.col(2) = Xv.subvec(n, 2 * n -1);
        arma::inv(xpxi, X.t() * X);
        xpxinv.slice(i) = xpxi;
    }

    arma::vec alphablend = arma::randu<arma::vec>(M);
    arma::vec alphacreep = arma::randu<arma::vec>(M);
    arma::vec alphamutat = arma::randu<arma::vec>(M);
    arma::vec alphabound = arma::randu<arma::vec>(M);
    arma::uvec swap(M);

    int iter = 0;
    while(iter < iterations)
    {
        arma::uvec second_parent = RcppArmadillo::sample(pidx, M, false);

        children = parents;

        for (int child=0; child<M; ++child)
        {
            grblend(children, child, parents, second_parent, alphablend(child));
            grcreep(children, child, .3, alphacreep(child));
            grmutat(children, child, alphamutat(child));
            grbound(children, child, alphabound(child));

            // NEED TO UPDATE THESE THINGS FOR THE STUFF
            // HARD CODE AN INTERCEPT FOR NOW.  ADD OPTION LATER.
            for (int j = 0; j < n; ++j)
            {
                X(j, 0) = 1;
            }

            X.col(1) = parents.col(child).subvec(0, n - 1);
            X.col(2) = parents.col(child).subvec(n, 2 * n - 1);

            for (int i = 0; i < n; ++i)
            {
                delta = 0;
                if (X(i, 1) != children(i, child) | X(i, 2) != children(i + n, child))
                {
                    arma::mat rowin(1, 3);
                    rowin(0, 0) = 1;
                    rowin(0, 1) = children(i, child);
                    rowin(0, 2) = children(i + n, child);
                    delta += get_delta_d(xpxinv.slice(child), rowin, X.row(i));
                    X.row(i) = rowin;
                }
            }
            swap(child) = 0;
            if (delta > 0)
            {
                if (arma::inv(xpxi, X.t() * X))
                {
                    parents.col(child) = children.col(child);
                    xpxinv.slice(child) = xpxi;
                    swap(child) = 1;
                }
            }
        }
        double newalphablend = adaptalpha(alphablend, swap);
        double newalphacreep = adaptalpha(alphacreep, swap);
        double newalphamutat = adaptalpha(alphamutat, swap);
        double newalphabound = adaptalpha(alphabound, swap);
        for (int child = 0; child < M; ++child)
        {
            if (swap(child))
            {
                alphablend(child) = newalphablend;
                alphacreep(child) = newalphacreep;
                alphamutat(child) = newalphamutat;
                alphabound(child) = newalphabound;
            }
        }

        iter++;
    }

    return parents;
}

// For each row of a child, with probability alpha switch that row to the row of
// the secondary parent
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

// For each entry of a child, with probability alpha perturb that entry with a
// normal distribution
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

// For each entry of a child, with probability alpha mutate to some new spot
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

// Jump to a boundary (-1 or 1) depending on sign of current value
void grbound(arma::mat & children, uint child, double alpha)
{
    for (int i = 0; i < children.n_rows; ++i)
    {
        if (R::runif(0, 1) < alpha)
        {
            if (children(i, child) < 0)
            {
                children(i, child) = -1;
            }
            else
            {
                children(i, child) = 1;
            }
        }
    }
}

arma::ivec primedecomp(int num)
{
    arma::ivec primes(num);

    int i = 2, idx = 0;
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
    arma::vec common = delta_common(xpxinv, row_in, row_out);
    return common(0);
}

double adaptalpha(arma::vec alpha, arma::uvec swap)
{
    double output;
    if (arma::sum(swap) > 0)
    {
        output = arma::dot(alpha, swap) / arma::sum(swap);
    }
    else
    {
        output = 0.5;
    }

    output = output + R::rnorm(0, 0.2);

    if (output > 1)
    {
        output = 1;
    }
    else if (output < 0)
    {
        output = 0;
    }

    return output;
}
