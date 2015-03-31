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
arma::vec getterm(int, arma::vec, int);
arma::ivec orderprimes(arma::ivec);

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

    // K + 1 for K in formula plus intercept
    arma::mat X(n, K + 1), xpxi(K + 1, K + 1);
    arma::cube xpxinv(K + 1, K + 1, M);

    // Change in criterion
    double delta;

    // hard code D criterion for now add other criteria later.
    for (int child = 0; child < M; ++child)
    {
        // hard code an intercept
        for (int j = 0; j < n; j++)
        {
            X(j, 0) = 1;
        }
        for (int term = 0; term < K; ++term)
        {
            X.col(term + 1) = getterm(formula(term), parents.col(child), n);
        }
        arma::inv(xpxi, X.t() * X);
        xpxinv.slice(child) = xpxi;
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

            // HARD CODE AN INTERCEPT FOR NOW.  ADD OPTION LATER.
            for (int j = 0; j < n; ++j)
            {
                X(j, 0) = 1;
            }

            // step over columns (terms in model) of X matrix and compute that
            // column
            for (int term = 0; term < K; ++term)
            {
                X.col(term + 1) = getterm(formula(term), parents.col(child), n);
            }

            // accumulate delta for all rows that changed
            delta = 0;
            for (int x_row = 0; x_row < n; ++x_row)
            {
                // flag to indicate this row has changed
                int flag = 0;
                for (int k = 0; k < K; ++k)
                {
                    // If any of the childs entry that goes into the nth row has
                    // changed
                    int idx = k * n + x_row;
                    if (parents(idx, child) != children(idx, child))
                    {
                        flag = 1;
                        k = K + 1;
                    }
                }
                if (flag)
                {
                    arma::mat rowin(1, K + 1);
                    rowin(0, 0) = 1;
                    for (int k = 0; k < K; ++k)
                    {
                        rowin(0, k + 1) = (getterm(formula(k), children.col(child), n))(x_row);
                    }
                    delta += get_delta_d(xpxinv.slice(child), rowin, X.row(x_row));
                    X.row(x_row) = rowin;
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
    for (int i=0; i<children.n_rows; ++i)
    {
        if (R::runif(0, 1) < alpha)
        {
            children(i, child) = parents(i, parent2[child]);
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

arma::ivec orderprimes(arma::ivec primes)
{
    arma::ivec indices(primes.n_elem);

    for (int i = 0; i < primes.n_elem; ++i)
    {
        if (primes(i) == 2)
        {
            indices(i) = 1;
        }
        else if (primes(i) == 3)
        {
            indices(i) = 2;
        } else if (primes(i) == 5)
        {
            indices(i) = 3;
        } else if (primes(i) == 7)
        {
            indices(i) = 4;
        }
        else if (primes(i) == 11)
        {
            indices(i) = 5;
        }
        else if (primes(i) == 13)
        {
            indices(i) = 6;
        }
        else if (primes(i) == 17)
        {
            indices(i) = 7;
        }
        else if (primes(i) == 19)
        {
            indices(i) = 8;
        }
        else if (primes(i) == 23)
        {
            indices(i) = 9;
        }
        else if (primes(i) == 29)
        {
            indices(i) = 10;
        }
        else
        {
            indices(i) = 1;
        }
    }
    return indices;
}

arma::vec getterm(int formula, arma::vec x, int n)
{
    arma::ivec primes = primedecomp(formula);
    arma::ivec factors = orderprimes(primes);

    arma::vec column(n);
    column.ones();

    for (int i = 0; i < factors.n_elem; ++i)
    {
        int start = (factors(i) - 1) * n;
        int stop = factors(i) * n - 1;
        column = column % x.subvec(start, stop);
    }
    return column;
}
