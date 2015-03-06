#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

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


