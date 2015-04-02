#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;

arma::vec delta_common(arma::mat xpxinv, arma::mat row_in, arma::mat row_out)
{
    arma::mat diim = row_in * xpxinv * row_in.t();
    arma::mat doom = row_out * xpxinv * row_out.t();
    arma::mat diom = row_in * xpxinv * row_out.t();

    // Fedorov values as double
    double dii = diim(0, 0);
    double doo = doom(0, 0);
    double dio = diom(0, 0);

    double delta = ((1 + dii) * (1 - doo) + dio * dio - 1);

    arma::vec common(4);
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

double get_delta_a(arma::mat xpxinv, arma::mat row_in, arma::mat row_out)
{
    arma::vec common = delta_common(xpxinv, row_in, row_out);

    double delta_d = common(0);
    double dii = common(1);
    double dio = common(2);
    double doo = common(3);

    arma::mat phiiim = row_in * xpxinv * xpxinv * row_in.t();
    arma::mat phiiom = row_out * xpxinv * xpxinv * row_in.t();
    arma::mat phioim = row_in * xpxinv * xpxinv * row_out.t();
    arma::mat phioom = row_out * xpxinv * xpxinv * row_out.t();

    // convert phi values to double
    double phiii = phiiim(0, 0);
    double phiio = phiiom(0, 0);
    double phioi = phioim(0, 0);
    double phioo = phioom(0, 0);

    double delta = ( (1 - doo) * phiii + dio * (phiio + phioi) - (1 + dii) * phioo ) / (1 + delta_d);

    return delta;
}

double get_delta_g(double g_crit_old, arma::mat X, arma::mat U_can)
{
    arma::vec s;
    arma::mat U, V;
    arma::svd_econ(U, s, V, X, "right");
    s = 1 / s;

    arma::mat Dinv = arma::diagmat(s);

    arma::mat svd_thing = U_can * V * Dinv;
    svd_thing = svd_thing % svd_thing;
    arma::vec leverages = arma::sum(svd_thing, 1);
    double g_crit = leverages.max();

    double delta = g_crit - g_crit_old;

    return delta;
}

double get_delta_i(arma::mat xpxinv, arma::mat row_in, arma::mat row_out, arma::mat B)
{
    arma::vec common = delta_common(xpxinv, row_in, row_out);

    double delta_d = common(0);
    double dii = common(1);
    double dio = common(2);
    double doo = common(3);

    arma::mat phiiim = row_in * xpxinv * B * xpxinv * row_in.t();
    arma::mat phiiom = row_out * xpxinv * B * xpxinv * row_in.t();
    arma::mat phioim = row_in * xpxinv * B * xpxinv * row_out.t();
    arma::mat phioom = row_out * xpxinv * B * xpxinv * row_out.t();

    // convert phi values to double
    double phiii = phiiim(0, 0);
    double phiio = phiiom(0, 0);
    double phioi = phioim(0, 0);
    double phioo = phioom(0, 0);

    double delta = ( (1 - doo) * phiii + dio * (phiio + phioi) - (1 + dii) * phioo ) / (1 + delta_d);

    return delta;
}
