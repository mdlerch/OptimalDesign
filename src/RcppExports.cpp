// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fedorovcpp
arma::ivec fedorovcpp(const arma::mat& xpxinv, const arma::mat& X, arma::ivec current, arma::ivec complete, int iter);
RcppExport SEXP OptimalDesign_fedorovcpp(SEXP xpxinvSEXP, SEXP XSEXP, SEXP currentSEXP, SEXP completeSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type xpxinv(xpxinvSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::ivec >::type current(currentSEXP );
        Rcpp::traits::input_parameter< arma::ivec >::type complete(completeSEXP );
        Rcpp::traits::input_parameter< int >::type iter(iterSEXP );
        arma::ivec __result = fedorovcpp(xpxinv, X, current, complete, iter);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// initDesign
arma::uvec initDesign(const arma::mat& X, int design_size);
RcppExport SEXP OptimalDesign_initDesign(SEXP XSEXP, SEXP design_sizeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP );
        Rcpp::traits::input_parameter< int >::type design_size(design_sizeSEXP );
        arma::uvec __result = initDesign(X, design_size);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
