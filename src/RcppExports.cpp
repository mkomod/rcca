// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sCCA_
Rcpp::List sCCA_(arma::mat X1, arma::mat X2, double l1, double l2, arma::colvec w2, u_int niter, double threshold);
RcppExport SEXP _rcca_sCCA_(SEXP X1SEXP, SEXP X2SEXP, SEXP l1SEXP, SEXP l2SEXP, SEXP w2SEXP, SEXP niterSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< double >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< u_int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(sCCA_(X1, X2, l1, l2, w2, niter, threshold));
    return rcpp_result_gen;
END_RCPP
}
// rCCA_
Rcpp::List rCCA_(arma::mat X1, arma::mat X2, double l1, double l2, u_int niter, double threshold, bool verbose);
RcppExport SEXP _rcca_rCCA_(SEXP X1SEXP, SEXP X2SEXP, SEXP l1SEXP, SEXP l2SEXP, SEXP niterSEXP, SEXP thresholdSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< double >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< u_int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rCCA_(X1, X2, l1, l2, niter, threshold, verbose));
    return rcpp_result_gen;
END_RCPP
}
// rCCA_opt
arma::colvec rCCA_opt(arma::mat X, arma::colvec w, arma::colvec c, double mu, double lambda, double l, u_int niter, double threshold);
RcppExport SEXP _rcca_rCCA_opt(SEXP XSEXP, SEXP wSEXP, SEXP cSEXP, SEXP muSEXP, SEXP lambdaSEXP, SEXP lSEXP, SEXP niterSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< u_int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(rCCA_opt(X, w, c, mu, lambda, l, niter, threshold));
    return rcpp_result_gen;
END_RCPP
}
// binary_search
double binary_search(arma::colvec v, double l, u_int niter, double threshold);
RcppExport SEXP _rcca_binary_search(SEXP vSEXP, SEXP lSEXP, SEXP niterSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< u_int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(binary_search(v, l, niter, threshold));
    return rcpp_result_gen;
END_RCPP
}
// soft_threshold
arma::colvec soft_threshold(arma::colvec v, double d);
RcppExport SEXP _rcca_soft_threshold(SEXP vSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_threshold(v, d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rcca_sCCA_", (DL_FUNC) &_rcca_sCCA_, 7},
    {"_rcca_rCCA_", (DL_FUNC) &_rcca_rCCA_, 7},
    {"_rcca_rCCA_opt", (DL_FUNC) &_rcca_rCCA_opt, 8},
    {"_rcca_binary_search", (DL_FUNC) &_rcca_binary_search, 4},
    {"_rcca_soft_threshold", (DL_FUNC) &_rcca_soft_threshold, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rcca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
