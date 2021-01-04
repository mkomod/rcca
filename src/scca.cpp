#include <Rcpp.h>
#include <RcppArmadillo.h>

#include "utils.h"

// [[Rcpp::export]]
Rcpp::List 
sCCA(arma::mat X1, arma::mat X2, double l1, double l2, 
	int niter = 1000, double threshold = 1.0e-6) 
{
    arma::colvec w1 = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec w1_old = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec w2 = arma::colvec(X2.n_cols, arma::fill::randu);
    arma::colvec w2_old = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec u = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec v = arma::colvec(X2.n_cols, arma::fill::randu);
    double d1 = 0;
    double d2 = 0;
    
    for (int i = 0; i < niter; ++i) {
	w1_old = w1;
	v = X1.t() * X2 * w2;
	d1 = norm(w1, 1) <= l1 ? 0 : binary_search(v, l1);
	w1 = normalise(soft_threshold(v, d1));

	w2_old = w2;
	u = X2.t() * X1 * w1;
	d2 = norm(w2, 1) <= l2 ? 0 : binary_search(u, l2);
	w2 = normalise(soft_threshold(u, d2));

	if (max(abs(w1_old - w1)) < threshold && max(abs(w2_old - w2)) < threshold)
	    break;
    }

    return Rcpp::List::create(
	    Rcpp::Named("u") = u,
	    Rcpp::Named("v") = v
    );
}

