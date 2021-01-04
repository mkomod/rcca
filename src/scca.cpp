// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "utils.h"

//' Sparse CCA
//' @export scca
// [[Rcpp::export]]
Rcpp::List 
sCCA(arma::mat X1, arma::mat X2, double l1, double l2, 
	int niter = 1000, double threshold = 1.0e-6) 
{
    double d1 = 0.0;
    double d2 = 0.0;
    arma::colvec w1 = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec w2 = normalise(arma::colvec(X2.n_cols, arma::fill::randu));
    arma::colvec u = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec v = arma::colvec(X2.n_cols, arma::fill::zeros);
    arma::colvec w1_old = w1;
    arma::colvec w2_old = w2;

    arma::mat V = X1.t() * X2;
    arma::mat U = X2.t() * X1;
    
    for (int i = 0; i < niter; ++i) {
	w1_old = w1;
	v =  V * w2;
	d1 = norm(v, 1) <= l1 ? 0 : binary_search(v, l1);
	w1 = arma::normalise(soft_threshold(v, d1));
    
	w2_old = w2;
	u = U * w1;
	d2 = norm(u, 1) <= l2 ? 0 : binary_search(u, l2);
	w2 = arma::normalise(soft_threshold(u, d2));

	if (arma::max(arma::abs(w1_old - w1)) < threshold && 
	    arma::max(arma::abs(w2_old - w2)) < threshold)
	    break;
    }

    return Rcpp::List::create(
	    Rcpp::Named("w1") = w1,
	    Rcpp::Named("w2") = w2
    );
}


