#include <RcppArmadillo.h>
#include "utils.h"

// [[Rcpp::export]]
double 
binary_search(arma::colvec v, double l, u_int niter, double threshold)
{
    double d_lower = 0.0;
    double d_upper = max(abs(v));
    double d = 0.0;
    arma::colvec w = arma::colvec(v.n_elem, arma::fill::zeros);

    for (int i = 0; i < niter; ++i) {
	d = (d_upper + d_lower) / 2.0;
	w = arma::normalise(soft_threshold(v, d));
	norm(w, 1) > l ? d_lower = d : d_upper = d;
	if (d_upper - d_lower <= threshold) 
	    return d;
    }

    Rcpp::Rcerr << "Failed to converge in " << niter << " iterations ";
    return d;
}

// [[Rcpp::export]]
arma::colvec 
soft_threshold(arma::colvec v, double d) 
{
    v.transform([d](double val) { 
	return std::abs(val) > d ? arma::sign(val) * (std::abs(val) - d) : 0;
    });
    return v;
}

