#include <Rcpp.h>
#include <RcppArmadillo.h>

#include "utils.h"


double 
binary_search(arma::colvec v, double l, u_int niter = 50, double threshold = 1.0e-6)
{
    double d_lower = 0.0;
    double d_upper = max(abs(v));
    double d = (d_upper + d_lower) / 2.0;
    arma::colvec w = arma::colvec(v.n_elem, arma::fill::zeros);

    for (int i = 0; i < niter; ++i) {
	w = normalise(soft_threshold(v, d));
	d = (d_upper + d_lower) / 2.0;
	norm(w, 1) < l ? d_upper = d : d_lower = d;
	if (d_upper - d_lower <= threshold) return d;
    }

    Rcpp::Rcerr << "Failed to converge in " << niter << " iterations ";
    return d;
}


arma::colvec 
soft_threshold(arma::colvec v, double d) 
{
    arma::colvec absv = abs(v);
    return arma::sign(v) % absv.clean(d);
}


