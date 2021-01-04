#ifndef PKG_UTILS_H
#define PKG_UTILS_H

double binary_search(arma::colvec v, double l, u_int niter = 1000,  
	double threshold = 1.0e-6);
arma::colvec soft_threshold(arma::colvec v, double d);

#endif
