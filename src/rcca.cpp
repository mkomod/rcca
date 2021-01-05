// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double binary_search(arma::colvec v, double l, u_int niter = 1000,  
	double threshold = 1.0e-6);
arma::colvec soft_threshold(arma::colvec v, double d);


//' Computer sparse cannonical correlation vectors
//' 
//' @param X1 Matrix of covariates.
//' @param X2 Matrix of covariates.
//' @param l1 Penalisation term.
//' @param l2 Penalisation term.
//' @param niter Number of iterations to run algorithm for (default = 1000).
//' @param threshold Stopping criterea threshold (default = 1e-6).
//' 
//' @export sCCA
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


//' @export
// [[Rcpp::export(.binary_search)]]
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

//' @export
// [[Rcpp::export(.soft_threshold)]]
arma::colvec 
soft_threshold(arma::colvec v, double d) 
{
    v.transform([d](double val) { 
	return std::abs(val) > d ? arma::sign(val) * (std::abs(val) - d) : 0;
    });
    return v;
}

