// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

double binary_search(arma::colvec v, double l, u_int niter = 1000,  
	double threshold = 1.0e-6);
arma::colvec 
rCCA_opt(arma::mat X, arma::colvec w, arma::colvec c, double mu, 
	double lambda, double l, u_int niter = 500, double threshold = 1.0e-6);
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
    const arma::mat B = X1.t() * X2;
    arma::mat U; arma::mat V; arma::vec s;
    arma::svd(U, s, V, B);
    arma::colvec w1 = U.col(1) / arma::norm(X1 * U.col(1));
    arma::colvec w2 = V.col(1) / arma::norm(X2 * V.col(1));
    arma::colvec u = arma::colvec(X1.n_cols, arma::fill::zeros);
    arma::colvec v = arma::colvec(X2.n_cols, arma::fill::zeros);
    arma::colvec w1_old = w1;
    arma::colvec w2_old = w2;

    
    for (int i = 0; i < niter; ++i) {
	w1_old = w1;
	v =  B * w2;
	d1 = norm(v, 1) <= l1 ? 0 : binary_search(v, l1);
	w1 = arma::normalise(soft_threshold(v, d1));
    
	w2_old = w2;
	u = B.t() * w1;
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

//' Computer sparse cannonical correlation vectors via Suo 2017 method
//' 
//' @param X1 Matrix of covariates.
//' @param X2 Matrix of covariates.
//' @param l1 Penalisation term.
//' @param l2 Penalisation term.
//' @param niter Number of iterations to run algorithm for (default = 1000).
//' @param threshold Stopping criterea threshold (default = 1e-6).
//' @param verbose Print debug information.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List
rCCA(arma::mat X1, arma::mat X2, double l1, double l2, u_int niter = 1000,
	double threshold = 1.0e-6, bool verbose = true)
{
    std::vector<double> loss;
    arma::mat B = X1.t() * X2;
    arma::mat U; arma::mat V; arma::vec s;
    arma::svd(U, s, V, B);
    arma::colvec w1 = U.col(1) / arma::norm(X1 * U.col(1));
    arma::colvec w2 = V.col(1) / arma::norm(X2 * V.col(1));
    arma::svd(s, X1); double mu1 = 1 / std::pow(arma::max(s), 2.0);
    arma::svd(s, X2); double mu2 = 1 / std::pow(arma::max(s), 2.0);
    
    for (u_int iteration = 0; iteration < niter; ++iteration) {
	// optimise the cannonical vectors
	w1 = rCCA_opt(X1, w1, B * w2    , mu1, 1, l1);
	w2 = rCCA_opt(X2, w2, B.t() * w1, mu2, 1, l1);

	arma::mat ls = (w1.t() * B * w2);
	loss.push_back(ls.at(0, 0));
	if (iteration > 0 && 
	    std::abs(loss.at(iteration) - loss.at(iteration-1)) < threshold) {
	    if (verbose) 
		Rcpp::Rcout << "Converged in " << iteration << "iterations";
	    break;
	}
    }
    
    // flip sign if largest weight negative
    if (arma::max(w1) < 0) w1 = -w1;
    if (arma::max(w2) < 0) w2 = -w2;

    return Rcpp::List::create(
	    Rcpp::Named("w1") = w1,
	    Rcpp::Named("w2") = w2,
	    Rcpp::Named("loss") = loss
    );
}

// [[Rcpp::export]]
arma::colvec 
rCCA_opt(arma::mat X, arma::colvec w, arma::colvec c, double mu, 
	double lambda, double l, u_int niter, double threshold)
{
    // initialise objects
    arma::colvec w_old = arma::colvec(X.n_cols, arma::fill::zeros);
    arma::colvec z = arma::colvec(X.n_rows, arma::fill::zeros);
    arma::colvec u = arma::colvec(X.n_rows, arma::fill::zeros);
    arma::colvec a = arma::colvec(X.n_rows, arma::fill::zeros);

    for (u_int i = 0; i < niter; ++i) {
	w_old = w;
	w = soft_threshold(w - (mu/lambda)* X.t() * (X * w - z + u) + mu * c,
			   mu * l);
	a = X * w;
	z = a + u;
	if (arma::norm(z, 2) > 1)
	    z = arma::normalise(z);
	u = a + u - z;
	if (arma::sum(arma::abs(w - w_old)) < threshold)
	    break;
    }
    return w;
}

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

