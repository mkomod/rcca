// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

arma::colvec soft_threshold(arma::colvec v, double d);
arma::colvec cca_opt(const arma::mat X, arma::colvec w, const arma::colvec c, 
	const double mu, const double lambda, const double l, 
	const u_int niter = 500, const double threshold = 1.0e-6);

// [[Rcpp::export]]
Rcpp::List
cca(const arma::mat X1, const arma::mat X2, const double l1, const double l2, 
	const u_int niter, const double threshold, const bool verbose)
{
    arma::mat B = X1.t() * X2;
    arma::mat U; arma::mat V; arma::vec s;
    arma::svd(U, s, V, B);
    arma::colvec w1 = U.col(0) / arma::norm(X1 * U.col(0));
    arma::colvec w2 = V.col(0) / arma::norm(X2 * V.col(0));
    arma::svd(s, X1); double mu1 = 1 / std::pow(arma::max(s), 2.0);
    arma::svd(s, X2); double mu2 = 1 / std::pow(arma::max(s), 2.0);
    
    std::vector<double> loss;
    for (u_int iteration = 0; iteration < niter; ++iteration) {
	w1 = cca_opt(X1, w1, B * w2    , mu1, 1, l1);
	w2 = cca_opt(X2, w2, B.t() * w1, mu2, 1, l2);

	arma::mat ls = (w1.t() * B * w2);
	loss.push_back(ls.at(0, 0));
	if (iteration > 0 && 
	    std::abs(loss.at(iteration) - loss.at(iteration-1)) < threshold) {
	    if (verbose) 
		Rcpp::Rcout << "Converged in " << iteration << " iterations\n" <<
			       "Loss: " << loss.at(iteration) << "\n";
	    break;
	}
    }
    
    return Rcpp::List::create(
	    Rcpp::Named("w1") = w1,
	    Rcpp::Named("w2") = w2,
	    Rcpp::Named("loss") = loss,
	    Rcpp::Named("corr") = loss.back()
    );
}

// [[Rcpp::export]]
arma::colvec 
cca_opt(const arma::mat X, arma::colvec w, const arma::colvec c, const double mu, 
	const double lambda, const double l, const u_int niter, 
	const double threshold)
{
    // initialise objects
    arma::colvec w_old = arma::colvec(X.n_cols, arma::fill::zeros);
    arma::colvec z = arma::colvec(X.n_rows, arma::fill::zeros);
    arma::colvec u = arma::colvec(X.n_rows, arma::fill::zeros);
    arma::colvec a = arma::colvec(X.n_rows, arma::fill::zeros);

    for (u_int i = 0; i < niter; ++i) {
	w_old = w;
	w = soft_threshold(w - (mu/lambda)*X.t()*(X*w - z + u) + mu*c, mu*l);
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

//' @export
// [[Rcpp::export]]
double
cca_permutation_validation(const arma::mat X1, const arma::mat X2, 
	const double l1, const double l2, const u_int permutations, 
	const u_int niter, const double threshold, const bool verbose)
{
    double loss_to_beat = cca(X1, X2, l1, l2, niter, threshold, verbose)[3];
    if (loss_to_beat == 0.0)
	return -1;
    
    double total = 0;
    arma::mat S1 = X1;
    for (int perm = 0; perm < permutations; ++perm) {
	S1 = shuffle(X1, 0);
	double loss = cca(S1, X2, l1, l2, niter, threshold, verbose)[3];
	if (loss >= loss_to_beat)
	    total += 1.0;
    }
    return total / permutations;
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
