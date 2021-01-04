#' Sparse CCA
#'
#' @param X1, X2 datasets
#' @param l1, l2 
#' @return List of weights w1, w2
#'
#' @references
#' Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
#' \emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
#' @export sCCA
sCCA <- function(X1, X2, l1, l2, iteration=1000, threshold=1e-6) {
    w1 <- double(length=ncol(X1))
    # initalise w2 to be a random vector with L2-norm equal to 1
    w2 <- normalise.vect(runif(ncol(X2)))
    
    # txtProgressBar() 
    while (iteration) {
	# w1
	w1.old <- w1
	v <- t(X1) %*% X2 %*% w2
	D1 <- ifelse(l1.norm(w1) <= l1, 0, sCCA.binary_search(v , l1))
	w1 <- normalise.vect(soft.thresh(v, D1))

	# w2
	w2.old <- w2
	u <- t(X2) %*% X1 %*% w1 
	D2 <- ifelse(l1.norm(w2) <= l2, 0, sCCA.binary_search(u , l2))
	w2 <- normalise.vect(soft.thresh(u, D2))
	
	# exit the loop if converged
	if (all(abs(w1.old - w1) < threshold) && 
	    all(abs(w2.old - w2) < threshold)) {
	    break
	}

	iteration  <- iteration - 1
    }
    return(list(w1=w1, w2=w2))
}


#' Binary search for choosing the value of d
#'
#' @param v 
#' @param l
#'
#' @return d
sCCA.binary_search <- function(v, l, iteration=150) {
    
    d.lower <- 0
    d.upper <- max(abs(v))             # soft threshold of this val is
    # <= 0 forall elements in w

    while (iteration) {
	w <- normalise.vect(soft.thresh(v, (d.upper + d.lower) / 2))
	if (l1.norm(w) > l) {
	   d.lower <-  (d.lower + d.upper) / 2
	} else {
	   d.upper <-  (d.lower + d.upper) / 2
	}

	if (d.upper - d.lower < 1e-6) return((d.upper + d.lower) / 2)
	iteration <- iteration - 1
    }

    return((d.upper + d.lower) / 2)
}
