
#' Sparse Canonical Correlation Analysis
#' 
#' @param X1,X2 Matrices of covariates.
#' @param l1,l2 Penalisation terms.
#' @param K Number of canonical vector pairs.
#' @param niter Number of iterations to run algorithm for (default = 1000).
#' @param threshold Stopping criterea threshold (default = 1e-6).
#' @param verbose Print debug information.
#' @return List with \code{K} canonical vectors \code{w1} and \code{w2}.
#'
#' @examples
#' # From PMA::CCA
#' u <- matrix(c(rep(1,25),rep(0,75)),ncol=1)
#' v1 <- matrix(c(rep(1,50),rep(0,450)),ncol=1)
#' v2 <- matrix(c(rep(0,50),rep(1,50),rep(0,900)),ncol=1)
#' X1 <- u%*%t(v1) + matrix(rnorm(100*500),ncol=500)
#' X2 <- u%*%t(v2) + matrix(rnorm(100*1000),ncol=1000)
#' X1 <- scale(X1, T, F)
#' X2 <- scale(X2, T, F)
#' cca <- sCCA(X1, X2)
#' 
#' @export
sCCA <- function(X1, X2, l1 = 0.3 * sqrt(ncol(X1)), l2 = 0.3 * sqrt(ncol(X2)),
                 K = 1, niter = 1000, threshold = 1.0e-6, verbose = TRUE)
{

    SVD <- svd(t(X1) %*% X2, nu=0, nv=K)
    W1 <- matrix(nrow=ncol(X1), ncol=0)
    W2 <- matrix(nrow=ncol(X2), ncol=0)
    for (i in 1:K) {
        w2  <- matrix(SVD$v[ , i], ncol=1)
        res <- sCCA_(X1, X2, l1, l2, w2, niter, threshold) # Cpp function
        X1 <- rbind(X1,  sqrt(res$d[1, 1]) * t(res$w1))
        X2 <- rbind(X2, -sqrt(res$d[1, 1]) * t(res$w2))
        W1 <- cbind(W1, res$w1)
        W2 <- cbind(W2, res$w2)
	if (verbose) cat("Loss: ", res$d, "\n")
    }
    return(list(w1=W1, w2=W2))
}


#' @rdname sCCA
#' @export
rCCA <- function(X1, X2, l1 = 0.3, l2 = 0.3, K = 1, niter = 1000, 
		 threshold = 1.0e-6, verbose = TRUE) 
{
    W1 <- matrix(nrow=ncol(X1), ncol=0)
    W2 <- matrix(nrow=ncol(X2), ncol=0)
    loss <- double(length=K)
    for (i in 1:K) {
	res <- rCCA_(rbind(X1, t(W1) %*% t(X1) %*% X1, t(W2) %*% t(X2) %*% X1), 
		     rbind(X2, t(W2) %*% t(X2) %*% X2, t(W1) %*% t(X1) %*% X2), 
		    l1, l2, niter, threshold, verbose)
	W1 <- cbind(W1, res$w1)
	W2 <- cbind(W2, res$w2)
	loss[i] <- res$loss[length(res$loss)]
    }
    return(list(w1=W1, w2=W2, loss=loss))
}

