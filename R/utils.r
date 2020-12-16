#' Soft thresholding 
soft.thresh <- function(a, b) {
    return(sign(a) * pmax(0, abs(a) - b))
}


#' L2 Norm
l2.norm <- function(x) {
    if (!is.matrix(x)) {
	x <- as.matrix(x)
    }
    return(norm(x, type="F"))
}


#' L1 Norm
l1.norm <- function(x) {
    if (!is.matrix(x)) {
	x <- as.matrix(x)
    }
    return(norm(x, type="O"))
}


#' Normalise a vector
normalise.vect <- function(x) {
    return(x / l2.norm(x))
}
