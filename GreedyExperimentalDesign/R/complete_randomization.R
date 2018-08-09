#' Implements complete balanced randomization
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization_with_balanced_ns = function(n, r){
	indicTs = matrix(NA, nrow = r, ncol = n)
	
	for (nsim in 1 : r){
		indicTs[nsim, ] = sample(c(rep(1, n / 2), rep(-1, n / 2)))
	}
	indicTs
}


#' Implements complete randomization
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization = function(n, r){
	indicTs = matrix(NA, nrow = r, ncol = n)
	
	for (nsim in 1 : r){
		indicTs[nsim, ] = 2 * (rbinom(n, 1, 0.5) - 0.5)
	}
	indicTs
}
