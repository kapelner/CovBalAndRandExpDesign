
#' Generates a design matrix with standardized predictors. Useful for debugging.
#' 
#' @param n					Number of rows in the design matrix 
#' @param p 				Number of columns in the design matrix
#' @param covariate_gen		The function to use to draw the covariate realizations (assumed to be iid).
#' 							This defaults to \code{rnorm} for $N(0,1)$ draws.
#' @param ...				Optional arguments to be passed to the \code{covariate_dist} function.
#' @return 					THe design matrix
#' 
#' @author Adam Kapelner
#' @export
generate_stdzied_design_matrix = function(n = 50, p = 1, covariate_gen = rnorm, ...){
	X = matrix(covariate_gen(n * p, ...), nrow = n, ncol = p)
	#now standardize the matrix to make things easier later
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}


#' Returns the objective value given a design vector as well an an objective function.
#' This is code duplication since this is implemented within Java. This is only to be
#' run if...
#' 
#' @param X 		 	The n x p design matrix
#' @param indic_T		The n-length binary allocation vector
#' @param objective		The objective function to use. Default is \code{abs_sum_diff}.
#' @param inv_cov_X		Optional: the inverse sample variance covariance matrix. Use this
#' 						argument if you will be doing many calculations since passing this
#' 						in will cache this data.
#' 
#' @author Adam Kapelner
#' @export
compute_objective_val = function(X, indic_T, objective = "abs_sum_diff", inv_cov_X = NULL){
	X_T = X[indic_T == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_C = X[indic_T == -1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_T_bar = colMeans(X_T)
	X_C_bar = colMeans(X_C)	
	
	if (objective == "abs_sum_diff"){
		s_j_s = apply(X, 2, sd)
		sum(abs((X_T_bar - X_C_bar) / s_j_s))
	} else if (objective == "mahal_dist"){
		#saves computation to pass it in if you're doing a lot of them in a row
		if (is.null(inv_cov_X)){
			inv_cov_X = solve(var(X))
		}	
		X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar) #need to matricize for next computation
		as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)
	}
}

#compute trace of matrix
tr = function(A){sum(diag(A))}
#compute squared Frobenius Norm of matrix
frob_norm_sq = function(A){sum(A^2)}
#checks for illegal arguments
optimal_rerandomization_argument_checks = function(W_base_object, estimator, q){
	if (class(W_base_object) != "W_base_object"){
		stop("W_base_object must be class type \"W_base_object\".")
	}
	if (!(estimator %in% c("linear", "difference_in_means"))){
		stop("estimator must be either \"linear\" or \"difference_in_means\".")
	}
	if (q <= 0 || q >= 1){
		stop("quantile must be between 0 and 1.")
	}
}
