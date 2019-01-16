optimal_rerandomization_sets = function(X, 
  max_designs = 25000,
  objective = "mahal_dist",
  estimator = "linear",
  c_grid = seq(from = 1, to = 5, by = 0.5),
  kappa_z_grid = seq(from = -1, to = 1, by = 0.5),
  ...){
  
  if (estimator == "linear"){
    all_results_a_star = array(NA, c(length(c_grid), length(kappa_z_grid)))
    dimnames(all_results_a_star) = list(c_grid, kappa_z_grid)
    all_results_Q_star = array(NA, c(length(c_grid), length(kappa_z_grid)))
    dimnames(all_results_Q_star) = list(c_grid, kappa_z_grid)
    all_results_W_star_size = array(NA, c(length(c_grid), length(kappa_z_grid)))
    dimnames(all_results_W_star_size) = list(c_grid, kappa_z_grid)
    for (i_c in 1 : length(c_grid)){
      c_val = c_grid[i_c]
      for (i_k in 1 : length(kappa_z_grid)){
        kappa_z = kappa_z_grid[i_k]
        
        results = optimal_rerandomization(
          X = X,
          max_designs = max_designs,
          objective = objective,
          estimator = estimator,
          c_val = c_val,
          kappa_z = kappa_z
        )
        
        all_results_a_star[i_c, i_k] = results$a_star
        all_results_Q_star[i_c, i_k] = results$Q_star
        all_results_W_star_size[i_c, i_k] = results$W_star_size
      }
    }
  }
  list(all_results_a_star = all_results_a_star, all_results_Q_star = all_results_Q_star, all_results_W_star_size = all_results_W_star_size)
  load( "different_c_kappa.RData")
  #save(ll, file = "different_c_kappa.RData")
}





#' Generate Base Assignments and Sorts
#' 
#' Generates the base vectors to be used when locating the optimal rerandomization threshold
#' 
#' @param X						The data as an \eqn{n \times p} matrix. 
#' @param max_designs 			The maximum number of designs
#' @param imbalance_function 	A string indicating the imbalance function. Currently,
#' 								"abs_sum_difference" and "mahal_dist" are the options with the
#' 								latter being the default.
#' @param r 					An experimental feature that adds lower imbalance vectors
#' 								to the base set using the \code{GreedyExperimentalDesign}
#' 								package. This controls the number of vectors to search through
#' 								on each iteration.
#' @param max_max_iters			An experimental feature that adds lower imbalance vectors
#' 								to the base set using the \code{GreedyExperimentalDesign}
#' 								package. The maximum number of iterations to use for the greedy search.
#' @return 						A list including all arguments plus a matrix \code{W_base_sorted}
#' 								whose \code{max_designs} rows are \code{n}-length allocation vectors
#' 								and the allocation vectors are in 
#' 		
#' 
#' @author Adam Kapelner
#' @export
generate_W_base_and_sort = function(X,
  max_designs = 25000,
  imbalance_function = "mahal_dist",
  r = 0,
  max_max_iters = 5
){
	n = nrow(X)
	#rewrite it as standardized
	X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
	
  	W_base = complete_randomization_with_forced_balance_plus_one_min_one(n, max_designs)  
  
  	if (r > 0){
	    library(GreedyExperimentalDesign)
	    
	    rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4)
	    res = resultsGreedySearch(rd, max_vectors = r)
	    all_vecs = t(res$ending_indicTs)
	    all_vecs[all_vecs == 0] = -1
	    W_base = rbind(all_vecs, W_base)
	    
		for (max_iters in 1 : max_max_iters){
			rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4, max_iters = max_iters)
			res = resultsGreedySearch(rd, max_vectors = r)
			all_vecs = t(res$ending_indicTs)
			all_vecs[all_vecs == 0] = -1
			W_base = rbind(all_vecs, W_base)
		}
	    
#	    all_mahal_obj = apply(W_base, 1, function(w){compute_objective_val(X, w, objective = "mahal_dist")})
#	    hist(all_mahal_obj_sort_log, br = 1000)
	    
	    
	    ##now create a representative sample
	    bin_width_in_log = 0.1
	    tab = table(round(all_mahal_obj_sort_log / bin_width_in_log))
	    
	    num_per_bin_star = floor(uniroot(function(num_per_bin){
	      sum(pmin(tab, num_per_bin)) - max_designs
	    }, interval = c(1, max_designs))$root)
	    sum(pmin(tab, num_per_bin_star))
	    
	    all_vecs_by_bin = list()
	    for (i in 1 : nrow(W_base)){
	      idx = as.character(round(all_mahal_obj_sort_log[i] / bin_width_in_log))
	      if (is.null(all_vecs_by_bin[[idx]])){
	        all_vecs_by_bin[[idx]] = matrix(NA, nrow = 0, ncol = n)
	      }
	      if (nrow(all_vecs_by_bin[[idx]]) < num_per_bin_star){
	        all_vecs_by_bin[[idx]] = rbind(all_vecs_by_bin[[idx]], W_base[i, ])
	      }
	    }
	    unlist(lapply(all_vecs_by_bin, nrow))
	    W_base = df <- do.call("rbind", all_vecs_by_bin)
	  }
  
	  if (imbalance_function == "mahal_dist"){
		  S_Xstd_inv = solve(var(X))
		  imbalance_by_w = apply(W_base, 1, function(w){compute_objective_val_plus_one_min_one_enc(X, w, "mahal_dist", S_Xstd_inv)})
	  } else {
		  imbalance_by_w = apply(W_base, 1, function(w){compute_objective_val_plus_one_min_one_enc(X, w, imbalance_function)})
	  }
	  sorted_idx = sort(imbalance_by_w, index.return = TRUE)$ix
	  	  
	  ll = list(
		X = X,
		n = n,
		imbalance_function = imbalance_function,
		W_base_sorted = W_base[sorted_idx, ], 
		max_designs = nrow(W_base), 
		imbalance_by_w_sorted = imbalance_by_w[sorted_idx]
	  )
	  class(ll) = "W_base_object"
	  ll
}

#' Plots a summary of the imbalances in a \code{W_base_object} object
#' 
#' @param x			The \code{W_base_object} object to be summarized in the plot
#' @param ...		\code{title}, \code{subtitle}, \code{xlab}, \code{bins} can 
#' 					be specified here to be passed to the ggplot plotting function.
#' 					Also \code{log10} can be set to \code{FALSE} to not log the x-axis.
#' 
#' @author 			Adam Kapelner
#' @method plot W_base_object
#' @export
plot.W_base_object = function(x, ...){
	dots = list(...)
	if (is.null(dots$title)){
		title = "Density of Imbalances in Base Strategy"
	} else {
		title = dots$title
	}
	if (is.null(dots$subtitle)){
		subtitle = ""
	} else {
		subtitle = dots$subtitle
	}
	if (is.null(dots$xlab)){
		xlab = x$imbalance_function
	} else {
		xlab = dots$xlab
	}
	if (is.null(dots$bins)){
		bins = x$max_designs / 10
	} else {
		bins = dots$bins
	}

	ggplot_obj = ggplot(data.frame(b = x$imbalance_by_w_sorted)) + 
			aes(x = b) + 
			ggtitle(title, subtitle = subtitle) +
			xlab(xlab) +
			geom_histogram(bins = bins)
	if (!isFALSE(dots$log10)){
		ggplot_obj = ggplot_obj + scale_x_log10()
	}
	
	plot(ggplot_obj)
}

#' Prints a summary of a \code{W_base_object} object
#' 
#' @param x			The \code{W_base_object} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print W_base_object
#' @export
print.W_base_object = function(x, ...){	
	cat("W base strategy with", x$max_designs, "assignments whose imbalances range from",
			round(min(x$imbalance_by_w_sorted), 3), "to", round(max(x$imbalance_by_w_sorted), 3), 
			"in", x$imbalance_function, "\n")
}

#' Prints a summary of a \code{W_base_object} object
#' 
#' @param object		The \code{W_base_object} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary W_base_object
#' @export
summary.W_base_object = function(object, ...){
	print(object, ...)
}


#' Find the Optimal Rerandomization Design Exactly
#' 
#' Finds the optimal rerandomization threshold based on a user-defined quantile
#' and a function that generates the non-linear component of the response
#' 
#' @param W_base_object			An object that contains the assignments to begin with sorted by 
#' @param estimator 			"linear" for the covariate-adjusted linear regression estimator (default)
#' 								or "difference_in_means".
#' @param q 					The tail criterion's quantile of MSE over z's. The default is 95\%. 
#' @param skip_search_length	In the exhaustive search, how many designs are skipped? Default is 1 for 
#' 								full exhaustive search through all assignments provided for in \code{W_base_object}.
#' @param binary_search			If \code{TRUE}, a binary search is employed to find the optimal threshold instead of 
#' 								an exhaustive search. Default is \code{FALSE}.
#' @param z_sim_fun 			This function returns vectors of numeric values of size \code{n}. No default is provided.
#' @param N_z 					The number of times to simulate z's within each strategy.
#' @param dot_every_x_iters		Print out a dot every this many iterations. The default is 100. Set to
#' 								\code{NULL} for no printout.
#' @return 						A list containing the optimal design threshold, strategy, and
#' 								other information.
#' 
#' @author Adam Kapelner
#' @export
optimal_rerandomization_exact = function(
		W_base_object = W_base_object,
		estimator = "linear",
		q = 0.95,
		skip_search_length = 1,
		binary_search = FALSE,
		z_sim_fun,
		N_z = 1000,
		dot_every_x_iters = 100
){
	optimal_rerandomization_argument_checks(W_base_object, estimator, q)

	n = W_base_object$n
	X = W_base_object$X
	W_base_sort = W_base_object$W_base_sort
	max_designs = W_base_object$max_designs
	imbalance_by_w_sorted = W_base_object$imbalance_by_w_sorted
	
	if (estimator == "linear"){
		Xt = t(X)
		XtXinv = solve(Xt %*% X)
		P = X %*% XtXinv %*% Xt
		I = diag(n)
		I_min_P = I - P
	}
	
	s_star = NULL
	Q_star = Inf
	Q_primes = array(NA, max_designs)
	rel_mse_zs = matrix(NA, nrow = max_designs, ncol = N_z)
	
	w_w_T_running_sum = matrix(0, n, n)
	if (estimator == "linear"){
		w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
	}
	for (s in seq(from = 1, to = max_designs, by = skip_search_length)){
		if (!is.null(dot_every_x_iters)){
			if (s %% dot_every_x_iters == 0){
				cat(".")
			}
		}
		w_s = W_base_sort[s, , drop = FALSE]
		w_s_w_s_T = t(w_s) %*% w_s
		w_w_T_running_sum = w_w_T_running_sum + w_s_w_s_T
		Sigma_W = 1 / (s / skip_search_length) * w_w_T_running_sum
		if (estimator == "linear"){			
			w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum + w_s_w_s_T %*% P %*% w_s_w_s_T
			D = 1 / (s / skip_search_length) * w_w_T_P_w_w_T_running_sum
			G = I_min_P %*% Sigma_W %*% I_min_P
			
			
			for (n_z in 1 : N_z){
				z = z_sim_fun()
				rel_mse_zs[s, n_z] = t(z) %*% (G + 2 / n * D) %*% z
			}
			Q_primes[s] = quantile(rel_mse_zs[s, ], q)
			
		} else if (estimator == "difference_in_means"){

		}
		
		
		if (Q_primes[s] < Q_star){
			Q_star = Q_primes[s]
			s_star = s
		}
	}
	cat("\n")
	# Q_primes[1:10000]
	
	all_data_from_run = data.frame(
		imbalance_by_w_sorted = imbalance_by_w_sorted, 
		rel_mse_zs = rel_mse_zs,
		Q_primes = Q_primes
	)
	
	ll = list(
		type = "exact",
		q = q,
		estimator = estimator,
		z_sim_fun = z_sim_fun,
		N_z = N_z,
		imbalance_function = W_base_object$imbalance_function,
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}

#' Find the Optimal Rerandomization Design Under the Gaussian Approximation
#' 
#' Finds the optimal rerandomization threshold based on a user-defined quantile
#' and a function that generates the non-linear component of the response
#' 
#' @param W_base_object			An object that contains the assignments to begin with sorted by 
#' @param estimator 			"linear" for the covariate-adjusted linear regression estimator (default)
#' 								or "difference_in_means".
#' @param q 					The tail criterion's quantile of MSE over z's. The default is 95\%. 
#' @param skip_search_length	In the exhaustive search, how many designs are skipped? Default is 1 for 
#' 								full exhaustive search through all assignments provided for in \code{W_base_object}.
#' @param binary_search			If \code{TRUE}, a binary search is employed to find the optimal threshold instead of 
#' 								an exhaustive search. Default is \code{FALSE}.
#' @param dot_every_x_iters		Print out a dot every this many iterations. The default is 100. Set to
#' 								\code{NULL} for no printout.
#' @return 						A list containing the optimal design threshold, strategy, and
#' 								other information.
#' 
#' @author Adam Kapelner
#' @export
optimal_rerandomization_normality_assumed = function(
		W_base_object = W_base_object,
		estimator = "linear",
		q = 0.95,
		skip_search_length = 1,
		binary_search = FALSE,
		dot_every_x_iters = 100){
	optimal_rerandomization_argument_checks(W_base_object, estimator, q)
	
	n = W_base_object$n
	X = W_base_object$X
	W_base_sort = W_base_object$W_base_sort
	max_designs = W_base_object$max_designs
	imbalance_by_w_sorted = W_base_object$imbalance_by_w_sorted
	
	if (estimator == "linear"){
		Xt = t(X)
		XtXinv = solve(Xt %*% X)
		P = X %*% XtXinv %*% Xt
		I = diag(n)
		I_min_P = I - P
	}
	
	s_star = NULL
	Q_star = Inf
	Q_primes = array(NA, max_designs)
	
	w_w_T_running_sum = matrix(0, n, n)
	if (estimator == "linear"){
		w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
	}
	for (s in seq(from = 1, to = max_designs, by = skip_search_length)){
		if (!is.null(dot_every_x_iters)){
			if (s %% dot_every_x_iters == 0){
				cat(".")
			}
		}
		w_s = W_base_sort[s, , drop = FALSE]
		w_s_w_s_T = t(w_s) %*% w_s
		w_w_T_running_sum = w_w_T_running_sum + w_s_w_s_T
		Sigma_W = 1 / (s / skip_search_length) * w_w_T_running_sum
		if (estimator == "linear"){			
			w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum + w_s_w_s_T %*% P %*% w_s_w_s_T
			D = 1 / (s / skip_search_length) * w_w_T_P_w_w_T_running_sum
			G = I_min_P %*% Sigma_W %*% I_min_P
			eigenvalues = eigen(G + 2 / n * D)$values
			Q_primes[s] = hall_buckley_eagleson_inverse_cdf(eigenvalues, q, n)			
		} else if (estimator == "difference_in_means"){
			
		}
		
		
		if (Q_primes[s] < Q_star){
			Q_star = Q_primes[s]
			s_star = s
		}
	}
	cat("\n")
	
	all_data_from_run = data.frame(
		imbalance_by_w_sorted = imbalance_by_w_sorted, 
		Q_primes = Q_primes
	)
	
	ll = list(
		type = "normal",
		estimator = estimator,
		q = q,
		imbalance_function = W_base_object$imbalance_function,
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}

#' Find the Optimal Rerandomization Design Under the Tail and Kurtosis Approximation
#' 
#' Finds the optimal rerandomization threshold based on a user-defined quantile
#' and kurtosis based on an approximation of tail standard errors
#' 
#' @param W_base_object			An object that contains the assignments to begin with sorted by 
#' @param estimator 			"linear" for the covariate-adjusted linear regression estimator (default)
#' 								or "difference_in_means".
#' @param q 					The tail criterion's quantile of MSE over z's. The default is 95\%. 
#' @param skip_search_length	In the exhaustive search, how many designs are skipped? Default is 1 for 
#' 								full exhaustive search through all assignments provided for in \code{W_base_object}.
#' @param binary_search			If \code{TRUE}, a binary search is employed to find the optimal threshold instead of 
#' 								an exhaustive search. Default is \code{FALSE}.
#' @param excess_kurtosis_z		An estimate of the excess kurtosis in the measure on z. Default is 0.
#' @param dot_every_x_iters		Print out a dot every this many iterations. The default is 100. Set to
#' 								\code{NULL} for no printout.
#' @return 						A list containing the optimal design threshold, strategy, and
#' 								other information.
#' 
#' @author Adam Kapelner
#' @export
optimal_rerandomization_tail_approx = function(
		W_base_object = W_base_object,
		estimator = "linear",
		q = 0.95,
		skip_search_length = 1,
		binary_search = FALSE,
		excess_kurtosis_z = 0,
		dot_every_x_iters = 100){
	
	optimal_rerandomization_argument_checks(W_base_object, estimator, q)
	
	c_val = qnorm(q)
	n = W_base_object$n
	X = W_base_object$X
	W_base_sort = W_base_object$W_base_sort
	max_designs = W_base_object$max_designs
	imbalance_by_w_sorted = W_base_object$imbalance_by_w_sorted
	
	if (estimator == "linear"){
	  Xt = t(X)
	  XtXinv = solve(Xt %*% X)
	  XtXinv_eigen <- eigen(XtXinv)
	  XtXinv_sqrt <- XtXinv_eigen$vectors %*% diag(sqrt(XtXinv_eigen$values)) %*% solve(XtXinv_eigen$vectors)
	  X_orth = X %*% XtXinv_sqrt
	  X_orth_T = t(X_orth)
	  P = X %*% XtXinv %*% Xt
	  I = diag(n)
	  I_min_P = I - P
	}

	s_star = NULL
	Q_star = Inf
	Q_primes = array(NA, max_designs)
	#diagram all terms
	balances = array(NA, max_designs)
	frob_norm_sqs = array(NA, max_designs)
	tr_gds = array(NA, max_designs)
	tr_d_sqs = array(NA, max_designs)
	r_i_sqs = array(NA, max_designs)
	
	w_w_T_running_sum = matrix(0, n, n)
	if (estimator == "linear"){
	  w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
	}
	for (s in seq(from = 1, to = max_designs, by = skip_search_length)){
		if (!is.null(dot_every_x_iters)){
			if (s %% dot_every_x_iters == 0){
				cat(".")
			}
		}

	  w_s = W_base_sort[s, , drop = FALSE]
	  w_s_w_s_T = t(w_s) %*% w_s
	  w_w_T_running_sum = w_w_T_running_sum + w_s_w_s_T
	  Sigma_W = 1 / (s / skip_search_length) * w_w_T_running_sum
	  if (estimator == "linear"){	    
	    w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum + w_s_w_s_T %*% P %*% w_s_w_s_T
	    D = 1 / (s / skip_search_length) * w_w_T_P_w_w_T_running_sum
	    G = I_min_P %*% Sigma_W %*% I_min_P
	    balances[s] = tr(X_orth_T %*% Sigma_W %*% X_orth)
	    frob_norm_sqs[s] = frob_norm_sq(I_min_P %*% Sigma_W)
	    tr_gds[s] = tr(G %*% D) / n
	    tr_d_sqs[s] = tr(D^2) / n^2
	    r_i_sqs[s] = sum((diag(G) + 2 * diag(D) / n)^2)
	    Q_primes[s] = balances[s] +
	      c_val * sqrt(
	        2 * frob_norm_sqs[s] + 
	          8 * tr_gds[s] +
	          8 * tr_d_sqs[s] + 
			  excess_kurtosis_z * r_i_sqs[s]
	      )
	  
	  } else if (estimator == "difference_in_means"){
	    lambda_max = max(eigen(Xt %*% Sigma_W %*% X)$values)
	    Q_primes[s] = d * lambda_max + 
	      c_val * sigsq_z * sqrt(
  	      n * excess_kurtosis_z + 
  	        2 * frob_norm_sq(Sigma_W) + 
  	        4 * d / sigsq_z * lambda_max
  	    )
	  }
	  
	  if (Q_primes[s] < Q_star){
	    Q_star = Q_primes[s]
	    s_star = s
	  }
	}
	cat("\n")

	all_data_from_run = data.frame(
	  imbalance_by_w_sorted = imbalance_by_w_sorted, 
	  Q_primes = Q_primes,
	  balances = balances,
	  frob_norm_sqs = frob_norm_sqs,
	  tr_gds = tr_gds,
	  tr_d_sqs = tr_d_sqs,
	  r_i_sqs = r_i_sqs
	)

	ll = list(
		type = "approx",
		estimator = estimator,
		q = q,
		imbalance_function = W_base_object$imbalance_function,
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}


#' Prints a summary of a \code{optimal_rerandomization_obj} object
#' 
#' @param x			The \code{optimal_rerandomization_obj} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print optimal_rerandomization_obj
#' @export
print.optimal_rerandomization_obj = function(x, ...){	
	cat("Optimal rerandomization found with", x$W_star_size, "assignments whose imbalances are smaller\nthan",
			round(x$a_star, 3), "in", x$imbalance_function, "using algorithm type", x$type, "at q =", x$q, "\n")
}

#' Prints a summary of a \code{optimal_rerandomization_obj} object
#' 
#' @param object		The \code{optimal_rerandomization_obj} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary optimal_rerandomization_obj
#' @export
summary.optimal_rerandomization_obj = function(object, ...){
	print(object, ...)
}

#' Plots a summary of a \code{optimal_rerandomization_obj} object
#' 
#' @param x			The \code{optimal_rerandomization_obj} object to be summarized in the plot
#' @param ...		The option \code{advanced = TRUE} can be passed here for optimal rerandomization 
#' 					results from algorithm type "approx" to see how all the terms in the criterion behave.
#' 					Also, \code{title}, \code{subtitle}, \code{xlab} and \code{ylab} can be passed here.
#' 
#' @author 			Adam Kapelner
#' @method plot optimal_rerandomization_obj
#' @export
plot.optimal_rerandomization_obj = function(x, ...){
	dots = list(...)
	if (is.null(dots$title)){
		title = "Optimal rerandomization by Tail Criterion and All Terms"
	} else {
		title = dots$title
	}
	if (is.null(dots$subtitle)){
		subtitle = "optimal indicated by green line"
	} else {
		subtitle = dots$subtitle
	}
	if (is.null(dots$xlab)){
		xlab = x$imbalance_function
	} else {
		xlab = dots$xlab
	}
	if (is.null(dots$ylab)){
		ylab = paste("Relative MSE Tail at q =", x$q)
	} else {
		ylab = dots$ylab
	}
	
	
	if (x$type == "approx" && isTRUE(dots$advanced)){
		plot(ggplot(x$all_data_from_run) +
			ggtitle(title, subtitle = subtitle) +
			xlab(xlab) +
			ylab(ylab) +
			scale_x_log10() +
			scale_y_log10() +
#			coord_trans(x = "log10", y = "log10") +
			geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
			geom_line(aes(x = imbalance_by_w_sorted, y = imbalance_by_w_sorted), col = "blue") +
			geom_line(aes(x = imbalance_by_w_sorted, y = frob_norm_sqs), col = "red") +
			geom_line(aes(x = imbalance_by_w_sorted, y = tr_gds), col = "orange") +
			geom_line(aes(x = imbalance_by_w_sorted, y = tr_d_sqs), col = "purple") + 
			geom_line(aes(x = imbalance_by_w_sorted, y = r_i_sqs), col = "yellow") +			
			geom_vline(xintercept = log(x$a_star), col = "green"))
	} else {
		plot(ggplot(x$all_data_from_run) +
			ggtitle(title, subtitle = subtitle) +
			xlab(xlab) +
			ylab(ylab) +
			scale_x_log10() +
			scale_y_log10() +
#			coord_trans(x = "log10", y = "log10") +
			geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
			geom_vline(xintercept = x$a_star, col = "green"))
	}

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
compute_objective_val_plus_one_min_one_enc = function(X, indic_T, objective = "abs_sum_diff", inv_cov_X = NULL){
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
#compute inverse CDF as a function of desired quantile using the hall-buckley-eagleson method
hall_buckley_eagleson_inverse_cdf = function(eigenvalues, q, sample_size, tol = 0.001){
	eigenvalues = eigenvalues[eigenvalues > 0] #only use non-zero eigenvalues
	if (any(eigenvalues > sample_size)){
		eigenvalues = c(sample_size) #theoretical upper limit
	}
	fun = function(x, eigenvalues){
		#cat("x =", x, "eigenvalues =", eigenvalues)
		hbe(coeff = eigenvalues, x = x) - q
	}
	uniroot(fun, eigenvalues, interval = c(tol, sample_size * qchisq(q, 1)) * 1.1, tol = tol)$root #1.1 is a fudge for numerical error.
}
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