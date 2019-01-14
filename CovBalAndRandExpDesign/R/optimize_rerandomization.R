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


max_designs = 25000
n = 100
p = 10
X = matrix(rnorm(n * p), nrow = n, ncol = p)

c_vals = seq(1, 10, by = 0.25)

res = list()

for (c_val in c_vals){
  res[[as.character(c_val)]] = optimal_rerandomization(X, max_designs, c_val = c_val)
}

#save(res, file = "res_thresholds.RData")
load("res_thresholds.RData")
unlist(lapply(res, function(l){l$Q_star}))
unlist(lapply(res, function(l){l$a_star}))
unlist(lapply(res, function(l){l$W_star_size}))


#' Generates the base vectors to be used when locating the optimal rerandomization threshold
#' 
#' @param X						The data as an \code{n \times p} matrix. 
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
#' @returnType 					W_base_object
#' @return 						A list including all arguments plus a matrix \code{W_base_sorted}
#' 								whose \code{max_designs} rows are \code{n}-length allocation vectors
#' 								and the allocation vectors are in 
#' 		
#' 
#' @author AdamLenovo
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
	
  	W_base = complete_randomization_with_forced_balance(n, max_designs)  
  
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
		  imbalance_by_w = apply(W_base, 1, function(w){compute_objective_val(X, w, "mahal_dist", S_Xstd_inv)})
	  } else {
		  imbalance_by_w = apply(W_base, 1, function(w){compute_objective_val(X, w, imbalance_function)})
	  }
	  sorted_idx = sort(imbalance_by_w, index.return = TRUE)$ix
	  W_base_sorted = W_base[sorted_idx, ]
	  imbalance_by_w_sorted = imbalance_by_w[sorted_idx]
	  
	  max_designs = nrow(W_base_sorted)
	  hist(log(imbalance_by_w), br = 1000)
	  
	  ll = list(
		X = X,
		n = n,
		imbalance_function = imbalance_function,
		r = r,
		W_base_sorted = W_base_sorted, 
		max_designs = max_designs, 
		imbalance_by_w_sorted = imbalance_by_w_sorted
	  )
	  class(ll) = "W_base_object"
	  ll
}


optimal_rerandomization_exact = function(
		W_base_object = W_base_object,
		estimator = "linear",
		q = 0.95,
		z_sim_fun,
		N_z = 1000
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
	
	w_w_T_running_sum = matrix(0, n, n)
	if (estimator == "linear"){
		w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
	}
	for (s in 1 : max_designs){
		if (s %% 100 == 0){
			cat(".")
		}
		W_s = W_base_sort[s : max_designs, , drop = FALSE]
		w_s = W_base_sort[s, , drop = FALSE]
		w_s_w_s_T = t(w_s) %*% w_s
		w_w_T_running_sum = w_w_T_running_sum + 
				w_s_w_s_T
		Sigma_W = 1 / s * w_w_T_running_sum
		if (estimator == "linear"){
			
			w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum +
					w_s_w_s_T %*% P %*% w_s_w_s_T
			D = 1 / s * w_w_T_P_w_w_T_running_sum
			G = I_min_P %*% Sigma_W %*% I_min_P
			
			rel_mse_zs = array(NA, N_z)
			for (n_z in 1 : N_z){
				z = z_sim_fun()
				rel_mse_zs[n_z] = t(z) %*% (G + 2 / n * D) %*% z
			}
			Q_primes[s] = quantile(rel_mse_zs, q)
			
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
		Q_primes = Q_primes,
		balances = balances,
		frob_norm_sqs = frob_norm_sqs,
		tr_gds = tr_gds,
		tr_d_sqs = tr_d_sqs,
		r_i_sqs = r_i_sqs
	)
	
	ll = list(
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		as = imbalance_by_w_sorted,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}


optimal_rerandomization_normality_assumed = function(
		W_base_object = W_base_object,
		estimator = "linear",
		q = 0.95){
	optimal_rerandomization_argument_checks(W_base_object, estimator, q)
	
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
	
	w_w_T_running_sum = matrix(0, n, n)
	if (estimator == "linear"){
		w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
	}
	for (s in 1 : max_designs){
		if (s %% 100 == 0){
			cat(".")
		}
		W_s = W_base_sort[s : max_designs, , drop = FALSE]
		w_s = W_base_sort[s, , drop = FALSE]
		w_s_w_s_T = t(w_s) %*% w_s
		w_w_T_running_sum = w_w_T_running_sum + 
				w_s_w_s_T
		Sigma_W = 1 / s * w_w_T_running_sum
		if (estimator == "linear"){
			
			w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum +
					w_s_w_s_T %*% P %*% w_s_w_s_T
			D = 1 / s * w_w_T_P_w_w_T_running_sum
			G = I_min_P %*% Sigma_W %*% I_min_P
			
			Q_prime_zs = array(NA, N_z)
			for (n_z in 1 : N_z){
				z = z_sim_fun()
				Q_prime_zs[n_z] = t(z) %*% (G + 2 / n * D) %*% z
			}
			Q_primes[s] = mean(Q_prime_zs)
			
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
			Q_primes = Q_primes,
			balances = balances,
			frob_norm_sqs = frob_norm_sqs,
			tr_gds = tr_gds,
			tr_d_sqs = tr_d_sqs,
			r_i_sqs = r_i_sqs
	)
	
	ll = list(
			W_star = W_base_sort[1 : s_star, ],
			W_star_size = s_star,
			as = imbalance_by_w_sorted,
			a_star = imbalance_by_w_sorted[s_star],
			a_stars = imbalance_by_w_sorted[1 : s_star],
			all_data_from_run = all_data_from_run,
			Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}


optimal_rerandomization_approx = function(
		W_base_object = W_base_object,
		estimator = "linear",
		q = 0.95,
		kappa_z = 0){
	
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
	for (s in 1 : max_designs){
	  if (s %% 100 == 0){
	    cat(".")
	  }
	  W_s = W_base_sort[s : max_designs, , drop = FALSE]
	  w_s = W_base_sort[s, , drop = FALSE]
	  w_s_w_s_T = t(w_s) %*% w_s
	  w_w_T_running_sum = w_w_T_running_sum + 
	    w_s_w_s_T
	  Sigma_W = 1 / s * w_w_T_running_sum
	  if (estimator == "linear"){
	    
	    w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum +
	      w_s_w_s_T %*% P %*% w_s_w_s_T
	    D = 1 / s * w_w_T_P_w_w_T_running_sum
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
	          kappa_z * r_i_sqs[s]
	      )
	  
	  } else if (estimator == "difference_in_means"){
	    lambda_max = max(eigen(Xt %*% Sigma_W %*% X)$values)
	    Q_primes[s] = d * lambda_max + 
	      c_val * sigsq_z * sqrt(
  	      n * kappa_z + 
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
	# Q_primes[1:10000]

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
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		as = imbalance_by_w_sorted,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}

#ggplot(all_data_from_run) + 
#		geom_vline(xintercept = log(imbalance_by_w_sorted[s_star]), col = "green") +
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(Q_primes)), lwd = 2)
#ggplot(all_data_from_run) + 
#		geom_vline(xintercept = log(imbalance_by_w_sorted[s_star]), col = "green") +
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(Q_primes)), lwd = 2) +
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(imbalance_by_w_sorted)), col = "blue") +
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(frob_norm_sqs)), col = "red") +
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(tr_gds)), col = "orange") +
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(tr_d_sqs)), col = "purple") + 
#		geom_line(aes(x = log(imbalance_by_w_sorted), y = log(r_i_sqs)), col = "yellow")
