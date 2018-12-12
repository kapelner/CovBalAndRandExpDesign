options(java.parameters = "-Xmx4000m")
pacman::p_load(ggplot2, GreedyExperimentalDesign)
NUM_CORES = 2

n = 100
w_0 = c(rep(1, n / 2), rep(-1, n/ 2))
p = 1

sigma_x = 1
mu_x = 0
mu_z = 0
sigma_zs = 10^seq(0.5, 2, by = 0.25) 

beta_T = 1
beta_0 = 1
bbeta = rep(1, p)

set.seed(0) #NOTE: the whole simulation is conditional on this one x
X = matrix(rnorm(n * p, mu_x, sigma_x), ncol = p)
X = X[order(X[,1 ]), , drop = FALSE] #we do this so the matching algorithm is easy: first with second, third with fourth, etc.
Sinv = solve(var(X))

N_z_draws = 1000
N_avg = 500
N_rand_test = 200 #must be strictly less than N_avg!
r_approx_opt_greedy_pair_sw = 10000 #must be strictly greater than N_avg!
alpha = 0.05

#we will be comparing three designs: CRFB, Matching, Optimal.
#For "optimal", we will run the greedy pair switch and take the best N_rand_test vectors
#For each sigma_z...
#For each design, 
#  For N_z_draws 
#     we will first simulate a z and then measure the average power over N_avg iterations.
#     For N_avg iterations
#       We draw a w vector from the design.
#       We compute the betaThat.
#       For N_rand_test randomizations, 
#          we compute a test statistic
#       We find the null distribution and calculate the p-value
#       If less than alpha, we tally a "find".
#     The proportion of finds is the power for this iteration
# The average power is computer for the design by averaging over the z_draws.

designs = c("CRFB", "matching", "more_optimal", "kallus_optimal")

all_results_by_sigma_z = list()

for (sigma_z in sigma_zs){
  
  W_z_properties = list()
  power_over_z_draws = matrix(NA, nrow = N_z_draws, ncol = length(designs))
  colnames(power_over_z_draws) = designs
  
  for (i_z in 1 : N_z_draws){
    cat("i_z:", i_z, "sigma_z", sigma_z, "\n")
    #draw one z
    z = rnorm(n, 0, sigma_z)
    
    #we now generate a whole bunch of w's specific to the design
    null_rejections = matrix(NA, nrow = N_avg, ncol = length(designs))
    colnames(null_rejections) = designs
    
    W_z_properties[[i_z]] = list()
    
    for (design in designs){
      W = matrix(NA, nrow = n, ncol = N_avg)
      if (design == "CRFB"){
        for (i_avg in 1 : N_avg){
          W[, i_avg] = sample(w_0)
        }
        W0 = W
        W0[W0 == -1] = 0
      } else if (design == "matching"){
        for (i_avg in 1 : N_avg){
          for (i_w in seq(2, n, by = 2)){
            W[c(i_w - 1, i_w), i_avg] = sample(c(-1, 1))
          }
        }
        W0 = W
        W0[W0 == -1] = 0
      } else if (design == "more_optimal"){
        rd = initGreedyExperimentalDesignObject(X, r_approx_opt_greedy_pair_sw, wait = TRUE, objective = "mahal_dist", num_cores = NUM_CORES)
        res = resultsGreedySearch(rd, max_vectors = N_avg)
        W0 = res$ending_indicTs[, 1 : N_avg]
        W = W0
        W[W == 0] = -1
        
        # mahal_dists = apply(W0, 2, function(w){compute_objective_val(X, w, objective = "mahal_dist")})
        # hist(log10(mahal_dists), br = 1000)
        # mean(mahal_dists)
        # sum(eigen(var(t(W)))$values^2)
      } else if (design == "kallus_optimal"){
        rd = initGreedyExperimentalDesignObject(X, r_approx_opt_greedy_pair_sw, wait = TRUE, objective = "mahal_dist", num_cores = NUM_CORES)
        res = resultsGreedySearch(rd, max_vectors = 1)
        W0 = res$ending_indicTs[, 1, drop = FALSE] #get the one best vector
        W = W0
        W[W == 0] = -1
      }
      
      #we now record properties of W and z

      W_z_properties[[i_z]][[design]] = list()
      if (design == "kallus_optimal"){
        sigma_W = W %*% t(W)
      } else {
        sigma_W = var(t(W))
      }
      
      lambdas = eigen(sigma_W)$values
      max_eigenval_sigma_W = max(lambdas)
      W_z_properties[[i_z]][[design]]$max_eigenval_sigma_W = max_eigenval_sigma_W
      W_z_properties[[i_z]][[design]]$frob_norm_sq = sum(lambdas^2)
      mahal_dists = apply(W0, 2, function(w){compute_objective_val(X, w, objective = "mahal_dist")})      
      W_z_properties[[i_z]][[design]]$mean_mahal = mean(mahal_dists)
      W_z_properties[[i_z]][[design]]$acc_bias = t(z) %*% sigma_W %*% z
      W_z_properties[[i_z]][[design]]$efron_bound = sum(z^2) * max_eigenval_sigma_W
      
      for (i_avg in 1 : N_avg){
        #get the experimental vector
        if (design == "kallus_optimal"){
          w_experimental = W[, 1]
        } else {
          w_experimental = W[, i_avg]
        }
        
        #"run" the experiment
        y = beta_0 + X %*% bbeta + w_experimental * beta_T + z
        
        yT = y[w_experimental == 1]
        yC = y[w_experimental == -1]
        
        #compute our estimator
        beta_hat_exp = (mean(yT) - mean(yC)) / 2
        
        if (design == "kallus_optimal"){
          #run the bootstrap test
          beta_hat_boot = array(NA, N_rand_test)
          for (i_rand in 1 : N_rand_test){
            #draw a bootstrap sample
            is_boot = sample(1 : n, n, replace = TRUE)
            Xboot = X[is_boot, , drop = FALSE]
            yboot = y[is_boot]
            
            rd = initGreedyExperimentalDesignObject(Xboot, r_approx_opt_greedy_pair_sw, wait = TRUE, objective = "mahal_dist", num_cores = NUM_CORES)
            res = resultsGreedySearch(rd, max_vectors = 1)
            w_boot = res$ending_indicTs[, 1, drop = FALSE] #get the one best vector
            w_boot[w_boot == 0] = -1
            
            yT_boot = y[w_boot == 1]
            yC_boot = y[w_boot == -1]
            
            #compute our estimator
            beta_hat_boot[i_rand] = (mean(yT_boot) - mean(yC_boot)) / 2
          }
          #accept or reject based on alpha level?
          if (sum(abs(beta_hat_boot) > abs(beta_hat_exp)) / N_rand_test < alpha){
            null_rejections[i_avg, design] = 1
          } else {
            null_rejections[i_avg, design] = 0
          }
        } else {
          #now run the randomization test
          #first approximate the null distribution
          Wrand = W[, sample(setdiff(1 : N_avg, i_avg), N_rand_test)]
          beta_hat_null = array(NA, N_rand_test)
          for (i_rand in 1 : N_rand_test){
            w_rand = Wrand[, i_rand]
            
            yT_rand = y[w_rand == 1]
            yC_rand = y[w_rand == -1]
            
            #compute our estimator
            beta_hat_null[i_rand] = (mean(yT_rand) - mean(yC_rand)) / 2
          }
          #accept or reject based on alpha level?
          if (beta_hat_exp > quantile(beta_hat_null, 1 - alpha / 2) || beta_hat_exp < quantile(beta_hat_null, alpha / 2)){
            null_rejections[i_avg, design] = 1
          } else {
            null_rejections[i_avg, design] = 0
          }
        }
      }
    }
    
    power_over_z_draws[i_z, ] = colMeans(null_rejections)
    cat("iteration\n")
    print(colMeans(null_rejections, na.rm = TRUE))
    print(W_z_properties[[i_z]])
    cat("running average\n")
    print(colMeans(power_over_z_draws, na.rm = TRUE))
    cat("running 5% quantile\n")
    print(apply(power_over_z_draws, 2, quantile, 0.05, na.rm = TRUE))
  }
  
  key = as.character(round(sigma_z, 2))
  all_results_by_sigma_z[[key]] = list()  
  all_results_by_sigma_z[[key]]$power = colMeans(power_over_z_draws)
  all_results_by_sigma_z[[key]]$power_over_z_draws = power_over_z_draws
  all_results_by_sigma_z[[key]]$W_z_properties = W_z_properties
  #in case computer restarts...
  save(all_results_by_sigma_z, file = "all_results_by_sigma_z.RData")
}


###################################