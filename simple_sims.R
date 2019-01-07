source("unique_perm.R")
pacman::p_load(ggplot2)

n = 20
p = 1

sigma_x = 1
mu_x = 0
sigma_z = 1.5 ###change this to give Figures 1,2,3 and comment our the darkorange
# sigma_z = 1
# sigma_z = 0.05 #R^2 = 1
mu_z = 0

beta_T = 1
beta_0 = 1
bbeta = rep(1, p)

set.seed(0) #NOTE: the whole simulation is conditional on this one x
X = matrix(rnorm(n * p, mu_x, sigma_x), ncol = p)
X = X[order(X[,1 ]), , drop = FALSE]
Sinv = solve(var(X))


#get all possible realizations
all_randomizations = uniqueperm2(rep(c(-1, 1), n / 2))
w_size = nrow(all_randomizations)



res_iter = data.frame(matrix(NA, nrow = w_size, ncol = 2))
colnames(res_iter) = c("i", "obs_imbalance")

for (i in 1 : w_size){
  indicT = all_randomizations[i, ]
  t_idx = indicT == 1
  c_idx = indicT == -1
  xT = X[t_idx]
  xC = X[c_idx]
  res_iter[i, ] = c(
    i,
    (mean(xT) - mean(xC)) %*% Sinv %*% (mean(xT) - mean(xC))
  )
}

#now let's order by observed imbalance
res_iter_ord_obs_imb = res_iter[order(res_iter$obs_imbalance), ]

###############START SIM

Nsim = 1000
Nsamprand = 300

#create place for simulation results to be stored and appropriately name it
colnames_results = c("mse_naive", "mse_regr")
rand_res_iter_obs = data.frame(matrix(NA, nrow = Nsim, ncol = length(colnames_results)))
colnames(rand_res_iter_obs) = colnames_results
match_res_iter_obs = data.frame(matrix(NA, nrow = Nsim, ncol = length(colnames_results)))
colnames(match_res_iter_obs) = colnames_results
opt_res_iter_obs = data.frame(matrix(NA, nrow = Nsim, ncol = 2))
colnames(opt_res_iter_obs) = colnames_results
worst_res_iter_obs = data.frame(matrix(NA, nrow = Nsim, ncol = 2))
colnames(worst_res_iter_obs) = colnames_results

r_sq_est = array(NA, Nsim)
for (nsim in 1 : Nsim){
  if (nsim %% 100 == 0){
    cat("random nsim: ", nsim, "\n")
  }
  
  #simulate the unobserved features
  z = rnorm(n, 0, sigma_z)
  
  #now sample for random
  tx_est = array(NA, Nsamprand)
  tx_est_regr = array(NA, Nsamprand)
  r_sq_est_int = array(NA, Nsamprand)
  for (i in 1 : Nsamprand){
    indicT = all_randomizations[sample(1 : w_size, 1), ]
    
    t_idx = indicT == 1
    c_idx = indicT == -1
    xT = X[t_idx]
    xC = X[c_idx]
    zT = z[t_idx]
    zC = z[c_idx]
    
    y = beta_0 + X %*% bbeta + z + indicT * beta_T
    # y = beta_0 + z + indicT * beta_T
    
    yT = y[t_idx]
    yC = y[c_idx]
    
    tx_est[i] = (mean(yT) - mean(yC)) / 2
    tx_est_regr[i] = coef(lm(y ~ X + indicT))[3]
    r_sq_est_int[i] = summary(lm(y ~ z))$r.squared
  }
  
  
  rand_res_iter_obs[nsim, ] = c(
    mean((tx_est - beta_T)^2),
    mean((tx_est_regr - beta_T)^2)
  )
  r_sq_est[nsim] = mean(r_sq_est_int)
}
mean(r_sq_est)

sigma_w_crfb = (1 + 1 / (n - 1)) * diag(n) - matrix(rep(1 / (n - 1), n^2), nrow = n)
t(X) %*% sigma_w_crfb %*% X / n^2
sigma_z^2 / n
t(X) %*% sigma_w_crfb %*% X / n^2 + sigma_z^2 / n

mean(rand_res_iter_obs$mse_naive)


for (nsim in 1 : Nsim){
  if (nsim %% 100 == 0){
    cat("matching nsim: ", nsim, "\n")
  }
  
  #simulate the unobserved features outside of the w loop
  z = rnorm(n, 0, sigma_z)

  tx_est = array(NA, Nsamprand)
  tx_est_regr = array(NA, Nsamprand)
  for (i in 1 : Nsamprand){ 
    
    #now sample for pairwise matching
    indicT = matrix(NA, n, 1)
    for (i_w in seq(2, n, by = 2)){
      indicT[c(i_w - 1, i_w), 1] = sample(c(-1, 1))
    }
    
    t_idx = indicT == 1
    c_idx = indicT == -1
    xT = X[t_idx]
    xC = X[c_idx]
    zT = z[t_idx]
    zC = z[c_idx]
    
    y = beta_0 + X %*% bbeta + z + indicT * beta_T
    # y = beta_0 + z + indicT * beta_T
    
    yT = y[t_idx]
    yC = y[c_idx]
    
    tx_est[i] = (mean(yT) - mean(yC)) / 2
    tx_est_regr[i] = coef(lm(y ~ X + indicT))[3]
  }
  
  match_res_iter_obs[nsim, ] = c(
    mean((tx_est - beta_T)^2),
    mean((tx_est_regr - beta_T)^2)
  )
}


sigma_w_matching = diag(n)
for (i in seq(from = 2, to = n, by = 2)){
  sigma_w_matching[i - 1, i] = -1
  sigma_w_matching[i, i - 1] = -1
}

t(X) %*% sigma_w_matching %*% X / n^2
sigma_z^2 / n
t(X) %*% sigma_w_matching %*% X / n^2 + sigma_z^2 / n

mean(match_res_iter_obs$mse_naive)

for (nsim in 1 : Nsim){
  if (nsim %% 100 == 0){
    cat("opt nsim: ", nsim, "\n")
  }
  
  #simulate the unobserved features
  z = rnorm(n, 0, sigma_z)
  
  #now sample for optimal
  tx_est = array(NA, 2)
  tx_est_regr = array(NA, 2)
  for (i in 1 : 2){ #there are only two optimal vectors!
    indicT = all_randomizations[res_iter_ord_obs_imb$i[i], ]
    
    t_idx = indicT == 1
    c_idx = indicT == -1
    xT = X[t_idx]
    xC = X[c_idx]
    zT = z[t_idx]
    zC = z[c_idx]
    
    y = beta_0 + X %*% bbeta + z + indicT * beta_T
    # y = beta_0 + z + indicT * beta_T
    
    yT = y[t_idx]
    yC = y[c_idx]
    
    tx_est[i] = (mean(yT) - mean(yC)) / 2
    tx_est_regr[i] = coef(lm(y ~ X + indicT))[3]
  } 
  
  opt_res_iter_obs[nsim, ] = c(
    mean((tx_est - beta_T)^2),
    mean((tx_est_regr - beta_T)^2)
  )
}



# w_star = all_randomizations[res_iter_ord_obs_imb$i[1], ]
# sigma_w_opt = w_star %*% t(w_star)
# 
# t(X) %*% sigma_w_opt %*% X / n^2
# sigma_z^2 / n
# t(X) %*% sigma_w_opt %*% X / n^2 + sigma_z^2 / n

#mean(opt_res_iter_obs$mse_naive)

r_worst = 500
res_iter_ord_obs_imb_worst = res_iter_ord_obs_imb[(w_size - r_worst + 1) : w_size, ]

for (nsim in 1 : Nsim){
  if (nsim %% 100 == 0){
    cat("worst nsim: ", nsim, "\n")
  }
  
  #simulate the unobserved features
  z = rnorm(n, 0, sigma_z)
  
  #now sample for worst
  tx_est = array(NA, r_worst)
  tx_est_regr = array(NA, r_worst)
  for (i in 1 : r_worst){
    indicT = all_randomizations[res_iter_ord_obs_imb_worst$i[i], ]
    
    t_idx = indicT == 1
    c_idx = indicT == -1
    xT = X[t_idx]
    xC = X[c_idx]
    zT = z[t_idx]
    zC = z[c_idx]
    
    y = beta_0 + X %*% bbeta + z + indicT * beta_T
    
    yT = y[t_idx]
    yC = y[c_idx]
    
    tx_est[i] = (mean(yT) - mean(yC)) / 2
    tx_est_regr[i] = coef(lm(y ~ X + indicT))[3]
  } 
  
  worst_res_iter_obs[nsim, ] = c(
    mean((tx_est - beta_T)^2),
    mean((tx_est_regr - beta_T)^2)
  )
}


#what happened?
mean(rand_res_iter_obs$mse_naive, na.rm = TRUE)
mean(match_res_iter_obs$mse_naive, na.rm = TRUE)
mean(opt_res_iter_obs$mse_naive, na.rm = TRUE)
mean(worst_res_iter_obs$mse_naive)
quantile(rand_res_iter_obs$mse_naive, 0.95, na.rm = TRUE)
quantile(match_res_iter_obs$mse_naive, 0.95, na.rm = TRUE)
quantile(opt_res_iter_obs$mse_naive, 0.95, na.rm = TRUE)
quantile(worst_res_iter_obs$mse_naive, 0.95, na.rm = TRUE)
#calculate the c constants
(quantile(rand_res_iter_obs$mse_naive, 0.95, na.rm = TRUE) - mean(rand_res_iter_obs$mse_naive, na.rm = TRUE)) / sd(rand_res_iter_obs$mse_naive, na.rm = TRUE)
(quantile(match_res_iter_obs$mse_naive, 0.95, na.rm = TRUE) - mean(match_res_iter_obs$mse_naive, na.rm = TRUE)) / sd(match_res_iter_obs$mse_naive, na.rm = TRUE)
(quantile(opt_res_iter_obs$mse_naive, 0.95, na.rm = TRUE) - mean(opt_res_iter_obs$mse_naive, na.rm = TRUE)) / sd(opt_res_iter_obs$mse_naive, na.rm = TRUE)
(quantile(worst_res_iter_obs$mse_naive, 0.95, na.rm = TRUE) - mean(worst_res_iter_obs$mse_naive, na.rm = TRUE)) / sd(opt_res_iter_obs$mse_naive, na.rm = TRUE)


ggplot(data.frame(rand_res_iter_obs)) + 
  geom_density(aes(mse_naive), alpha = 0.3, fill = "red") + 
  geom_density(aes(mse_naive), alpha = 0.3, fill = "blue", data = match_res_iter_obs) +
  geom_density(aes(mse_naive), alpha = 0.3, fill = "green", data = opt_res_iter_obs) + 
  geom_density(aes(mse_naive), alpha = 0.3, fill = "darkorange", data = worst_res_iter_obs) + 
  xlim(0, 0.6) + xlab("MSE") +
  geom_vline(xintercept = mean(rand_res_iter_obs$mse_naive), col = "red", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = mean(match_res_iter_obs$mse_naive), col = "blue", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = mean(opt_res_iter_obs$mse_naive), col = "green", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = mean(worst_res_iter_obs$mse_naive), col = "darkorange", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = quantile(rand_res_iter_obs$mse_naive, .95), col = "red", alpha = 0.3, lwd = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(match_res_iter_obs$mse_naive, .95), col = "blue", alpha = 0.3, lwd = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(opt_res_iter_obs$mse_naive, .95), col = "green", alpha = 0.3, lwd = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(worst_res_iter_obs$mse_naive, .95), col = "darkorange", alpha = 0.3, lwd = 1, linetype = "dashed")

max(rand_res_iter_obs$mse_naive)
max(match_res_iter_obs$mse_naive)
max(opt_res_iter_obs$mse_naive)
max(worst_res_iter_obs$mse_naive)

#conclusion: matching wins

### investigate regression estimator


mean(rand_res_iter_obs$mse_regr, na.rm = TRUE)
mean(match_res_iter_obs$mse_regr, na.rm = TRUE)
mean(opt_res_iter_obs$mse_regr, na.rm = TRUE)
mean(worst_res_iter_obs$mse_regr)
quantile(rand_res_iter_obs$mse_regr, 0.95, na.rm = TRUE)
quantile(match_res_iter_obs$mse_regr, 0.95, na.rm = TRUE)
quantile(opt_res_iter_obs$mse_regr, 0.95, na.rm = TRUE)
quantile(worst_res_iter_obs$mse_regr, 0.95, na.rm = TRUE)


ggplot(data.frame(rand_res_iter_obs)) + 
  geom_density(aes(mse_regr), alpha = 0.3, fill = "red") + 
  geom_density(aes(mse_regr), alpha = 0.3, fill = "blue", data = match_res_iter_obs) +
  geom_density(aes(mse_regr), alpha = 0.3, fill = "green", data = opt_res_iter_obs) + 
  geom_density(aes(mse_regr), alpha = 0.3, fill = "darkorange", data = worst_res_iter_obs) + 
  xlim(0, 0.6) + xlab("MSE") +
  geom_vline(xintercept = mean(rand_res_iter_obs$mse_regr), col = "red", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = mean(match_res_iter_obs$mse_regr), col = "blue", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = mean(opt_res_iter_obs$mse_regr), col = "green", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = mean(worst_res_iter_obs$mse_regr), col = "darkorange", alpha = 0.3, lwd = 1) +
  geom_vline(xintercept = quantile(rand_res_iter_obs$mse_regr, .95), col = "red", alpha = 0.3, lwd = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(match_res_iter_obs$mse_regr, .95), col = "blue", alpha = 0.3, lwd = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(opt_res_iter_obs$mse_regr, .95), col = "green", alpha = 0.3, lwd = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(worst_res_iter_obs$mse_regr, .95), col = "darkorange", alpha = 0.3, lwd = 1, linetype = "dashed")

max(rand_res_iter_obs$mse_regr)
max(match_res_iter_obs$mse_regr)
max(opt_res_iter_obs$mse_regr)
max(worst_res_iter_obs$mse_regr)
