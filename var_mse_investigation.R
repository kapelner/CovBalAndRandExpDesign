rm(list = ls())


pacman::p_load(ggplot2)

set.seed(1821)

##Section 2.1
z_orthog = FALSE
n = 16
w_size = choose(n, n / 2)
beta_T = 1

sigma_x = 1
sigma_z = 1
KEEP_SIGMA = FALSE
x = rnorm(n, 0, sigma_x)

source("unique_perm.R")

Nsim = 1000

res_all_sims = matrix(NA, nrow = Nsim, ncol = w_size)
all_Rsqs = array(NA, Nsim)

for (nsim in 1 : Nsim){
  cat("nsim:", nsim, "\n")
  source("_one_bound_sim.R")
  colnames(res_all_sims) = res_iter_ord_obs_imb$lambda_max
  g = ggplot(res_iter_ord_obs_imb) + 
    geom_point(aes(x = lambda_max, y = mse)) + 
    # geom_point(aes(x = lambda_max, y = true_variance), col = "blue") +
    geom_point(aes(x = lambda_max, y = var_obj_bound), col = "red") +
    xlab("Maximum Eigenvalue") + ylab("MSE")
  plot(g)
  
  res_all_sims[nsim, ] = res_iter_ord_obs_imb$mse
  print(res_all_sims[1 : nsim, c(2, 248)])
  all_Rsqs[nsim] = mean(res_iter$Rsq)
  cat("running avg Rsq of y ~ x:", mean(all_Rsqs, na.rm = T), "\n")
  cat("mse avgs: ", colMeans(res_all_sims[1 : nsim, c(2, 248), drop = FALSE]), "\n")
  cat("mse sds: ", apply(res_all_sims[1 : nsim, c(2, 248), drop = FALSE], 2, sd), "\n")
}

g = ggplot(res_iter_ord_obs_imb) + 
  # geom_point(aes(x = lambda_max, y = true_variance), col = "blue") +
  # geom_line(aes(x = lambda_max, y = var_obj_bound), col = "red", lwd = 1.3) +
  xlab("Maximum Eigenvalue") + ylab("MSE")
for (nsim in 1 : 400){
  res_all_sims_nsim = data.frame(
    lambda_max = res_iter_ord_obs_imb$lambda_max[2 : w_size],
    mse = res_all_sims[nsim, 2 : w_size]
  )
  g = g + geom_line(aes(x = lambda_max, y = mse), lwd = 0.2, alpha = 0.3, data = res_all_sims_nsim)
}
#now plot the avg
res_all_sims_avg = data.frame(
  lambda_max = res_iter_ord_obs_imb$lambda_max[2 : w_size],
  mse = colMeans(res_all_sims)[2 : w_size]
)
g = g + geom_line(aes(x = lambda_max, y = mse), data = res_all_sims_avg, col = "green", lwd = 1.3)
#now plot the 95%ile
res_all_sims_quantile = data.frame(
  lambda_max = res_iter_ord_obs_imb$lambda_max[2 : w_size],
  qu = apply(res_all_sims[, 2 : w_size], 2, quantile, 0.95)
)
g = g + geom_line(aes(x = lambda_max, y = qu), data = res_all_sims_quantile, col = "blue", lwd = 1.3)
g = g + ylim(0, 0.52)
g = g + geom_line(aes(x = lambda_max, y = var_obj_bound), data = res_iter_ord_obs_imb, col = "red", lwd = 1.3)
plot(g)


#where is the kink in the avg MSE?
ggplot(res_all_sims_avg[res_all_sims_avg$lambda_max < 3, ]) + geom_line(aes(x = lambda_max, y = mse))




# ggplot(data.frame(res_all_sims)) + geom_line(aes(x = res_iter_ord_obs_imb$lambda_max, y = res_iter_ord_obs_imb$var_obj_bound), col = "red")
# ggplot(data.frame(order_num = 100 : (w_size - 1), lambda_max_lag_1_diff = diff(res_iter_ord_obs_imb$lambda_max)[100 : (w_size - 1)])) + geom_point(aes(x = order_num, y = lambda_max_lag_1_diff))
