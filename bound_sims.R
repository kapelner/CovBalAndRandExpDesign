rm(list = ls())


pacman::p_load(xtable, ggplot2, ggExtra, latex2exp)

set.seed(1821)

##Section 2.1
z_orthog = FALSE
n = 16
beta_T = 1

sigma_x = 1
sigma_z = 1
KEEP_SIGMA = TRUE

x = rnorm(n, 0, sigma_x)

source("unique_perm.R")
source("_one_bound_sim.R")

w_start = 2
w_end = 3000



#checks to make sure our variance expression is correct
ggplot(res_iter_ord_obs_imb) + 
  geom_point(aes(x = lambda_max, y = true_variance), col = "yellow") + 
  geom_point(aes(x = lambda_max, y = mse), col = "blue", alpha = 0.4, lwd = 0.3) +
  xlab(TeX("$\\lambda_{max}$")) + ylab("MSE & true variance")


# ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = w_start : w_size, y = mse)) + xlab("in increasing obs imb")
# ggplot(res_iter_ord_obs_imb[w_start : w_end, ]) + 
#   geom_point(aes(x = w_start : w_end, y = mse)) + xlab("in increasing obs imb")


# ggplot(res_iter_ord_unobs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = w_start : w_size, y = mse)) + xlab("in increasing unobs imb")

# ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = w_start : w_size, y = lambda_max)) + xlab("in increasing obs imb")
# ggplot(res_iter_ord_unobs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = w_start : w_size, y = lambda_max)) + xlab("in increasing unobs imb")


ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
  geom_point(aes(x = obs_imbalance, y = mse))
ggplot(res_iter_ord_obs_imb[w_start : w_end, ]) + 
  geom_point(aes(x = obs_imbalance, y = mse))

ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
  geom_point(aes(x = obs_imbalance, y = lambda_max))
# ggplot(res_iter_ord_unobs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = unobs_imbalance, y = lambda_max))

ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
  geom_point(aes(x = lambda_max, y = mse))

# ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
#   geom_smooth(aes(x = lambda_max, y = mse))

# ggplot(res_iter_ord_unobs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = lambda_max, y = mse))

res_iter_ord_obs_imb[res_iter_ord_obs_imb$var_obj_bound == min(res_iter_ord_obs_imb$var_obj_bound, na.rm = TRUE), ]$i
i_star = res_iter_ord_obs_imb[res_iter_ord_obs_imb$var_obj_bound == min(res_iter_ord_obs_imb$var_obj_bound, na.rm = TRUE), ]$i
round(sigma_w_obs_imb[[i_star[2]]], 2)
#find largest entry of the varcov matrix that is not a one
max_entry = max(abs(sigma_w_obs_imb[[i_star[2]]] - diag(n)))
max_entry
entries = which(abs(sigma_w_obs_imb[[i_star[2]]] - diag(n)) == max_entry, arr.ind = TRUE)[1, ]
x[entries]
sort(x)

res_iter_ord_obs_imb[res_iter_ord_obs_imb$var_obj_bound == 
  min(res_iter_ord_obs_imb$var_obj_bound, na.rm = TRUE), ]$lambda_max
res_iter_ord_obs_imb[res_iter_ord_obs_imb$var_obj_bound == 
                       min(res_iter_ord_obs_imb$var_obj_bound, na.rm = TRUE), ]$i
res_iter_ord_obs_imb[res_iter_ord_obs_imb$mse == 
  min(res_iter_ord_obs_imb$mse, na.rm = TRUE), ]$lambda_max

# ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = lambda_max, y = var_obj_x_term))
# ggplot(res_iter_ord_obs_imb[w_start : 100, ]) + 
#   geom_point(aes(x = lambda_max, y = var_obj_x_term))
# 
# 
# ggplot(res_iter_ord_obs_imb[w_start : w_size, ]) + 
#   geom_point(aes(x = lambda_max, y = var_obj_x_term_prop))
# ggplot(res_iter_ord_obs_imb[w_start : 100, ]) + 
#   geom_point(aes(x = lambda_max, y = var_obj_x_term_prop))


ggplot(res_iter_ord_obs_imb) + 
  geom_point(aes(x = lambda_max, y = mse)) + 
  # geom_point(aes(x = lambda_max, y = true_variance), col = "blue") +
  geom_point(aes(x = lambda_max, y = var_obj_bound), col = "red") +
  xlab("Maximum Eigenvalue") + ylab("MSE")


ggplot(res_iter_ord_obs_imb[res_iter_ord_obs_imb$lambda_max < 1.1, ]) +
  geom_point(aes(x = lambda_max, y = mse)) +
  geom_point(aes(x = lambda_max, y = var_obj_bound), col = "red") +
  xlab("Maximum Eigenvalue") + ylab("")


#how do the pieces of the variance bound contribute?
ggplot(res_iter_ord_obs_imb) +
  geom_point(aes(x = lambda_max, y = var_from_observed_cov), col = "darkgreen") +
  geom_point(aes(x = lambda_max, y = var_bound_from_unobserved_cov), col = "brown") +
  xlab("Maximum Eigenvalue") + ylab("Variance Contribution")

ggplot(res_iter_ord_obs_imb) +
  geom_point(aes(x = lambda_max, y = var_from_observed_cov), col = "darkgreen") +
  geom_point(aes(x = lambda_max, y = var_bound_from_unobserved_cov), col = "brown") +
  xlab("Maximum Eigenvalue") + ylab("Variance Contribution")






