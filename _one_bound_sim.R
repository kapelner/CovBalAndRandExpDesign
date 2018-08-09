
#simulate x above
z = rnorm(n, 0, sigma_z)


if (z_orthog){
  z = z - (x %*% z) / sum(x^2) * x
}

z_norm_sq = sum(z^2)

#under complete randomization with forced balance, here all all possible randomizations
all_randomizations = uniqueperm2(rep(c(-1, 1), n / 2))
w_size = nrow(all_randomizations)

#create place for simulation results to be stored and appropriately name it
res_iter = data.frame(matrix(NA, nrow = w_size, ncol = 6))
colnames(res_iter) = c("i", "obs_imbalance", "unobs_imbalance", "tx_estimate", "Rsq", "varY")

for (i in 1 : w_size){
  indicT = all_randomizations[i, ]
  
  t_idx = indicT == 1
  c_idx = indicT == -1
  xT = x[t_idx]
  xC = x[c_idx]
  zT = z[t_idx]
  zC = z[c_idx]
  
  y = x + z + indicT * beta_T
  yT = y[t_idx]
  yC = y[c_idx]
  res_iter[i, ] = c(
    i,
    (mean(xT) - mean(xC))^2, 
    (mean(zT) - mean(zC))^2, 
    (mean(yT) - mean(yC)) / 2,
    summary(lm(y ~ x))$r.squared,
    var(y)
  )
}
# 
# sqrt(mean(res_iter$Rsq))

#now let's order by observed imbalance
res_iter_ord_obs_imb = res_iter[order(res_iter$obs_imbalance), ]
# res_iter_ord_unobs_imb = res_iter[order(res_iter$unobs_imbalance), ]
# ggplot(res_iter_ord_obs_imb) + 
#   geom_histogram(aes(obs_imbalance), fill = "red", alpha = 0.5) + 
#   geom_histogram(aes(unobs_imbalance), fill = "blue", alpha = 0.5)
# 
# ggplot(res_iter_ord_obs_imb) + 
#   geom_point(aes(x = 1 : w_size, y = obs_imbalance), col = "red") + 
#   geom_point(aes(x = 1 : w_size, y = unobs_imbalance), col = "blue", lwd = 0.1, alpha = 0.3) + 
#   geom_smooth(aes(x = 1 : w_size, y = unobs_imbalance), col = "blue")

# g = ggplot(res_iter_ord_obs_imb, aes(x = obs_imbalance, y = unobs_imbalance)) + 
#   geom_point(lwd = 0.1, alpha = 0.3) +
#   # stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE) +
#   xlab("observed imbalance") + ylab("unobserved imbalance")
# # plot(g)
# ggMarginal(g, type = "histogram", binwidth = 0.02)

# ggplot(res_iter_ord_unobs_imb) + 
#   geom_point(aes(x = 1 : w_size, y = unobs_imbalance), col = "blue") + 
#   geom_smooth(aes(x = 1 : w_size, y = obs_imbalance), col = "red")

#now let's calculate mse and lambda_max
res_iter_ord_obs_imb$mse = array(NA, w_size)
res_iter_ord_obs_imb$true_variance = array(NA, w_size)
res_iter_ord_obs_imb$lambda_max = array(NA, w_size)
# res_iter_ord_obs_imb$var_obj_new = array(NA, w_size)
# res_iter_ord_obs_imb$var_obj_x_term = array(NA, w_size)
# res_iter_ord_obs_imb$var_obj_x_term_prop = array(NA, w_size)
res_iter_ord_obs_imb$var_obj_bound = array(NA, w_size)
res_iter_ord_obs_imb$var_from_observed_cov = array(NA, w_size)
res_iter_ord_obs_imb$var_bound_from_unobserved_cov = array(NA, w_size)
# res_iter_ord_unobs_imb$lambda_max = array(NA, w_size)
# res_iter_ord_unobs_imb$var_obj = array(NA, w_size)
# res_iter_ord_unobs_imb$var_obj_bound = array(NA, w_size)
sigma_w_obs_imb = list()
# sigma_w_obs_unimb = list()

for (i in seq(from = 2, to = w_size)){
  res_iter_eff = res_iter_ord_obs_imb[1 : i, ]
  res_iter_ord_obs_imb$mse[i] = mean((res_iter_eff$tx_estimate - beta_T)^2)
  sigma_w = var(all_randomizations[res_iter_eff$i, ]) * (i - 1) / i
  if (KEEP_SIGMA){
    sigma_w_obs_imb[[i]] = sigma_w
  }
  
  # res_iter_ord_obs_imb$true_variance[i] = (1 / n^2) * 
    # (t(x) %*% sigma_w %*% (x) + 2 * t(x) %*% sigma_w %*% (z) + t(z) %*% sigma_w %*% (z))
  # res_iter_ord_obs_imb$var_obj_new[i] = t(x) %*% sigma_w %*% (x) + 2 * t(x) %*% sigma_w %*% (z) + t(z) %*% sigma_w %*% (z)
  # res_iter_ord_obs_imb$var_obj_x_term[i] = 2 * t(x) %*% sigma_w %*% (z)
  # res_iter_ord_obs_imb$var_obj_x_term_prop[i] = t(x) %*% sigma_w %*% (z) / t(z) %*% sigma_w %*% (z)
  lambda_max = max(eigen(sigma_w)$values)
  res_iter_ord_obs_imb$var_obj_bound[i] = (1 / n^2) * (t(x) %*% sigma_w %*% (x) + n * sigma_z^2 * lambda_max)
  # res_iter_ord_obs_imb$var_from_observed_cov[i] = t(x) %*% sigma_w %*% (x)
  # res_iter_ord_obs_imb$var_bound_from_unobserved_cov[i] = n * sigma_z^2 * lambda_max
  res_iter_ord_obs_imb$lambda_max[i] = lambda_max
  
  # res_iter_eff = res_iter_ord_unobs_imb[1 : i, ]
  # res_iter_ord_unobs_imb$mse[i] = mean((res_iter_eff$tx_estimate - 2 * betaT_over_two)^2)
  # sigma_w = var(all_randomizations[res_iter_eff$i, ]) * (i - 1) / i
  # sigma_w_obs_unimb[[i]] = sigma_w
  # res_iter_ord_unobs_imb$var_obj[i] = t(x) %*% sigma_w %*% (x) + t(z) %*% sigma_w %*% (z)
  # lambda_max = max(eigen(sigma_w)$values)
  # res_iter_ord_unobs_imb$var_obj_bound[i] = t(x) %*% sigma_w %*% (x) + z_norm_sq * lambda_max 
  # res_iter_ord_unobs_imb$lambda_max[i] = lambda_max
}