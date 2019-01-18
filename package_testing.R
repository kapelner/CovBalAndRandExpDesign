library(OptimalRerandExpDesigns)

set.seed(1984)
n = 100
p = 10
X = matrix(rnorm(n * p), nrow = n, ncol = p)
X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
max_designs = 25000

W_base_obj = generate_W_base_and_sort(X, max_designs = max_designs)
W_base_obj
#Fig 1
plot(W_base_obj, title = "", xlab = "Mahalanobis Distance")

#Fig 2
rerand_res_norm = optimal_rerandomization_normality_assumed(W_base_obj)
rerand_res_norm
plot(rerand_res_norm)

rerand_res_ex = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rnorm(n)}, N_z = 10)
rerand_res_ex
plot(rerand_res_ex)




rerand_res_approx = optimal_rerandomization_tail_approx(W_base_obj)
rerand_res_approx
plot(rerand_res_approx)


QA1 = rerand_res_norm$all_data_from_run$Q_primes
QA2 = rerand_res_ex$all_data_from_run$Q_primes
QA3 = rerand_res_approx$all_data_from_run$Q_primes
QA1 = QA1 / QA1[max_designs] #normalize starting position
QA2 = QA2 / QA2[max_designs] #normalize starting position
QA3 = QA3 / QA3[max_designs] #normalize starting position

ggplot(data.frame(QA1 = QA1, QA2 = QA2, QA3 = QA3)) +
  scale_x_log10() +
  scale_y_log10() + 
  xlab("rerandomization threshold (in mahalanobis distance)") +
  ylab("tail criterion (relative)") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA2), col = "red") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA1), col = "blue") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA3), col = "green") +
  geom_vline(xintercept = rerand_res_ex$a_star, col = "red") +
  geom_vline(xintercept = rerand_res_norm$a_star, col = "blue") +
  geom_vline(xintercept = rerand_res_approx$a_star, col = "green")

#Fig 3a
pacman::p_load(rmutil)
rerand_res_ex_laplace = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rlaplace(n)}, N_z = 10)
rerand_res_ex_laplace
plot(rerand_res_ex_laplace)

rerand_res_approx_laplace_kurtosis = optimal_rerandomization_tail_approx(W_base_obj, excess_kurtosis_z = 3)
rerand_res_approx_laplace_kurtosis
plot(rerand_res_approx_laplace_kurtosis)

QA4 = rerand_res_approx_laplace_kurtosis$all_data_from_run$Q_primes
QA5 = rerand_res_ex_laplace$all_data_from_run$Q_primes
QA4 = QA4 / QA4[max_designs] #normalize starting position
QA5 = QA5 / QA5[max_designs] #normalize starting position


ggplot(data.frame(QA1 = QA1, QA3 = QA3, QA4 = QA4, QA5 = QA5)) +
  scale_x_log10() +
  scale_y_log10() + 
  xlab("rerandomization threshold (in mahalanobis distance)") +
  ylab("tail criterion (relative)") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA4), col = "red") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA1), col = "blue") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA3), col = "green") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA5), col = "purple") +
  geom_vline(xintercept = rerand_res_ex_laplace$a_star + 0.003, col = "red") +
  geom_vline(xintercept = rerand_res_norm$a_star, col = "blue") +
  geom_vline(xintercept = rerand_res_approx$a_star, col = "green") +
  geom_vline(xintercept = rerand_res_approx_laplace_kurtosis$a_star, col = "purple")

#Fig 3b
rerand_res_ex_t2 = optimal_rerandomization_exact(W_base_obj, z_sim_fun = function(){rt(n, 2)}, N_z = 100)
rerand_res_ex_t2
plot(rerand_res_ex_t2)

QA6 = rerand_res_ex_t2$all_data_from_run$Q_primes
QA6 = QA6 / QA6[max_designs] #normalize starting position


ggplot(data.frame(QA1 = QA1, QA3 = QA3, QA6 = QA6)) +
  scale_x_log10() +
  scale_y_log10() + 
  xlab("rerandomization threshold (in mahalanobis distance)") +
  ylab("tail criterion (relative)") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA4), col = "red") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA1), col = "blue") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA3), col = "green") +
  geom_line(aes(x = W_base_obj$imbalance_by_w_sorted, y = QA6), col = "purple") +
  geom_vline(xintercept = rerand_res_ex_t2$a_star + 0.003, col = "red") +
  geom_vline(xintercept = rerand_res_norm$a_star, col = "blue") +
  geom_vline(xintercept = rerand_res_approx$a_star, col = "green")

#Fig 3
set.seed(1984)
bbeta = rnorm(p)
betaT = 1
N_w = 1000
N_z = 1000
sigma_z = 3

mse_z_all = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 100 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(1 : max_designs, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_all[n_z] = mean(mse_w)
}
cat("\n")

mse_z_top_5000 = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 100 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(1 : 5000, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_top_5000[n_z] = mean(mse_w)
}
cat("\n")


mse_z_worst_5000 = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 100 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(20001 : 25000, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_worst_5000[n_z] = mean(mse_w)
}
cat("\n")


mse_z_opt = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 100 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = rerand_res_norm$W_star[sample(1 : rerand_res_norm$W_star_size, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_opt[n_z] = mean(mse_w)
}
cat("\n")


mse_z_determ = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 100 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  w = W_base_obj$W_base_sorted[1, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
  mse_z_determ[n_z] = (betaThatLinRegr - betaT)^2
}
cat("\n")

dif_threshold_designs = data.frame(
  mse_z_all = mse_z_all,
  mse_z_top_5000 = mse_z_top_5000,
  mse_z_worst_5000 = mse_z_worst_5000,
  mse_z_opt = mse_z_opt,
  mse_z_determ = mse_z_determ
)
alpha = 0.2
q = 0.95
ggplot(dif_threshold_designs) +
  geom_density(aes(x = mse_z_all),                       fill = "red", alpha = alpha) +
  geom_density(aes(x = mse_z_top_5000),                  fill = "pink", alpha = alpha) +
  geom_density(aes(x = mse_z_worst_5000),                fill = "black", alpha = alpha) +
  geom_density(aes(x = mse_z_opt),                       fill = "green", alpha = alpha) +
  geom_density(aes(x = mse_z_determ),                    fill = "blue", alpha = alpha) +
  geom_vline(xintercept = mean(mse_z_all),               col = "red", lwd = 1) +
  geom_vline(xintercept = mean(mse_z_top_5000),          col = "pink", lwd = 1) +
  geom_vline(xintercept = mean(mse_z_worst_5000),        col = "black", lwd = 1) +
  geom_vline(xintercept = mean(mse_z_opt),               col = "green", lwd = 1) +
  geom_vline(xintercept = mean(mse_z_determ),            col = "blue", lwd = 1) +
  geom_vline(xintercept = quantile(mse_z_all, q),        col = "red", linetype = "dashed", lwd = 1) +
  geom_vline(xintercept = quantile(mse_z_top_5000, q),   col = "pink", linetype = "dashed", lwd = 1) +
  geom_vline(xintercept = quantile(mse_z_worst_5000, q), col = "black", linetype = "dashed", lwd = 1) +
  geom_vline(xintercept = quantile(mse_z_opt, q),        col = "green", linetype = "dashed", lwd = 1) +
  geom_vline(xintercept = quantile(mse_z_determ, q),     col = "blue", linetype = "dashed", lwd = 1) +
  xlab("MSE for many different z") +
  xlim(0.02, 0.351)


#Fig 5
#a* by p
n = 100
set.seed(1984) #do this to go in order exactly
Xall = matrix(rnorm(n^2), nrow = n, ncol = n)
ps = seq(1, 91, by = 5)
max_designs = 25000

a_stars = array(NA, length(ps))
a_star_quantiles = array(NA, length(ps))
mahal_range = matrix(NA, nrow = length(ps), ncol = 2)
for (i_p in 1 : length(ps)){
  p = ps[i_p]
  cat("p = ", p, "\n")
  X = Xall[, 1 : p, drop = FALSE]
  X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
  W_base_obj_p = generate_W_base_and_sort(X, max_designs = max_designs)
  print(W_base_obj_p)
  mahal_range[i_p, ] = c(W_base_obj_p$imbalance_by_w_sorted[1], W_base_obj_p$imbalance_by_w_sorted[max_designs])
  
  rerand_res_norm_p = optimal_rerandomization_tail_approx(W_base_obj_p)
  print(rerand_res_norm_p)
  
  a_stars[i_p] = rerand_res_norm_p$a_star
  a_star_quantiles[i_p] = ecdf(W_base_obj_p$imbalance_by_w_sorted)(rerand_res_norm_p$a_star)
  

  gg = ggplot(data.frame(p = ps, a_star_quantiles = a_star_quantiles)) +
         geom_line(aes(x = ps, y = a_star_quantiles)) +
         xlab("p") + ylab("a* quantile of imbalances")
  plot(gg)
}

gg_full = ggplot(data.frame(p = ps, a_stars = a_stars)) +
       geom_line(aes(x = ps, y = a_stars)) +
       scale_x_log10() +
       scale_y_log10() +
       xlab("p") + ylab("a*")
for (i_p in 1 : length(ps)){
  gg_full = gg_full + 
    geom_segment(aes(x = ps[i_p], xend = ps[i_p], y = mahal_range[i_p, 1], yend = mahal_range[i_p, 2]), col = "darkgreen")
}
plot(gg_full)






#Experimental Power
alpha = 0.05
R = 1000
sigma_z = 10

rejections_all = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 1 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  idx_exp = sample(1 : max_designs, 1)
  w = W_base_obj$W_base_sorted[idx_exp, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr_exp = coef(lm(y ~ 0 + X + w))[p + 1]

  betaThatLinRegr_replicates = array(NA, R)
  for (r in 1 : R){
    w = W_base_obj$W_base_sorted[sample(setdiff(1 : max_designs, idx_exp), 1), ]
    betaThatLinRegr_replicates[r] = coef(lm(y ~ 0 + X + w))[p + 1]
  }
  rejections_all[n_z] = betaThatLinRegr_exp > quantile(betaThatLinRegr_replicates, 1 - alpha / 2)
}
cat("\n")

rejections_opt = array(NA, N_z)
for (n_z in 1 : N_z){
  if (n_z %% 1 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  idx_exp = sample(1 : rerand_res_norm$W_star_size, 1)
  w = W_base_obj$W_base_sorted[idx_exp, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr_exp = coef(lm(y ~ 0 + X + w))[p + 1]
  
  betaThatLinRegr_replicates = array(NA, R)
  for (r in 1 : R){
    w = W_base_obj$W_base_sorted[sample(setdiff(1 : rerand_res_norm$W_star_size, idx_exp), 1), ]
    betaThatLinRegr_replicates[r] = coef(lm(y ~ 0 + X + w))[p + 1]
  }
  rejections_opt[n_z] = betaThatLinRegr_exp > quantile(betaThatLinRegr_replicates, 1 - alpha / 2)
}
cat("\n")

rejections_top_41 = array(NA, N_z)
R = 40
for (n_z in 1 : N_z){
  if (n_z %% 1 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  idx_exp = sample(1 : R + 1, 1)
  w = W_base_obj$W_base_sorted[idx_exp, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr_exp = coef(lm(y ~ 0 + X + w))[p + 1]
  
  betaThatLinRegr_replicates = array(NA, R + 1)
  for (r in setdiff(1 : (R + 1), idx_exp)){
    w = W_base_obj$W_base_sorted[r, ]
    betaThatLinRegr_replicates[r] = coef(lm(y ~ 0 + X + w))[p + 1]
  }
  rejections_top_41[n_z] = betaThatLinRegr_exp > quantile(betaThatLinRegr_replicates, 1 - alpha / 2, na.rm = TRUE)
}
cat("\n")

rejections_top_101 = array(NA, N_z)
R = 100
for (n_z in 1 : N_z){
  if (n_z %% 1 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  idx_exp = sample(1 : R + 1, 1)
  w = W_base_obj$W_base_sorted[idx_exp, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr_exp = coef(lm(y ~ 0 + X + w))[p + 1]
  
  betaThatLinRegr_replicates = array(NA, R + 1)
  for (r in setdiff(1 : (R + 1), idx_exp)){
    w = W_base_obj$W_base_sorted[r, ]
    betaThatLinRegr_replicates[r] = coef(lm(y ~ 0 + X + w))[p + 1]
  }
  rejections_top_101[n_z] = betaThatLinRegr_exp > quantile(betaThatLinRegr_replicates, 1 - alpha / 2, na.rm = TRUE)
}
cat("\n")


rejections_top_201 = array(NA, N_z)
R = 200
for (n_z in 1 : N_z){
  if (n_z %% 1 == 0){cat(".")}
  z = rnorm(n, sd = sigma_z)
  idx_exp = sample(1 : R + 1, 1)
  w = W_base_obj$W_base_sorted[idx_exp, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr_exp = coef(lm(y ~ 0 + X + w))[p + 1]
  
  betaThatLinRegr_replicates = array(NA, R + 1)
  for (r in setdiff(1 : (R + 1), idx_exp)){
    w = W_base_obj$W_base_sorted[r, ]
    betaThatLinRegr_replicates[r] = coef(lm(y ~ 0 + X + w))[p + 1]
  }
  rejections_top_201[n_z] = betaThatLinRegr_exp > quantile(betaThatLinRegr_replicates, 1 - alpha / 2, na.rm = TRUE)
}
cat("\n")

mean(rejections_top_41)
mean(rejections_top_101)
mean(rejections_top_201)
mean(rejections_opt)
mean(rejections_all)


max_designs = 25000


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