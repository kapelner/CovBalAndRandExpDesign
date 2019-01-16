library(OptimalRerandExpDesigns)

set.seed(1984)
n = 100
p = 10
X = matrix(rnorm(n * p), nrow = n, ncol = p)
X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})

W_base_obj = generate_W_base_and_sort(X)
W_base_obj
#Fig 1
plot(W_base_obj, title = "", xlab = "Mahalanobis Distance")


rerand_res_norm = optimal_rerandomization_normality_assumed(W_base_obj)
rerand_res_norm
plot(rerand_res_norm)

#Fig 2
set.seed(1984)
bbeta = rnorm(p)
betaT = 1
N_w = 1000
N_z = 1000
sigma_z = 3

mse_z_all = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = W_base_obj$W_base_sorted[sample(1 : 25000, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_all[n_z] = mean(mse_w)
}


mse_z_top_5000 = array(NA, N_z)
for (n_z in 1 : N_z){
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


mse_z_worst_5000 = array(NA, N_z)
for (n_z in 1 : N_z){
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


mse_z_opt = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  mse_w = array(NA, N_w)
  for (n_w in 1 : N_w){
    w = rerand_res_norm$W_star[sample(1 : rerand_res_norm$w_star_size, 1), ]
    y = X %*% bbeta + z + betaT * w
    betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
    mse_w[n_w] = (betaThatLinRegr - betaT)^2
  }
  mse_z_opt[n_z] = mean(mse_w)
}


mse_z_determ = array(NA, N_z)
for (n_z in 1 : N_z){
  z = rnorm(n, sd = sigma_z)
  w = W_base_obj$W_base_sorted[1, ]
  y = X %*% bbeta + z + betaT * w
  betaThatLinRegr = coef(lm(y ~ 0 + X + w))[p + 1]
  mse_z_determ[n_z] = (betaThatLinRegr - betaT)^2
}


rerand_res_ex = optimal_rerandomization_exact(W_base_sorted, z_sim_fun = function(){rnorm(n)})
rerand_res_ex
plot(rerand_res_ex)

rerand_res_approx = optimal_rerandomization_tail_approx(W_base_sorted)
rerand_res_approx
plot(rerand_res_approx)


ggplot(data.frame(Qs = rerand_res_approx$all_data_from_run$Q_primes[1 : 100])) + aes(x = Qs) + geom_density()
ggplot(data.frame(Qs = rerand_res_approx$all_data_from_run$Q_primes[1 : 200])) + aes(x = Qs) + geom_density()
ggplot(data.frame(Qs = rerand_res_approx$all_data_from_run$Q_primes[1 : 300])) + aes(x = Qs) + geom_density()
ggplot(data.frame(Qs = rerand_res_approx$all_data_from_run$Q_primes[1 : 400])) + aes(x = Qs) + geom_density()
ggplot(data.frame(Qs = rerand_res_approx$all_data_from_run$Q_primes[1 : 500])) + aes(x = Qs) + geom_density()
ggplot(data.frame(Qs = rerand_res_approx$all_data_from_run$Q_primes[1 : 600])) + aes(x = Qs) + geom_density()


# ggplot(x$all_data_from_run) +
#   ggtitle(title, subtitle = subtitle) +
#   xlab(xlab) +
#   ylab(ylab) +
#   geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
#   geom_vline(xintercept = x$a_star, col = "green") +
#   scale_x_log10(   
#     # breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
#     # labels = scales::trans_format("log10", scales::scientific_format(digits = 2))
#     ) +
#   scale_y_log10()
#   # coord_trans(y = "log10")










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