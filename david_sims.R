source("unique_perm.R")
if (!require("pacman")){install.packages("pacman")}
pacman::p_load(lmenssp, ggplot2)

n = 20
p = 1

sigma_x = 1
mu_x = 0
sigma_z = 1.5
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

#calculate all possible balances
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

#then grab the optimal vector
w_star = all_randomizations[res_iter_ord_obs_imb$i[1], ]
t_idx = w_star == 1
c_idx = w_star == -1

Nsim = 2000
beta_hat_lin_reg_adj = array(NA, Nsim)
beta_hat_classic = array(NA, Nsim)
reject_lin_reg = array(NA, Nsim)
reject_classic = array(NA, Nsim)
alpha = 0.05

#since we are testing the null distribution
beta_T = 0
for (nsim in 1 : Nsim){
  #pull a random z out of a hat
  z = rnorm(n, 0, sigma_z)
  #generate the outcome
  y = beta_0 + X %*% bbeta + z + w_star * beta_T
  #pull the estimate of beta_T from the linear regression
  lm_mod = lm(y ~ X + w_star)
  beta_hat_lin_reg_adj[nsim] =  coef(lm_mod)[3] / 2
  yT = y[t_idx]
  yC = y[c_idx]
  beta_hat_classic[nsim] = (mean(yT) - mean(yC)) / 2
  reject_lin_reg[nsim] = coef(summary(lm_mod))[3, 4] < alpha
  reject_classic = t.test(yT, yC)$p.value < alpha
}

df = n - 3
qqplot.t(beta_hat_lin_reg_adj, df)
theoretically_predicted_null_distr = rt(Nsim, df) 
ks.test(beta_hat_lin_reg_adj, theoretically_predicted_null_distr)

ggplot(data.frame(theoretically_predicted_null_distr = theoretically_predicted_null_distr, beta_hat_lin_reg_adj = beta_hat_lin_reg_adj)) + 
  geom_density(aes(theoretically_predicted_null_distr), alpha = 0.3, fill = "red") + 
  geom_density(aes(beta_hat_lin_reg_adj), alpha = 0.3, fill = "blue")

mean(reject_lin_reg)
mean(reject_classic)

