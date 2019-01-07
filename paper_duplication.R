source("unique_perm.R")
pacman::p_load(xtable)

##Table 1 and Introduction
n = 8
x = rep(seq(1, n / 2), 1, each = 2) * 4
x = x - mean(x)
zval1 = 6
z = rep(c(zval1, -zval1), n / 2)
betaT = 1 #the PATE


#under complete randomization with forced balance
all_randomizations = uniqueperm2(rep(c(-1, 1), n / 2))
#under matching restriction
restricted_randomizations = all_randomizations[
  all_randomizations[, 1] + all_randomizations[, 2] == 0 &
  all_randomizations[, 3] + all_randomizations[, 4] == 0 &
  all_randomizations[, 5] + all_randomizations[, 6] == 0 &
  all_randomizations[, 7] + all_randomizations[, 8] == 0,
]

#Table 1
xtable(rbind(1 : n, x, z, ifelse(restricted_randomizations[10, ] == 1, "T", "C")))

#first two rows
res_iter = data.frame(matrix(NA, nrow = nrow(all_randomizations) + nrow(restricted_randomizations), ncol = 6))
colnames(res_iter) = c("procedure", "obs_imbalance", "unobs_imbalance", "naive_tx_estimate", "regression_estimate", "Rsq")

for (i in 1 : nrow(res_iter)){
  if (i <= nrow(all_randomizations)){
    indicT = all_randomizations[i, ]
    res_iter[i, 1] = "random"
  } else {
    indicT = restricted_randomizations[i - nrow(all_randomizations), ]
    res_iter[i, 1] = "restricted"
  }
  t_idx = indicT == 1
  c_idx = indicT == -1
  xT = x[t_idx]
  xC = x[c_idx]
  zT = z[t_idx]
  zC = z[c_idx]
  y = x + z + indicT * betaT
  yT = y[t_idx]
  yC = y[c_idx]
  res_iter[i, 2 : 6] = c(
    (mean(xT) - mean(xC))^2, 
    (mean(zT) - mean(zC))^2, 
    (mean(yT) - mean(yC)) / 2,
    coef(lm(y ~ 0 + x + indicT))[2],
    summary(lm(y ~ x))$r.squared
  )
}

res = data.frame(
  "avg_obs_imbalance" = c(
    mean(res_iter[res_iter$procedure == "random", "obs_imbalance"]),
    mean(res_iter[res_iter$procedure == "restricted", "obs_imbalance"])
  ),
  "avg_unobs_imbalance" = c(
    mean(res_iter[res_iter$procedure == "random", "unobs_imbalance"]),
    mean(res_iter[res_iter$procedure == "restricted", "unobs_imbalance"])
  ),
  "mse_naive" = c(
    mean((res_iter[res_iter$procedure == "random", "naive_tx_estimate"] - betaT)^2),
    mean((res_iter[res_iter$procedure == "restricted", "naive_tx_estimate"] - betaT)^2)
  ),
  "mse_regression" = c(
    mean((res_iter[res_iter$procedure == "random", "regression_estimate"] - betaT)^2),
    mean((res_iter[res_iter$procedure == "restricted", "regression_estimate"] - betaT)^2)
  )
)
rownames(res) = c("random", "restricted")
#table 2
xtable(res)
#now let's see check the $R^2$
mean(res_iter$Rsq)

