source("unique_perm.R")
pacman::p_load(xtable)

##Table 1 and Introduction
n = 6
x = rep(seq(1, n / 2), 1, each = 2)
zval1 = 1.5
zval2 = -1.5
z = rep(c(zval1, zval2), n / 2)
betaT = 1 #the PATE

#under complete randomization with forced balance
all_randomizations = uniqueperm2(rep(c(-1, 1), n / 2))
#under matching restriction
restricted_randomizations = all_randomizations[
  all_randomizations[, 1] + all_randomizations[, 2] == 0 &
  all_randomizations[, 3] + all_randomizations[, 4] == 0 &
  all_randomizations[, 5] + all_randomizations[, 6] == 0,
]

#first two rows
res_iter = data.frame(matrix(NA, nrow = nrow(all_randomizations) + nrow(restricted_randomizations), ncol = 5))
colnames(res_iter) = c("procedure", "obs_imbalance", "unobs_imbalance", "tx_estimate", "Rsq")

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
  res_iter[i, 2 : 5] = c(
    (mean(xT) - mean(xC))^2, 
    (mean(zT) - mean(zC))^2, 
    (mean(yT) - mean(yC)) / 2,
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
  "mse_tx_est" = c(
    mean((res_iter[res_iter$procedure == "random", "tx_estimate"] - betaT)^2),
    mean((res_iter[res_iter$procedure == "restricted", "tx_estimate"] - betaT)^2)
  )
)
rownames(res) = c("random", "restricted")
xtable(res)
#now let's see the $R^2$
mean(res_iter$Rsq)







#last two rows
# Nrep = 200
# res_iter = data.frame(matrix(NA, nrow = Nrep * (nrow(all_randomizations) + nrow(restricted_randomizations)), ncol = 4))
# colnames(res_iter) = c("procedure", "obs_imbalance", "tx_estimate", "Rsq")
# 
# for (i in 1 : nrow(res_iter)){
#   if (i <= Nrep * nrow(all_randomizations)){
#     indicT = all_randomizations[ceiling(i / Nrep), ]
#     res_iter[i, 1] = "random"
#   } else {
#     indicT = restricted_randomizations[ceiling(i / Nrep) - nrow(all_randomizations), ]
#     res_iter[i, 1] = "restricted"
#   }
#   t_idx = indicT == 1
#   c_idx = indicT == -1
#   xT = x[t_idx]
#   xC = x[c_idx]
#   y = x + sample(c(zval1, -zval1), n, replace = TRUE) + indicT * betaT
#   yT = y[t_idx]
#   yC = y[c_idx]
#   
#   res_iter[i, 2 : 4] = c(
#     (mean(xT) - mean(xC))^2,
#     mean(yT) - mean(yC),
#     summary(lm(y ~ x))$r.squared
#   )
# }
# 
# res = data.frame(
#   "avg_obs_imbalance" = c(
#     mean(res_iter[res_iter$procedure == "random", "obs_imbalance"]),
#     mean(res_iter[res_iter$procedure == "restricted", "obs_imbalance"])
#   ),
#   "mse_tx_est" = c(
#     mean((res_iter[res_iter$procedure == "random", "tx_estimate"] - betaT)^2),
#     mean((res_iter[res_iter$procedure == "restricted", "tx_estimate"] - betaT)^2)
#   )
# )
# rownames(res) = c("random", "restricted")
# xtable(res)
# #now let's see the $R^2$
# mean(res_iter[res_iter$procedure == "random", "Rsq"])



