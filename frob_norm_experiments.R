source("unique_perm.R")
pacman::p_load(ggplot2)

n = 16
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

frob_norm_sqs = array(NA, w_size)
max_entry_sq = array(NA, w_size)
entry_1_sq = array(NA, w_size)
entry_2_sq = array(NA, w_size)
entry_3_sq = array(NA, w_size)
entry_4_sq = array(NA, w_size)
entry_5_sq = array(NA, w_size)
entry_6_sq = array(NA, w_size)
entry_7_sq = array(NA, w_size)
entry_8_sq = array(NA, w_size)
min_entry_sq = array(NA, w_size)

for (last_i in w_size : 1){
  if (last_i %% 1000 == 0){
    cat(".")
  }
  idxs = res_iter_ord_obs_imb[1 : last_i, ]$i
  vecs = all_randomizations[idxs, , drop = FALSE]
  sigma_w = t(vecs) %*% vecs / nrow(vecs)
  entry_1_sq[last_i] = sum(sigma_w[1, 16]^2)
  entry_2_sq[last_i] = sum(sigma_w[2, 15]^2)
  entry_3_sq[last_i] = sum(sigma_w[3, 14]^2)
  entry_4_sq[last_i] = sum(sigma_w[4, 13]^2)
  entry_5_sq[last_i] = sum(sigma_w[5, 12]^2)
  entry_6_sq[last_i] = sum(sigma_w[6, 11]^2)
  entry_7_sq[last_i] = sum(sigma_w[7, 10]^2)
  entry_8_sq[last_i] = sum(sigma_w[8, 9]^2)
  sigma_w_no_diag = sigma_w
  sigma_w_no_diag = sigma_w_no_diag - diag(n)
  max_entry_sq[last_i] = max(sigma_w_no_diag^2)
  min_entry_sq[last_i] = min(sigma_w^2)
  
  frob_norm_sqs[last_i] = sum(sigma_w^2)
}
cat("\n")

ggplot(data.frame(size = 1 : w_size, frob_norm_sqs = frob_norm_sqs)) +
  geom_point(aes(x = size, y = frob_norm_sqs))

ggplot(data.frame(size = 1 : w_size, frob_norm_sqs = frob_norm_sqs)) +
  geom_point(aes(x = size, y = frob_norm_sqs)) + 
  xlim(100, 800) + ylim(20, 30)

ggplot(data.frame(imbalance_threshold = res_iter_ord_obs_imb$obs_imbalance, frob_norm_sqs = frob_norm_sqs)) +
  geom_point(aes(x = log(imbalance_threshold), y = frob_norm_sqs)) +
  ylim(25, 50) #min(frob_norm_sqs)

delta = 0.005
ggplot(data.frame(imbalance_threshold = res_iter_ord_obs_imb$obs_imbalance, entry_1_sq = entry_1_sq, entry_2_sq = entry_2_sq, entry_3_sq = entry_3_sq)) +
  aes(x = log(imbalance_threshold)) +
  xlab("entry squared value") + 
  geom_line(aes(y = frob_norm_sqs / n^2), col = "grey", linetype = "dashed", lwd = 2) + 
  geom_line(aes(y = entry_1_sq), col = "blue") + 
  geom_line(aes(y = entry_2_sq + delta), col = "red") + 
  geom_line(aes(y = entry_3_sq + 2 * delta), col = "green") + 
  geom_line(aes(y = entry_4_sq + 3 * delta), col = "purple") + 
  geom_line(aes(y = entry_5_sq + 4 * delta), col = "orange") + 
  geom_line(aes(y = entry_6_sq + 5 * delta), col = "darkgreen") +
  geom_line(aes(y = entry_7_sq + 6 * delta), col = "grey") + 
  geom_line(aes(y = entry_8_sq + 7 * delta), col = "yellow") +
  geom_line(aes(y = max_entry_sq), col = "black", linetype = "dashed", lwd = 2) +
  geom_line(aes(y = min_entry_sq), col = "white", linetype = "dashed", lwd = 2)



