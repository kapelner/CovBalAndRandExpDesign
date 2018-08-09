

n = 200

sigma_x = 1
mu_x = 0
sigma_z = 2
mu_z = 0

set.seed(0)
p = 1
X = matrix(rnorm(n * p, mu_x, sigma_x), ncol = p)
Sigma_X_inv = solve(var(X))

#the whole simulation is conditional on this one x

#results for linear case
beta_0 = 1
beta = rep(1, p)
w_0_forced_balance = c(rep(1, n / 2), rep(-1, n / 2))

Nsim = 100
Nrerand = 1e6

W = matrix(NA, nrow = n, ncol = Nrerand)
mahal_dists_scaled = array(NA, Nrerand)
for (w_i in seq(from = 2, to = Nrerand, by = 2)){
  #draw a random vector
  W[, w_i - 1] = sample(w_0_forced_balance) #rbinom(n, 1, 0.5) *2 - 1
  w_tr_X = t(X) %*% W[, w_i - 1]
  mahal_dists_scaled[w_i - 1] = t(w_tr_X) %*% Sigma_X_inv %*% w_tr_X
  #mirror it
  W[, w_i] = -W[, w_i - 1]
  w_tr_X = t(X) %*% W[, w_i]
  mahal_dists_scaled[w_i] = t(w_tr_X) %*% Sigma_X_inv %*% w_tr_X
}

#now we arrange by order and stop at Frobenius norm sqd 2n
order_of_mahal_dists_scaled = order(mahal_dists_scaled)

frob_normsq = array(NA, Nrerand)
for (w_i in seq(from = 2, to = Nrerand, by = 2)){
  W_0 = W[, order_of_mahal_dists_scaled[1 : w_i]]
  sigma_W_0 = var(t(W_0))
  frob_normsq[w_i] = sum(eigen(sigma_W_0)$values^2)
  if (w_i == 2 && frob_normsq[w_i] < 2 * n){
    stop("Never found the right Frobenius norm sq")
  }
  if (frob_normsq[w_i] < 2 * n){
    plot(1 : w_i, frob_normsq[1:w_i], ylim = c(0, n * 2 * 4)); abline(h = 2 * n, col = "green")
    cat("total vectors:", w_i, "\n")
    break
  }
}

# for (nsim in 1 : Nsim){
#   z = rnorm(n, mu_z, sigma_z)
#  
# }