
pacman::p_load(ggplot2)

n = 100

sigma_x = 1
mu_x = 0
sigma_z = 0.0001
mu_z = 0

set.seed(0)
p = 1
X = matrix(rnorm(n * p, mu_x, sigma_x), ncol = p)
Sigma_X_inv = solve(var(X))

#the whole simulation is conditional on this one x

#results for linear case
beta_0 = 1
beta_X = rep(1, p)
w_0_forced_balance = c(rep(1, n / 2), rep(-1, n / 2))

Nsim = 100
Nrerand = 1e4

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


obj_function = array(NA, Nrerand)
frob_norm_sqs_function = array(NA, Nrerand)
# APPROX = 2000

#trace out our objective function
f = X %*% beta_X
c_const = 1

for (w_i in seq(from = Nrerand, to = 2, by = -10)){
  # if (w_i > APPROX){
  #   W_0 = W[, sample(order_of_mahal_dists_scaled[1 : w_i], APPROX)]
  # } else {
    W_0 = W[, order_of_mahal_dists_scaled[1 : w_i]]
  # }
  
  sigma_W_0 = var(t(W_0))
  frob_norm_sqs_function[w_i] = sum(eigen(sigma_W_0)$values^2)
  obj_function[w_i] = 
    1 / n^2 * t(f) %*% sigma_W_0 %*% f +
    1 / n * sigma_z^2 + 
    c_const * sqrt(
      0 +
      2 * sigma_z^2^2 / n^4 * frob_norm_sqs_function[w_i] +
      4 * sigma_z^2^2 / n^4 * t(f) %*% sigma_W_0^2 %*% f
    )
  if (w_i %% 1000 == 0){
    cat("w_i", w_i, "\n")
  }
  # if (frob_norm_sqs_function[w_i] > 2 * n){
  #   break
  # }
}

  ggplot(data.frame(frob_norm_sqs_function = frob_norm_sqs_function, obj_function = obj_function)) +
    geom_line(aes(x = frob_norm_sqs_function, y = obj_function)) #+ xlim(100, 300) #+ ylim(0.01,0.02)


  which.min(obj_function)
  frob_norm_sqs_function[10]
# hist(frob_norm_sqs_function, br = 1000,xlim = c(0, 1000))
