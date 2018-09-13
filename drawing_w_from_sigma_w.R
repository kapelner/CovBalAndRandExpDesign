pacman::p_load(mvtnorm)

###first try a basic CRFB example
n = 4
D = (1 + 1 / (n-1)) * diag(n) - 1 / (n - 1)
# small perturbations???
# for (i in 1 : n){
#   for (j in 1 : i){
#     D[i, j] = D[i, j] + rnorm(1, 0, 0.001)
#     D[j, i] = D[i, j]
#   }
# }

###try an optimal example
w_star = c(rep(1, n / 2), rep(-1, n / 2))
D = w_star %*% t(w_star)  

###try matching
D = diag(n) 
for (i in seq(2, n, by = 2)){
  D[i, i - 1] = -1
  D[i - 1, i] = -1
}

###sample many CRFB vectors
Nvec = 10000
D = matrix(0, n, n)
for (nvec in 1 : Nvec){
  w = sample(c(rep(1, n / 2), rep(-1, n / 2)))
  D = D + w %*% t(w)
}
D = 1 / Nvec * D
D

###sample many CR vectors
Nvec = 500
D = matrix(0, n, n)
for (nvec in 1 : Nvec){
  w = (rbinom(n, 1, 0.5) - 0.5) * 2
  D = D + w %*% t(w)
}
D = 1 / Nvec * D
D

G = diag(n + 1)
G[1 : n, 1 : n] = D
H = sin(pi * G / 2)


Nsim = 1000

sigma_w_sim = matrix(0, nrow = n, ncol = n)
for (nsim in 1 : Nsim){
  a = rmvnorm(n = 1, mean = rep(0, n + 1), sigma = H)
  b = sign(a)
  w = (b * b[n + 1])[1 : n]
  sigma_w_sim = sigma_w_sim + w %*% t(w)
}
sigma_w_sim = 1 / Nsim * sigma_w_sim
D
round(sigma_w_sim, 3)
sum(abs(sigma_w_sim - D)) / n^2
