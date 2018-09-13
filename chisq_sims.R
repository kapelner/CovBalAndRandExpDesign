source("unique_perm.R")
source("_one_bound_sim.R")
pacman::p_load(ggplot2, rmutil)



n = 16
sigma_x = 1
sigma_z = 1
mu_z = 0
mu_x = 0
KEEP_SIGMA = TRUE

set.seed(1984)
x = rnorm(n, mu_x, sigma_x)

#let's take a sigma_w in the middle with around lambda_max = 2
sigma_w = sigma_w_obs_imb[[250]]
max(eigen(sigma_w)$values)
sum(eigen(sigma_w)$values^2)

#what does the quadratic form look like over z?
Nsim = 10000
vals = array(NA, Nsim)
for (nsim in 1 : Nsim){
  # z = rnorm(n, mu_z, sigma_z)
  z = runif(n, -5, 5)
  # z = rexp(n, 5) - 1/5
  # z = c(rnorm(n / 2, 0, 1), rnorm(n / 2, 10, 1)) - 5
  # z = rlaplace(n)
  vals[nsim] = t(x + z) %*% sigma_w %*% (x + z)
}

ggplot(data.frame(quad_form_val = vals)) + geom_histogram(aes(x = quad_form_val), bins = 100)
# ggplot(data.frame(z = z)) + geom_histogram(aes(x = z), bins = 100)


