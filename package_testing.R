library(CovBalAndRandExpDesign)



















max_designs = 25000
n = 100
p = 10
X = matrix(rnorm(n * p), nrow = n, ncol = p)

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