

options(java.parameters = "-Xmx4000m")
pacman::p_load(ggplot2, GreedyExperimentalDesign)

n = 100
r = 1200000

set.seed(1984)
x = rnorm(n)


cr_w_s = complete_randomization(n, r)

bals = cr_w_s %*% x
c(min(bals), max(bals))

num_bins_each_dir_from_zero = 50
bin_width = sd(bals) / num_bins_each_dir_from_zero * 1.2

vecs_in_bins = list()
num_in_bins = list()
for (bin_num in -num_bins_each_dir_from_zero : (num_bins_each_dir_from_zero - 1)){
  key = paste(as.character(bin_num), "...", as.character(bin_num + 1))
  
  a = bin_num * bin_width
  b = (bin_num + 1) * bin_width
  w_i_s = which(a <= bals & bals < b)
  cr_w_s_bucket = cr_w_s[w_i_s, ]
  num_in_bin = nrow(cr_w_s_bucket)
  num_in_bins[[key]] = num_in_bin
  
  vecs_in_bins[[key]] = 
    t(cr_w_s_bucket) %*% cr_w_s_bucket / num_in_bin
}

res = data.frame(
  i = -num_bins_each_dir_from_zero : (num_bins_each_dir_from_zero - 1) + 0.5,
  num_in_bin = unlist(num_in_bins),
  frob_norm_sq = unlist(lapply(vecs_in_bins, function(W){sum(diag(W %*% W))}))
)
res

quartic_model = lm(frob_norm_sq ~ poly(i, 4, raw = TRUE), res)

quartic_model_fun = function(x){
  coef(quartic_model)[1] +
    coef(quartic_model)[2] * x +
    coef(quartic_model)[3] * x^2 +
    coef(quartic_model)[4] * x^3 +
    coef(quartic_model)[5] * x^4
}

ggplot(res, aes(x = i, y = frob_norm_sq)) + geom_point() + 
  stat_function(fun = quartic_model_fun) 
#+ 
 # ylim(n, max(res$frob_norm_sq))




# ggplot(data.frame(bals = bals)) + geom_histogram(aes(x = bals), binwidth = 0.5)



