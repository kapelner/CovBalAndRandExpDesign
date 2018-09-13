####calculate some Frob norm sq:
r = 10000
rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4)
res = resultsGreedySearch(rd, max_vectors = r)

all_vecs_opt = res$ending_indicTs

r = 10000
rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4)
res = resultsGreedySearch(rd, max_vectors = r)

all_vecs_mid = res$ending_indicTs
#make one section random
prob_randomize = 0.01
length_randomize = 2
for (m in 1 : r){
  indicT = all_vecs_mid[, m]
  for (i_w in length_randomize : n){
    if (rbinom(1, 1, prob_randomize) == 1 & sum(indicT[(i_w - length_randomize) : i_w]) == length_randomize / 2){
      indicT[(i_w - length_randomize) : i_w] = sample(indicT[(i_w - length_randomize) : i_w])
    }    
  }
  all_vecs_mid[, m] = indicT
}


num_match = 10000
all_vecs_match = matrix(NA, nrow = n, ncol = num_match)
for (m in 1 : num_match){
  indicT = array(NA, n)
  for (i_w in seq(2, n, by = 2)){
    indicT[c(i_w - 1, i_w)] = sample(c(0, 1))
  }
  all_vecs_match[, m] = indicT
}
all_vecs_rand = t(complete_randomization_with_forced_balanced(n, num_match))
all_vecs_rand[all_vecs_rand == -1] = 0
all_vecs = cbind(all_vecs_opt, all_vecs_opt2, all_vecs_match, all_vecs_rand)
all_vecs = cbind(all_vecs_mid)
all_vecs = cbind(all_vecs_mid, all_vecs_opt, all_vecs_match, all_vecs_rand)

all_mahal_obj = apply(all_vecs, 2, function(w){compute_objective_val(X, w, objective = "mahal_dist")})
all_vecs = all_vecs[, order(all_mahal_obj)]
all_mahal_obj_sort_log10 = sort(log10(all_mahal_obj))
hist(all_mahal_obj_sort_log10, br = 1000)


##now create a representative sample
num_per_order_of_mag = 200
round_all_mahal_obj_sort_log10 = round(all_mahal_obj_sort_log10)
all_vecs_rep = matrix(NA, nrow = n, ncol = 0)
for (level_log in min(round_all_mahal_obj_sort_log10) : max(round_all_mahal_obj_sort_log10)){
  idx = which(round_all_mahal_obj_sort_log10 == level_log)
  if (length(idx) > num_per_order_of_mag){
    idx = sample(idx, num_per_order_of_mag)
  }
  all_vecs_rep = cbind(all_vecs_rep, all_vecs[, idx])
}

all_mahal_obj_rep = apply(all_vecs_rep, 2, function(w){compute_objective_val(X, w, objective = "mahal_dist")})
all_vecs_rep = all_vecs_rep[, order(all_mahal_obj_rep)]
all_mahal_obj_rep = sort(all_mahal_obj_rep)
hist(log10(all_mahal_obj_rep), br = 1000)

quantile(all_mahal_obj_rep, seq(0, 1, by = 0.1))
num_rep_vecs = ncol(all_vecs_rep)
pct_num_rep_vecs = (1 : num_rep_vecs) / num_rep_vecs

min_criterions_sigma_z_i = list()
sigma_zs = c(0.1, 0.3, 1, 3, 10)
for (sigma_z in sigma_zs){
  if (sigma_z == sigma_zs[1]){
    R_terms = array(NA, num_rep_vecs)
    B1_terms = array(NA, num_rep_vecs)
    B2_terms = array(NA, num_rep_vecs)
  }
  criterion_values = array(NA, num_rep_vecs)
  for (i in seq(1, num_rep_vecs, by = 10)){
    if (sigma_z == sigma_zs[1]){
      sigma_w_optimals = matrix(0, n, n)
      for (j in 1 : i){
        w = all_vecs_rep[, j]
        w[w == 0] = -1
        sigma_w_optimals = sigma_w_optimals + 1 / i * w %*% t(w)
      }
      R_terms[i] = sum(eigen(sigma_w_optimals)$values^2)
      B1_terms[i] = t(X) %*% sigma_w_optimals %*% X
      B2_terms[i] = t(X) %*% sigma_w_optimals %*% sigma_w_optimals %*% X
    }
    criterion_values[i] = 1 / n^2 * (B1_terms[i] + n * sigma_z^2 + 2 * sqrt(2 * sigma_z^4 *  R_terms[i] + 4 * sigma_z^2 * B2_terms[i]))
  }
  B1_terms[seq(1, num_rep_vecs, by = 10)]
  B2_terms[seq(1, num_rep_vecs, by = 10)]
  R_terms[seq(1, num_rep_vecs, by = 10)]
  criterion_values[seq(1, num_rep_vecs, by = 10)]
  
  
  sigma_z = as.character(sigma_z)
  min(criterion_values, na.rm = TRUE)
  min_criterions_sigma_z_i[[sigma_z]] = which(criterion_values == min(criterion_values, na.rm = TRUE)) / num_rep_vecs
  
  
  ggplot(data.frame(pct_vecs = pct_num_rep_vecs, B1 = B1_terms)) + 
    geom_point(aes(pct_vecs, B1), col = "red") + 
    geom_vline(xintercept = min_criterions_sigma_z_i[[sigma_z]], col = "black")
  ggplot(data.frame(pct_vecs = pct_num_rep_vecs, B2 = B2_terms)) + 
    geom_point(aes(pct_vecs, B2), col = "red") + 
    geom_vline(xintercept = min_criterions_sigma_z_i[[sigma_z]], col = "black")
  ggplot(data.frame(pct_vecs = pct_num_rep_vecs, R = R_terms)) + 
    ylim(0, 1000) +
    geom_point(aes(pct_vecs, R), col = "blue") + 
    geom_vline(xintercept = min_criterions_sigma_z_i[[sigma_z]], col = "black")
  ggplot(data.frame(pct_vecs = pct_num_rep_vecs, tail_criterion = criterion_values)) + 
    # ylim(0, 0.07) +
    geom_point(aes(pct_vecs, tail_criterion), col = "purple") + 
    geom_vline(xintercept = min_criterions_sigma_z_i[[sigma_z]], col = "black")
  
}









sigma_z = 10






It seems with "reasonable" levels of sigma_z (relative to sigma_x = 1), the answer is always around 70%. It is important to restate what this means again. 70% of the vectors means from optimal all the way up to about 10^-6. This means designs consisting of optimal all the way up to the best of matching considering uniformly spaced vectors on the log10 scale of imbalance (which admittedly is a strange design procedure). The most restricted it gets is when sigma_z = 0.001 (basically R^2 = 100%) and here we see that we don't demand anything better than 43% of the vectors i.e. up to 10^-11 or better in imbalance. This imbalance better than matching but much worse than greedy pair switching. The least restricted is when sigma_z = 10 (basically R^2 = 0%) where we want 94% of the vectors or about 10^-3 in imbalance or better than the worst matching / best CRFB's.




I find this to be an interesting observation.



Adam