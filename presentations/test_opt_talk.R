
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

NUM_CORES = 3
MINN = 4
MAXN = 26
AVG_GREEDY = 100
reps = 100
ns = seq(from = MINN, to = MAXN, by = 2)
ps = c(1) #c(1, 2, 5, 10, 20)
rs = c(1) #c(1, 5, 10, 100, 1000, 10000)
opt_obj_vals = matrix(NA, length(ns), reps)
greedy_obj_vals = matrix(NA, length(ns), reps)

for (rep in 1 : reps){
	for (n_i in 1 : length(ns)){
		for (i in 1 : length(ps)){
			X = generate_stdzied_design_matrix(n = ns[n_i], p = ps[i])
			
			for (r in 1 : length(rs)){
				ged = initGreedyExperimentalDesignObject(X, max_designs = AVG_GREEDY, num_cores = NUM_CORES, wait = TRUE)
				greedy_obj_vals[n_i, rep] = mean(resultsGreedySearch(ged, max_vectors = 0)$obj_vals)
			}
			
			oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "abs_sum_diff", wait = TRUE)
			opt_obj_vals[n_i, rep] = resultsOptimalSearch(oed)$obj_val
		}	
	}
}


log_greedy_obj_vals = log(rowMeans(greedy_obj_vals)) / log(10)
log_opt_obj_vals = log(rowMeans(opt_obj_vals)) / log(10)


plot(ns, log_greedy_obj_vals,
		ylim = c(min(log_greedy_obj_vals, log_opt_obj_vals), max(log_greedy_obj_vals, log_opt_obj_vals)),
		ylab = "log10 obj function", xlab = "n", main = paste("greedy switch vs optimal"), type = "o", col = "blue")
points(ns, log_opt_obj_vals - 0.02, type = "o", col = "green")
