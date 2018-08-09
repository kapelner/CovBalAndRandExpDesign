
options(java.parameters = "-Xmx4000m")
library(GreedyExperimentalDesign)

NUM_CORES = 3

ns = c(10, 100, 1000)
ps = c(1, 2, 5, 10)
x_reps = 5
g_reps = 6

greedy_obj_vals = array(NA, c(length(ns), length(ps), x_reps))
num_switches = array(NA, c(length(ns), length(ps), x_reps))

for (rep in 1 : x_reps){
	for (n_i in 1 : length(ns)){
		for (p_i in 1 : length(ps)){
			cat("n", ns[n_i], "p", ps[p_i], "rep", rep, "\n")
			X = generate_stdzied_design_matrix(n = ns[n_i], p = ps[p_i])
			
			ged = initGreedyExperimentalDesignObject(X, max_designs = g_reps, num_cores = NUM_CORES, wait = TRUE)
			res = resultsGreedySearch(ged, max_vectors = 2)
			greedy_obj_vals[n_i, p_i, rep] = mean(res$obj_vals)
			num_switches[n_i, p_i, rep] = mean(res$num_iters)
			
			save.image("talk_data.RData")
		}	
	}
}

num_switches_avg = rowMeans(num_switches, dim = 2, na.rm = TRUE)

greedy_obj_vals_avg_log = log(rowMeans(greedy_obj_vals, dim = 2, na.rm = TRUE))

plot(log(ns), greedy_obj_vals_avg_log[, 1], col = "blue", type = "o",
		ylim = c(min(greedy_obj_vals_avg_log), max(greedy_obj_vals_avg_log)),
		xlab = "log(n)", ylab = "log(O)", main = "log(O) by log(n) for p = 1, 2, 5, 10")
points(log(ns), greedy_obj_vals_avg_log[, 2], col = "red", type = "o")
points(log(ns), greedy_obj_vals_avg_log[, 3], col = "green", type = "o")
points(log(ns), greedy_obj_vals_avg_log[, 4], col = "purple", type = "o")

mod_p_1 = lm(greedy_obj_vals_avg_log[, 1] ~ log(ns))
summary(mod_p_1)
mod_p_2 = lm(greedy_obj_vals_avg_log[, 2] ~ log(ns))
summary(mod_p_2)
mod_p_5 = lm(greedy_obj_vals_avg_log[, 3] ~ log(ns))
summary(mod_p_5)
mod_p_10 = lm(greedy_obj_vals_avg_log[, 4] ~ log(ns))
summary(mod_p_10)

manip_ns = ns

plot(manip_ns, num_switches_avg[, 1], col = "blue", type = "o",
		ylim = c(min(num_switches_avg), max(num_switches_avg)),
		xlab = "n", ylab = "# switches", main = "#switches by n for p = 1, 2, 5, 10")
points(manip_ns, num_switches_avg[, 2], col = "red", type = "o")
points(manip_ns, num_switches_avg[, 3], col = "green", type = "o")
points(manip_ns, num_switches_avg[, 4], col = "purple", type = "o")

mod_p_1 = lm(num_switches_avg[, 1] ~ manip_ns)
summary(mod_p_1)
mod_p_2 = lm(num_switches_avg[, 2] ~ manip_ns)
summary(mod_p_2)
mod_p_5 = lm(num_switches_avg[, 3] ~ manip_ns)
summary(mod_p_5)
mod_p_10 = lm(num_switches_avg[, 4] ~ manip_ns)
summary(mod_p_10)

