#' This method creates an object of type optimal_experimental_design and will immediately initiate
#' a search through $1_{T}$ space.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param objective			The objective function to use when greedily searching design space. This is a string
#' 							"\code{abs_sum_diff}" (default) or "\code{mahal_dist}."
#' @param wait				Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 							deafult is \code{FALSE}.
#' @param start				Should we start searching immediately (default is \code{TRUE}).
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @return					An object of type \code{optimal_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initOptimalExperimentalDesignObject = function(X,
		objective = "abs_sum_diff", 
		wait = FALSE, 
		start = TRUE,
		num_cores = 1){
	#get dimensions immediately
	n = nrow(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	p = ncol(X)
	
	if (objective == "abs_sum_diff"){
		#standardize it -- much faster here
		Xstd = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
	}
	if (objective == "mahal_dist"){
		if (p < n){
			SinvX = solve(var(X))
		}
	}
	
	#we are about to construct a OptimalExperimentalDesign java object. First, let R garbage collect
	#to clean up previous objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("OptimalExperimentalDesign.OptimalExperimentalDesign")
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setNandP", as.integer(n), as.integer(p))
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}	
	
	#feed in the data
	for (i in 1 : n){	
		if (objective == "abs_sum_diff"){
			.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), Xstd[i, , drop = FALSE]) #java indexes from 0...n-1
		} else {
			.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
		}
	}
	
	#feed in the inverse var-cov matrix
	if (objective == "mahal_dist"){
		if (p < n){
			for (j in 1 : p){
				.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
			}
		}
	}
		
	#now return information as an object (just a list)
	optimal_experimental_design_search = list()
	optimal_experimental_design_search$start = start
	optimal_experimental_design_search$wait = wait
	optimal_experimental_design_search$X = X
	optimal_experimental_design_search$n = n
	optimal_experimental_design_search$p = p
	optimal_experimental_design_search$objective = objective
	optimal_experimental_design_search$java_obj = java_obj
	class(optimal_experimental_design_search) = "optimal_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(optimal_experimental_design_search)
	}
	#return the final object
	optimal_experimental_design_search
}

#' Prints a summary of a \code{optimal_experimental_design_search} object
#' 
#' @param x			The \code{optimal_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print optimal_experimental_design_search
#' @export
print.optimal_experimental_design_search = function(x, ...){
	progress = .jcall(x$java_obj, "D", "progress")
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the OptimalExperimentalDesign. Did you run \"startOptimalSearch?\"\n")
	} else if (progress < 1){
		cat("The search is", round(progress * 100, 2), "% complete.\n")
	} else {
		cat("The search is complete.\n")
	}
}

#' Prints a summary of a \code{optimal_experimental_design_search} object
#' 
#' @param object		The \code{optimal_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary optimal_experimental_design_search
#' @export
summary.optimal_experimental_design_search = function(object, ...){
	print(object, ...)
}

#' Returns the results (thus far) of the optimal design search
#' 
#' @param obj 			The \code{optimal_experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
resultsOptimalSearch = function(obj){
	opt_obj_val = .jcall(obj$java_obj, "D", "getOptObjectiveVal")
	indicT = .jcall(obj$java_obj, "[I", "getOptIndicT", .jevalArray)
	all_obj_vals = .jcall(obj$java_obj, "[D", "getAllObjectiveVals", .jevalArray)
	list(
		opt_obj_val = opt_obj_val,
		indicT = indicT,
		all_obj_vals = all_obj_vals
	)
}
