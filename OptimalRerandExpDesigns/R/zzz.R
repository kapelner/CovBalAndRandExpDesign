
.onAttach = function(libname, pkgname){
  packageStartupMessage(
	  paste("Welcome to OptimalRerandExpDesigns v", 
		  utils::packageVersion("OptimalRerandExpDesigns"), 
		  ".", sep = "")
	)
}