
.onAttach = function(libname, pkgname){
  packageStartupMessage(
	  paste("Welcome to CovBalAndRandExpDesign v", 
		  utils::packageVersion("CovBalAndRandExpDesign"), 
		  ".\n", sep = "")
	)
}