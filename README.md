# CovBalAndRandExpDesign

An R package for computing near optimal experimental designs for a two-arm experiment with covariates hedging for the covariates you don't observe.

The repository contains a testscript as well for working through some examples.

This is joint work with Abba M. Krieger of the Wharton School of the University of Pennsylvania, Uri Shalit and David Azriel of The Technion.

To load the package, make sure `rJava' is installed and properly configured! There are some issues with the latest Java 10. Make sure you read about this online. For now, Java 7 or 8 is recommended. Then:

	options(java.parameters = "-Xmx4000m")
	library(CovBalAndRandExpDesign)
	
And if you want to use the optimization feature via `Gurobi', install Gurobi first and then run something like:

	.jaddClassPath("/gurobi801/win64/lib/gurobi.jar")
		
Note that you have to register your machine and run `grbgetkey` with the proper license number and ensure the `.lic` file is installed. Expect `ClassNotFoundException`s if this is not setup properly!
