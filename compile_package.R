##@S Contains code to compile the netcompLib package (to be run from R)


# Update Documentation (if needed/possible) -------------------------------

## This section only needs to be run if the source code is not necessarily up to date (ie new functions were written, or if parameters were added). It will also only be run if the package actually exists. 
if (require(codeProcessing)) {
  ## This is a package I've written to speed up my own coding efficiency to use when rewriting code / package writing
  
  ## Update documentation -- this looks for new parameters and creates documentation lines for it
  update_fx_documentation(FD = FilesDescription(dirlist = "netcompLib/R/"), test_run = FALSE, regexp_nodocu = "[.]")
}

# Compile Package ---------------------------------------------------------

if (require(roxygen2)) {
  ## Generate the documentation -- THIS MUST be run before building packages (since the documentation files are not version-controlled, as the version-controlled version is in the raw source code)
  roxygenise("netcompLib/", clean = TRUE)
  
  ## Install the package
  system("R CMD INSTALL netcompLib")
} else {
  stop("Package 'roxygen2' is not installed: the netcompLib package is not compilable.")
}

if (FALSE) {
  ## Code to load the library...
  library(netcompLib)
}



## TODO: Create function in codeProcessing that reorders generics & makes templates and such. 
