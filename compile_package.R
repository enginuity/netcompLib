## Contains code for testing and compiling package.

library(testthat)

## Run tests

library(metacode)
library(stringr)

library(roxygen2)

update_fx_documentation(FD = FilesDescription(dirlist = "netcomplib/"), fill_emptyparam = FALSE)
roxygenise("netcomplib/", clean = TRUE)



## ONLY IF PASSES TEST -> INSTALL PACKAGE

system("R CMD INSTALL netcomplib")

