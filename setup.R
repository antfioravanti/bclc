script_root <- dirname(rstudioapi::getActiveDocumentContext()$path)
proj_root   <- normalizePath(file.path(script_root, ".."), winslash = "/")
setwd(script_root)

setwd("bclc")
Rcpp::compileAttributes()

# 1. Set your working directory TO bclc/, not the parent
setwd("bclc")   # or open bclc/ as its own RStudio project

# 2. Add only what's missing, without re-creating the whole scaffold
usethis::use_description(fields = list(
  Title = "Bias-Corrected Lattice Covariance Estimation",
  Version = "0.1.0"
))
usethis::use_mit_license("Antonio Fioravanti")
usethis::use_rcpp_armadillo()   # adds Rcpp/RcppArmadillo to DESCRIPTION
Rcpp::compileAttributes()   # generates R/RcppExports.R and src/RcppExports.cpp
devtools::document()        # generates NAMESPACE and man/ from roxygen comments
devtools::check()           # validates everything
