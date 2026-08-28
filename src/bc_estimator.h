#pragma once
#include <RcppArmadillo.h>
using namespace Rcpp;

// Bias-corrected estimator, defined in bc_estimator.cpp.
arma::mat psd_project(arma::mat M);
List      bias_corrected_estimator(NumericMatrix X, int m1, int m2,
                                   bool Nh_normalize);
List      fit_to_Sigma(List fit, std::string which, bool psd);
