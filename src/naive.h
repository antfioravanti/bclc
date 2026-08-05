#pragma once
#include <RcppArmadillo.h>
using namespace Rcpp;

// Naive (sample) autocovariance estimators, defined in naive.cpp.
double        naive_sample_cov(NumericMatrix X, NumericVector hvec,
                               bool Nh_normalize);
NumericMatrix build_Gamma_naive(NumericMatrix X, bool Nh_normalize);
NumericMatrix build_Sigma_naive(NumericMatrix X, bool Nh_normalize);
