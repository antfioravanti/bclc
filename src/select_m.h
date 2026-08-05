#pragma once
#include <RcppArmadillo.h>
using namespace Rcpp;

// Data-driven selection of the truncation cut-off (defined in select_m.cpp).
List select_m(NumericMatrix X,
              std::string method,
              bool pooled,
              double kappa,
              double eta,
              int K,
              int m_max,
              bool Nh_normalize,
              IntegerVector cv_grid,
              int cv_block,
              int cv_shift,
              double cv_nugget);
