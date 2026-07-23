#pragma once
#include <RcppArmadillo.h>
using namespace Rcpp;

int  A_count(int n, int h, int r);
int  B_count(int n, int h, int r);
int  D_count(int n, int r);
List build_count_matrices(int n_k);

arma::mat build_Pm(int L, int m);
arma::mat build_Pmc(int L, int m);
arma::mat build_Sm(int L1, int L2, int m1, int m2);
List build_weights(List mat1, List mat2, int m1, int m2, bool Nh_normalize);
