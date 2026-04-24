#pragma once
#include <Rcpp.h>
using namespace Rcpp;

int  A_count(int n, int h, int r);
int  B_count(int n, int h, int r);
int  D_count(int n, int r);
List build_count_matrices(int n_k, bool Nh_normalize = true);
