#ifndef LATTICEBIASCORRECTION_COVARIANCE_H
#define LATTICEBIASCORRECTION_COVARIANCE_H

//-----------------------------------------------------------------------------
// covariance.h
// Core autocovariance function for the modified exponential model.
// Defined inline so that any .cpp file can #include this header and call
// true_autocov() directly without duplication.
//
// Model:
//   C(h1,h2) = sigma2 * exp(-(|h1|^alpha1 / lambda1
//                              + |h2|^alpha2 / lambda2
//                              + beta * |h1 - h2| / lambda1))
//
// Parameters:
//   sigma2 = root of variance parameter (sigma)
//   alpha1, alpha2 = smoothness parameters (positive)
//   lambda1, lambda2 = scale parameters (positive)
//   beta = separability parameter (beta=0: separable; beta>0: non-separable)
//-----------------------------------------------------------------------------

#include <cmath>

inline double true_cov(double h1, double h2,
                           double beta,
                           double lambda1, double lambda2,
                           double alpha1, double alpha2,
                           double sigma) {
    double sigma2 = sigma * sigma;
    return sigma2 * std::exp(-(std::pow(std::abs(h1), alpha1) / lambda1 +
                             std::pow(std::abs(h2), alpha2) / lambda2 +
                             beta * std::abs(h1 - h2) / lambda1));
}

#endif // ending here LATTICEBIASCORRECTION_COVARIANCE_H
