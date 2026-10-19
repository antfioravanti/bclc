# bclc — Bias-Corrected Lattice Covariance

On an `n1 × n2` rectangular lattice the naive sample autocovariance `Chat(h)` is
**biased** when the sample mean is used in the estimation (case of unknown true mean).
 `bclc` removes that bias with a linear correction: the truncated vector of sample autocovariances is multiplied by
the inverse of a deterministic **weight matrix** `W*_m` that depends only on the grid
geometry (`n1, n2, m`). 


The correction requires the tuning of a hyperparameter `m`. As the full weight matrix `W^*` is rank deficient, 
we adapt the work from [Vogelsang and Yang (2016)](http://dx.doi.org/10.1111/jtsa.12184) and [Li (2023)](http://dx.doi.org/10.1111/jtsa.12701) developed for time series to the spatial lattice setup, and take
a central submatrix of the weight matrix `W^*`. This submatrix is `W*_m` and `m` controls how many lags are 
included and how many are discarded. In practice, the rule is that we include all lags that are
```
|h| \leq m.
```
Then, the vector of bias-corrected sample covariances is given by
```
C*_m = (W*_m)^{-1} Chat_m
```
## Installation

The package uses Rcpp / RcppArmadillo, so a working C++ toolchain is required
(Rtools on Windows, Xcode CLT on macOS, `r-base-dev` on Linux).

```r
install.packages(c("Rcpp", "RcppArmadillo", "MASS", "remotes"))   # core deps
remotes::install_github("antfioravanti/bclc")
library(bclc)
```
## Example

**1) Simulate one Gaussian realisation on a 20 × 20 lattice** — a separable
exponential covariance (`beta = 0`) with medium spatial dependence (`lambda = 4`):

```r
library(bclc)

sim <- simulate_lattice_process(
  n1 = 20, n2 = 20,
  sigma   = 1,
  alpha1  = 1, alpha2  = 1,
  lambda1 = 4, lambda2 = 4,
  beta    = 0,
  seed    = 1
)
```

The result holds the data matrix and the parameters used to generate it:

- `sim$X` — the 20 × 20 data matrix
- `sim$Sigma` — the true N × N covariance used to generate `X`
- `sim$params` — the parameters above

**2) Bias-corrected estimate**, truncating lags to `|h| \leq m` with `m1 = 3`
estimating with bias correction:

```r
fit <- bias_corrected_estimator(sim$X, m1 = 3)   # m2=m1 by default
```

This returns both the naive and the corrected autocovariances:

- `fit$Chat_m` — naive sample autocovariances (vectorised, column-major)
- `fit$Chat_star_m` — bias-corrected autocovariances (same ordering)
- `fit$W_star_m` — the weight matrix that was inverted
- `fit$m1`, `fit$m2`, `fit$n1`, `fit$n2` — the cut-offs and grid dimensions

Both vectors are named `C(h1,h2)` in column-major lag order, so you can compare the
naive and corrected values lag by lag:

```r
head(data.frame(naive = fit$Chat_m, corrected = fit$Chat_star_m), 9)
```

## Turning the estimate into an N×N covariance matrix

Two helpers take the fit object directly. `fit_to_Gamma()` reshapes the corrected lag
vector into a lag-indexed matrix; `fit_to_Sigma()` builds the plug-in field covariance
matrix (zero outside the `m`-window) together with its nearest-PSD projection:

```r
Gamma <- fit_to_Gamma(fit)          # (2m1+1) x (2m2+1) lag matrix, rows = h1, cols = h2

Sig <- fit_to_Sigma(fit)            # psd = TRUE by default
Sig$Sigma_hat_m                     #   N x N field covariance, PSD-projected
is_psd(Sig$Sigma_hat_m)             # TRUE

Sigma_hat <- fit_to_Sigma(fit, "corrected", psd = FALSE)
Sigma_hat$Sigma_hat_m                     
Sigma_hat$psd_applied                     
```

## Main entry points

| Function | Purpose |
|---|---|
| `simulate_lattice_process()` | Simulate one Gaussian realisation on the lattice. |
| `bias_corrected_estimator()` | Main estimator: data to naive and corrected autocovariances. |
| `fit_to_Gamma()` | Fit to `(2m1+1) × (2m2+1)` lag matrix. |
| `fit_to_Sigma()` | Fit to N×N field covariance, optionally with its nearest-PSD projection (`psd`). |
| `build_count_matrices()` | Count of overlap matrices for one grid dimension. |
| `build_weights()` | Build `W_m` / `W*_m` from the count matrices. |

## License

MIT © Antonio Fioravanti
