# bclc — Bias-Corrected Lattice Covariance Estimation

On an `n1 × n2` rectangular lattice the naive sample autocovariance `Chat(h)` is
**biased**, and the bias grows with the lag `|h|`. `bclc` removes that bias with a
linear correction: the truncated vector of sample autocovariances is multiplied by
the inverse of a deterministic **weight matrix** `W*_m` that depends only on the grid
geometry (`n1, n2, m`) — not on the data or the unknown parameters:

```
C*_m = (W*_m)^{-1} Chat_m
```

The package also ships analytical (true-parameter) MSE formulas under the modified
exponential covariance model, so the correction can be evaluated without Monte-Carlo
simulation.

## Installation

The package uses Rcpp / RcppArmadillo, so a working C++ toolchain is required
(Rtools on Windows, Xcode CLT on macOS, `r-base-dev` on Linux).

```r
# from the parent directory that contains bclc/
install.packages(c("Rcpp", "RcppArmadillo", "MASS"))   # core deps
devtools::install("bclc")
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

**2) Bias-corrected estimate**, truncating lags to `{-m,...,m}^2` with `m = 3`:

```r
fit <- bias_corrected_estimator(sim$X, m = 3)
```

This returns both the naive and the corrected autocovariances:

- `fit$Chat_m` — naive sample autocovariances (vectorised, column-major)
- `fit$Chat_star_m` — bias-corrected autocovariances (same ordering)
- `fit$W_star_m` — the weight matrix that was inverted

Both vectors are named `C(h1,h2)` in column-major lag order, so you can compare the
naive and corrected values lag by lag:

```r
head(data.frame(naive = fit$Chat_m, corrected = fit$Chat_star_m), 9)
```

## Turning the estimate into an N×N covariance matrix

`reshape_Sigma_hat_m()` reshapes the corrected lag vector into a lag-indexed matrix
and a plug-in field covariance matrix (zero outside the `m`-window):

```r
res <- reshape_Sigma_hat_m(fit$Chat_star_m, m = 3, n1 = 20, n2 = 20)

res$Gamma_hat_corr_m   # (2m+1) x (2m+1) lag matrix, rows = h1, cols = h2
res$Sigma_hat_m        # N x N covariance matrix

# The reshaped matrix is not guaranteed PSD; project if you need a valid covariance:
Sigma_psd <- psd_project(res$Sigma_hat_m)
is_psd(Sigma_psd)      # TRUE
```

## Does the correction actually help here? (no simulation needed)

Because `W*_m` and the analytical MSE are both available in closed form, you can ask
whether the correction reduces MSE for a given grid and parameter set before touching
any data. The ratio `tr(MSE_corrected) / tr(MSE_standard)` is **< 1 when the
correction helps**:

```r
mse_ratio(
  n1 = 20, n2 = 20, m = 3,
  beta = 0, lambdas = c(4, 4), alphas = c(1, 1), sigma = 1
)
# e.g. 0.78  ->  ~22% lower total MSE than the naive estimator

# Full breakdown (per-lag MSE matrices for both estimators):
cmp <- compute_mse_comparison(
  n1 = 20, n2 = 20, m = 3,
  beta = 0, lambdas = c(4, 4), alphas = c(1, 1), sigma = 1
)
cmp$scalar_MSE_standard      # tr(MSE) of the naive estimator
cmp$scalar_MSE_corrected     # tr(MSE) of the corrected estimator
```

## Choosing the truncation level m

`W*_m` must be invertible for the correction to exist. `check_invertibility()` maps
out the safe range of `m` by reporting the reciprocal condition number at each level:

```r
m1 <- build_count_matrices(20)
m2 <- build_count_matrices(20)
check_invertibility(m1, m2)   # data.frame: m, rcond, singular
```

## Notation and conventions

- **Lag matrix.** For a grid of size `n` per axis there are `L = 2n - 1` lags,
  running `{-(n-1), …, 0, …, n-1}`. The truncated window keeps `{-m,…,m}^2`.
- **`h1` is the vertical (row) lag, `h2` the horizontal (column) lag** — fixed
  everywhere.
- **Vectorisation is column-major** (R / Armadillo default): lag `(h1,h2)` sits at
  position `k = (h2 + m)(2m+1) + (h1 + m)` (0-based) in the truncated vector.
- **Normalisation.** `Nh_normalize = TRUE` (the default) uses the `N_h`
  normalisation throughout — keep it consistent between the estimator and any MSE
  comparison.

## The covariance model

All analytical results assume the **modified exponential** autocovariance:

```
C(h1, h2) = sigma^2 * exp( -( |h1|^alpha1 / lambda1
                            + |h2|^alpha2 / lambda2
                            + beta * |h1 - h2| / lambda1 ) )
```

with `beta = 0` giving a separable model and `beta > 0` a non-separable one.

## Main entry points

| Function | Purpose |
|---|---|
| `simulate_lattice_process()` | Simulate one Gaussian realisation on the lattice. |
| `bias_corrected_estimator()` | Main estimator: data → naive and corrected autocovariances. |
| `reshape_Sigma_hat_m()` | Corrected lag vector → lag matrix + N×N covariance. |
| `mse_ratio()` / `compute_mse_comparison()` | Analytical MSE comparison (no simulation). |
| `check_invertibility()` | Safe range of the truncation parameter `m`. |
| `build_weight_matrices()` | Build `W_m` / `W*_m` directly from `n1, n2, m`. |

## License

MIT © Antonio Fioravanti
