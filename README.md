# bclc — Bias-Corrected Lattice Covariance Estimation

On an `n1 × n2` rectangular lattice the naive sample autocovariance `Chat(h)` is
**biased**, and the bias grows with the lag `|h|`. `bclc` removes that bias with a
linear correction: the truncated vector of sample autocovariances is multiplied by
the inverse of a deterministic **weight matrix** `W*_m` that depends only on the grid
geometry (`n1, n2, m`) — not on the data or the unknown parameters:

```
C*_m = (W*_m)^{-1} Chat_m
```

The weight matrix is derived under the modified exponential covariance model, and the
correction itself is fully non-parametric at estimation time.

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

**2) Bias-corrected estimate**, truncating lags to `{-m1,...,m1}^2` with `m1 = 3`
estimating with bias correction:

```r
fit <- bias_corrected_estimator(sim$X, m1 = 3)   # m2 defaults to m1
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

raw <- fit_to_Sigma(fit, "corrected", psd = FALSE)
raw$Sigma_hat_m                     #   same estimator, NOT projected (may be indefinite)
raw$psd_applied                     # FALSE
```

The key is always `Sigma_hat_m`: `psd` changes what that matrix *is*, never which name
to read. Never branch on `psd` to choose a field — pass the flag through and read the
one name. `psd_applied` echoes the argument back for a caller holding only the result.

Both helpers default to the corrected estimate; pass `which = "truncated"` to reshape the
uncorrected (tapered) lag vector instead. The projection is an `O(N^3)`
eigendecomposition, so pass `psd = FALSE` when the raw plug-in is what you are
studying; `psd_project()` and `is_psd()` remain exported if you need them directly.

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
- **Separable Covariance.** `beta = 0` giving a separable model and `beta > 0` 
a non-separable one.

## Main entry points

| Function | Purpose |
|---|---|
| `simulate_lattice_process()` | Simulate one Gaussian realisation on the lattice. |
| `bias_corrected_estimator()` | Main estimator: data → naive and corrected autocovariances. |
| `fit_to_Gamma()` | Fit → `(2m1+1) × (2m2+1)` lag matrix. |
| `fit_to_Sigma()` | Fit → N×N field covariance, optionally with its nearest-PSD projection (`psd`). |
| `build_count_matrices()` | 1D overlap-count matrices for one grid dimension. |
| `build_weights()` | Build `W_m` / `W*_m` from the count matrices (per-dimension cut-offs). |

## License

MIT © Antonio Fioravanti
