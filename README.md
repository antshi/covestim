
# covestim

<!-- badges: start -->
<!-- badges: end -->

The package covestim includes various functions for estimating the covariance matrix of stock returns, especially in a high-dimensional setting. An exemplary dataset of 200 monthly stock returns from the S&P 500 index is provided for testing purposes.

Here is a short overview of the estimation methods:

* Sample
* Maximum-Likelihood
* EWMA
* Bayes-Stein
* Stein-Haff
* Eigenvalue Clipping
  * Marcenko-Pastur Edge
  * Bouchaud-Potters
* Linear Shrinkage
  * target with constant variance equal to 1 and correlation equal to 0
  * target with constant variance equal to the average of variances and correlation equal to 0
  * target with constant correlation equal to the average of correlations
* Nonlinear Shrinkage
  * Ledoit-Wolf
* Factor Models
  * Exact Factor Model
  * Single Market Factor Model (for stock returns)
  * Approximate Factor Model
* PCA
* POET
* Ridge Regularization
* Graphical Lasso Regularization
* t-distributed Lasso Regularization

## Installation

You can install this version of the package from github with:

``` r
install.packages("devtools")
library(devtools)
install_github("antshi/covestim")
library(covestim)
```

## Examples

These are basic examples which show you how to use the wrapper function cov_estim_wrapper. First, let's prepare with

``` r
# include the package
library(covestim)

# data
data(rets_m)
str(rets_m)
```

### Example 1

Maximum-Likelihood Estimator

```r
sigma_ml <- cov_estim_ml(rets_m)[[1]]

sigma_ml <- cov_estim_wrapper(rets_m, est_func = cov_estim_ml)

sigma_ml <- cov_estim_wrapper2(rets_m, est_type = "ML")

```

### Example 2

Linear Shrinkage towards the Constant Correlation Covariance Matrix

```r
sigma_lwcc_results <- cov_estim_lwcc(rets_m)
sigma_lwcc <- sigma_lwcc_results[[1]]
sigma_lwcc_param <- sigma_lwcc_results[[2]]

# The argument res_all = FALSE (default) in the wrapper functions cov_estim_wrapper and cov_estim_wrapper2 returns only the covariance estimator.
sigma_lwcc <- cov_estim_wrapper(rets_m, cov_estim_lwcc, res_all = FALSE)
```

You can also parse the user-defined tuning parameters such as shrinking intensity to the estimation function.

```r
sigma_lwcc_results <- cov_estim_wrapper(rets_m, cov_estim_lwcc, res_all = TRUE, shrink_int = 0.3)
str(sigma_lwcc_results)
```

### Example 3

Some more complicated application would be the Approximate Factor Model with a nonlinear shrinkage estimation on the covariance matrix of residuals

```r
sigma_afm_lwnl_results <- cov_estim_wrapper(rets_m, est_func = cov_estim_afm, resid_est_func = cov_estim_lwnl)
str(sigma_afm_lwnl_results)
```

### Example 4

While the function cov_estim_wrapper can take up any estimation function from the package, the function cov_estim_wrapper2 is interesting for practical applications.
For example, if you want to calculate different covariance estimators on the same dataset for comparison purposes, cov_estim is easier to work with.

```r
# Maximum-Likelihood, Exact Factor Model, Ledoit-Wolf linear shrinkage and Approximate Factor Model with nonlinear shrinkage on the residuals:
est_types <- c("ML", "EFM", "LW-IDENT", "LWNL-AFM") 
sigma_list <- lapply(est_types, cov_estim_wrapper2, data = rets_m)
```
