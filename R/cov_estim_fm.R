#' Exact Factor Model Covariance Estimation
#'
#' Computes the Exact Factor Model (EFM) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param factors a nxf matrix with factors.
#' Default value is NULL and the factor is equal to the cross-sectional average of all the variables in the data.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here NA.
#' }
#'
#' @details The EFM covariance estimator is calculated with the following formula:
#' \deqn{\hat{\Sigma} = B\hat{\Sigma}_F B' + \hat{\Sigma_u},}
#' where \eqn{\hat{\Sigma}_F} is the sample covariance matrix of the common factors and
#' \eqn{\hat{\Sigma}_u} is the covariaance matrix of residuals, assumed to have zero correlation.
#'
#' @examples
#' data(rets_m)
#' sigma_efm <- cov_estim_efm(rets_m)[[1]]
#'
#' @export cov_estim_efm
#'
cov_estim_efm <- function(data,
                          factors = NULL,
                          zeromean_log = FALSE) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  if (is.null(factors)) {
    if (!zeromean_log) {
      factors <- as.matrix(rowMeans(apply(data, 2, function(x) {
        x - mean(x)
      })))
    } else {
      factors <- as.matrix(rowMeans(data))
    }
  }
  sigma_factors <-
    stats::var(factors, use = "pairwise", na.rm = TRUE)
  fm <- stats::lm(data ~ factors - 1)
  rm(data, factors)
  gc()
  factor_betas <- t(fm$coefficients)
  sigma_res <- diag(diag(stats::var(fm$residuals)))
  sigma_fm <-
    as.matrix(factor_betas %*% sigma_factors %*% t(factor_betas))
  sigma_mat <- as.matrix(sigma_fm + sigma_res)

  colnames(sigma_mat) <- names_data
  rownames(sigma_mat) <- names_data

  list(sigma_mat, NA)
}

#' Approximate Factor Model Covariance Estimation
#'
#' Computes the Approximate Factor Model (AFM) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param factors a nxf matrix with factors.
#' Default value is NULL and the factor is equal to the cross-sectional average of all the variables in the data.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @param resid_est_func a covariance estimation function, applied to the residuals covariance matrix.
#' @param ... further arguments to be parsed to resid_est_func
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, depending on the estimation function type.
#' }
#'
#' @details The AFM covariance estimator is calculated with the following formula:
#' \deqn{\hat{\Sigma} = B\hat{\Sigma}_F B' + \hat{\Sigma_u},}
#' where \eqn{\hat{\Sigma}_F} is the sample covariance matrix of the common factors and
#' \eqn{\hat{\Sigma}_u} is the residuals covariance matrix, estimated with the user-sapplied estim_func.
#'
#' @examples
#' data(rets_m)
#' # using the Ledoit-Wolf nonlinear shrinkage estimator
#' sigma_afm <- cov_estim_afm(rets_m, resid_est_func = cov_estim_lwnl)[[1]]
#' # using the Ledoit-Wolf linear shrinkage estimator with shrinkage intensity 0.1
#' sigma_afm <- cov_estim_afm(rets_m, resid_est_func = cov_estim_lwone, shrink_int = 0.1)[[1]]
#'
#' @export cov_estim_afm
#'
cov_estim_afm <-
  function(data,
           factors = NULL,
           zeromean_log = FALSE,
           resid_est_func,
           ...) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    if (is.null(factors)) {
      if (!zeromean_log) {
        factors <- as.matrix(rowMeans(apply(data, 2, function(x) {
          x - mean(x)
        })))
      } else {
        factors <- as.matrix(rowMeans(data))
      }
    }

    sigma_factors <-
      stats::var(factors, use = "pairwise", na.rm = TRUE)
    fm <- stats::lm(data ~ factors - 1)
    rm(data, factors)
    gc()
    sigma_fm <-
      t(fm$coefficients) %*% sigma_factors %*% t(t(fm$coefficients))
    res_estim <-
      cov_estim_wrapper(fm$residuals, res_all = TRUE, resid_est_func, ...)
    sigma_res <- res_estim[[1]]
    param_res <- res_estim[[2]]

    sigma_mat <- as.matrix(sigma_fm + sigma_res)

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    list(sigma_mat, param_res)
  }
