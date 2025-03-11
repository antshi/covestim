#' Wrapper Function for Covariance Estimation I
#'
#' Estimates the covariance matrix of a dataset
#' according to the user-defined function.
#'
#' @param data an nxp data matrix.
#' @param est_func an estimation function.
#' @param res_all a logical, defining the return object.
#' If FALSE, only the estimated covariance matrix is returned.
#' If TRUE, a list with two entries is returned. The first entry is the estimated covariance matrix.
#' The second entry is the estimation specific (tuning) parameter. Default value is FALSE.
#' @param ... additional parameters, parsed to est_func
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific (tuning) parameter.
#' }
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_ml <- cov_estim_wrapper(sp_rets, est_func = cov_estim_ml)[[1]]
#'
#' @export cov_estim_wrapper
#'
cov_estim_wrapper <-
  function(data, est_func, res_all = FALSE, ...) {
    result <- est_func(data, ...)

    if (res_all) {
      return(result)
    } else {
      return(result[[1]])
    }
  }

#' Wrapper Function for Covariance Estimation II
#'
#' Estimates the covariance matrix of a dataset
#' according to the user-defined character string.
#'
#' @param data an nxp data matrix.
#' @param est_type a character string, defining the estimation method.
#' @param param a double, setting an estimation specific (tuning) parameter.
#' @param factors an nxm data matrix of factors.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @param res_all a logical, defining the return object.
#' If FALSE, only the estimated covariance matrix is returned.
#' If TRUE, a list with two entries is returned. The first entry is the estimated covariance matrix.
#' The second entry is the estimation specific (tuning) parameter. Default value is FALSE.
#'
#' @return a list with following entries
#' \itemize{
#' \item a pxp estimated covariance matrix
#' \item an estimation specific tuning parameter
#' }
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_ml <- cov_estim_wrapper2(sp_rets, "ML")[[1]]
#' sigma_lwcc <- cov_estim_wrapper2(sp_rets, "LW-CC", param = 0.3, res_all = TRUE)[[1]]
#'
#' @export cov_estim_wrapper2
#'
cov_estim_wrapper2 <-
  function(data,
           est_type,
           param = NULL,
           factors = NULL,
           zeromean_log = FALSE,
           res_all = FALSE) {
    if (est_type == "SAMPLE") {
      result <- cov_estim_sample(data)
    } else if (est_type == "ML") {
      result <- cov_estim_ml(data)
    } else if (est_type == "BS") {
      result <- cov_estim_bs(data)
    } else if (est_type == "STEIN-HAFF") {
      result <- cov_estim_sh(data)
    } else if (est_type == "EVC-MP") {
      result <- cov_estim_evc_mp(data)
    } else if (est_type == "EVC-BP") {
      result <- cov_estim_evc_bp(data, param)
    } else if (est_type == "PCA") {
      result <- cov_estim_pca(data, param)
    } else if (est_type == "RIDGE") {
      result <- cov_estim_ridge(data, param)
    } else if (est_type == "LW-IDENT") {
      result <- cov_estim_lwident(data, param)
    } else if (est_type == "LW-ONE") {
      result <- cov_estim_lwone(data, param)
    } else if (est_type == "LW-CC") {
      result <- cov_estim_lwcc(data, param)
    } else if (est_type == "LW-NONLIN") {
      result <- cov_estim_lwnl(data, param)
    } else if (est_type == "POET") {
      result <- cov_estim_poet(data, param)
    } else if (est_type == "GLASSO") {
      result <- cov_estim_glasso(data, param)
    } else if (est_type == "EWMA") {
      result <- cov_estim_ewma(data, param)
    } else if (est_type == "EFM") {
      result <- cov_estim_efm(data, factors)
    } else if (est_type == "LW-CC-SF") {
      result <- cov_estim_lwcc_sf(data, param)
    } else if (est_type == "LWNL-SF") {
      result <-
        cov_estim_precond_sf(data, zeromean_log, cov_estim_lwnl,
          bandwidth_speed = param
        )
    } else if (est_type == "LWNL-AFM") {
      result <-
        cov_estim_afm(data, factors, zeromean_log, cov_estim_lwnl)
    }

    if (res_all) {
      return(result)
    } else {
      return(result[[1]])
    }
  }
