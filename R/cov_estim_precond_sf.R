#' Covariance Estimation after Preconditioning with a Single Factor Model
#'
#' Computes a specific estimator of the covariance matrix after preconditioning
#' the data with a Single Factor (SF) model - the implicit market portfolio.
#'
#' @param data an nxp data matrix.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @param precond_est_func a function for estimating the precondtioned covariance matrix.
#' @param ... further arguments to be parsed to precond_est_func
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the bandwidth speed.
#' }
#'
#' @examples
#' data(rets_m)
#' sigma_lwnl_sf <- cov_estim_precond_sf(rets_m,
#'   precond_est_func = cov_estim_lwnl, bandwidth_speed = NULL
#' )[[1]]
#' sigma_glasso_sf <- cov_estim_precond_sf(rets_m,
#'   precond_est_func = cov_estim_glasso, rho = 0.01
#' )[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_precond_sf
#'
cov_estim_precond_sf <-
  function(data,
           zeromean_log = FALSE,
           precond_est_func = NULL,
           ...) {
    data <- as.matrix(data)
    names_data <- colnames(data)

    sqrt_sigma_mat_fm <-
      sqrt_root_mat_calc(cov_estim_efm(data, zeromean_log = zeromean_log)[[1]])
    data_precond <- data %*% solve(sqrt_sigma_mat_fm)
    rm(data)
    gc()
    results_precond <-
      cov_estim_wrapper(
        data = data_precond,
        est_func = precond_est_func,
        res_all = TRUE,
        ...
      )
    sigma_mat_precond <- results_precond[[1]]
    param_precond <- results_precond[[2]]

    sigma_mat <-
      sqrt_sigma_mat_fm %*% sigma_mat_precond %*% sqrt_sigma_mat_fm

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    list(sigma_mat, param_precond)
  }
