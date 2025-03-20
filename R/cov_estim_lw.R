#' Ledoit-Wolf Covariance Estimation (Linear Shrinkage) I
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the identity matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity.
#' Default is the optimal shrinkage intensity as in \insertCite{ledoit2003identity;textual}{covestim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the identity matrix
#' is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix,
#' s is the user-supplied or optimal shrinkage intensity and \eqn{\Sigma_{T}} is a pxp identity matrix.
#' This covariance estimator assumes a zero correlation and variances of one
#' as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under
#' \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(rets_m)
#' sigma_lwident <- cov_estim_lwident(rets_m)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_lwident
#'
cov_estim_lwident <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    p <- dim(data)[2]

    if (!zeromean_log) {
      centered <- apply(data, 2, function(x) {
        x - mean(x)
      })
      n <- dim(data)[1] - 1
    } else {
      centered <- data
      n <- dim(data)[1]
    }
    rm(data)
    gc()

    sigma_sample <- t(centered) %*% centered / n
    sigma_target <- diag(1, p)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        t(centered^2) %*% (centered^2) / n - 2 * t(centered) %*% (centered) *
          sigma_sample / n + sigma_sample^2
      asyvar <- sum(asyvar_mat)
      gamma <- sum((sigma_target - sigma_sample)^2)
      kappa <- asyvar / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }

    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    list(sigma_mat, shrink_int)
  }


#' Ledoit-Wolf Covariance Estimation (Linear Shrinkage) II
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the one-parameter matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity.
#' Default is the optimal shrinkage intensity as in \insertCite{ledoit2004oneparam;textual}{covestim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards
#' the diagonal matrix of equal variances is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and
#' \eqn{\Sigma_{T}} is a diagonal matrix with the average sample variance \eqn{\bar{\sigma}^2} on the diagonal.
#' This covariance estimator assumes a zero correlation and equal variances
#' as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under
#' \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(rets_m)
#' sigma_lwone <- cov_estim_lwone(rets_m)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_lwone
#'
cov_estim_lwone <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    p <- dim(data)[2]
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x) {
        x - mean(x)
      })
      n <- dim(data)[1] - 1
    } else {
      centered <- data
      n <- dim(data)[1]
    }
    rm(data)
    gc()

    sigma_sample <- t(centered) %*% centered / n
    sigma_target <- mean(diag(sigma_sample)) * diag(1, p)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        t(centered^2) %*% (centered^2) / n - 2 * t(centered) %*% (centered) *
          sigma_sample / n + sigma_sample^2
      asyvar <- sum(asyvar_mat)
      gamma <- sum((sigma_target - sigma_sample)^2)
      kappa <- asyvar / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }

    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    list(sigma_mat, shrink_int)
  }

#' Ledoit-Wolf Covariance Estimation (Linear Shrinkage) III
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix
#' towards the constant correlation covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity.
#' Default is the optimal shrinkage intensity as in \insertCite{ledoit2004cc;textual}{covestim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards
#' the constant correlation covariance matrix is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and
#' \eqn{\Sigma_{T}} is the constant correlation covariance matrix.
#' This covariance estimator assumes a constant correlation and individual variances
#' as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under
#' \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(rets_m)
#' sigma_lwcc <- cov_estim_lwcc(rets_m)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_lwcc
#'
cov_estim_lwcc <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    p <- dim(data)[2]
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x) {
        x - mean(x)
      })
      n <- dim(data)[1] - 1
    } else {
      centered <- data
      n <- dim(data)[1]
    }
    rm(data)
    gc()
    sigma_sample <- t(centered) %*% centered / n
    corr_mat <- stats::cor(centered)
    rbar <- (sum(colSums(corr_mat)) - p) / (p * (p - 1))
    sigma_target <- rbar * sigma_sample / corr_mat
    diag(sigma_target) <- diag(sigma_sample)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        t(centered^2) %*% (centered^2) / n - 2 * t(centered) %*% (centered) *
          sigma_sample / n + sigma_sample^2
      asyvar <- sum(asyvar_mat)
      term1 <- t((centered)^3) %*% centered
      term2 <- diag(sigma_sample) * (t(centered) %*% centered)
      term3 <-
        sigma_sample * (t(centered^2) %*% matrix(rep(1, p * nrow(centered)),
          ncol =
            p
        ))
      term4 <- (diag(sigma_sample) %o% rep(1, p)) * sigma_sample
      term_all <- (term1 - term2 - term3 + term4) / n
      ratios <-
        (diag(sigma_sample) %o% diag(sigma_sample)^-1)^0.5
      rhos <-
        0.5 * rbar * (ratios * term_all + t(ratios) * t(term_all))
      rho <-
        sum(diag(asyvar_mat), na.rm = TRUE) + sum(rhos[lower.tri(rhos)],
          na.rm =
            TRUE
        ) + sum(rhos[upper.tri(rhos)], na.rm = TRUE)
      gamma <- sum((sigma_target - sigma_sample)^2)
      kappa <- (asyvar - rho) / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }

    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    list(sigma_mat, shrink_int)
  }


#' Ledoit-Wolf Covariance Estimation (Linear Shrinkage) IV
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards
#' the one-factor covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity.
#' Default is the optimal shrinkage intensity as in \insertCite{ledoit2003identity;textual}{covestim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE).
#' Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards
#' the one-factor covariance matrix is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix,
#' s is the user-supplied or optimal shrinkage intensity and
#' \eqn{\Sigma_{T}} is the covariance matrix estimator, given by a one-factor model,
#' where the factor is equal to the cross-sectional average of all the variables.
#' This covariance estimator assumes a zero correlation and variances of one
#' as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under
#' \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(rets_m)
#' sigma_lwcc_sf <- cov_estim_lwcc_sf(rets_m)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_lwcc_sf
#'
cov_estim_lwcc_sf <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    p <- dim(data)[2]

    if (!zeromean_log) {
      centered <- apply(data, 2, function(x) {
        x - mean(x)
      })
    } else {
      centered <- data
    }
    factors <- rowMeans(centered)

    data_all <- cbind(data, factors)
    rm(data)
    gc()
    if (!zeromean_log) {
      centered_all <- apply(data_all, 2, function(x) {
        x - mean(x)
      })
      centered <- centered_all[, 1:p]
      n <- dim(data_all)[1] - 1
    } else {
      centered_all <- apply(data_all, 2, function(x) {
        x - mean(x)
      })
      centered <- centered_all[, 1:p]
      n <- dim(data_all)[1]
    }

    sigma_sample_all <- t(centered_all) %*% centered_all / n
    sigma_factors <- sigma_sample_all[1:p, p + 1]
    var_factors <- as.numeric(sigma_sample_all[p + 1, p + 1])
    sigma_sample <- sigma_sample_all[-(p + 1), -(p + 1)]
    sigma_target <- (sigma_factors %o% sigma_factors) / var_factors
    diag(sigma_target) <- diag(sigma_sample)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        (t(centered^2) %*% (centered^2)) / n - sigma_sample^
          2
      asyvar <- sum(asyvar_mat)
      rhomat_diag <- sum(diag(asyvar_mat))
      term1 <-
        1 / n * t(centered^2) %*% (centered * factors) - sigma_factors * sigma_sample
      term2 <-
        1 / n * t(centered * factors) %*% (centered * factors) - var_factors *
          sigma_sample
      rhomat1 <-
        sum(term1 * t(matrix(
          rep(sigma_factors, p),
          ncol = p, byrow = FALSE
        )) / var_factors) - sum(diag(term1) * sigma_factors / var_factors)
      rhomat2 <-
        sum(term2 * (sigma_factors %o% sigma_factors) / var_factors^2) - sum(diag(term2) *
          (sigma_factors^2) / var_factors^2)
      rho <- 2 * rhomat1 - rhomat2 + rhomat_diag
      gamma <- norm(sigma_target - sigma_sample, "F")^2
      kappa <- (asyvar - rho) / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }
    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    list(sigma_mat, shrink_int)
  }
