#' Bayes-Stein Covariance Estimation
#'
#' Computes the Bayes-Stein (BS) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#'
#' @details The Bayes-Stein estimator of the covariance matrix is computed according
#' to \insertCite{jorion1986bayes;textual}{covestim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_bs <- cov_estim_bs(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_bs
#'
cov_estim_bs <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  mu <- colMeans(data)
  centered <- apply(data, 2, function(x) {
    x - mean(x)
  })
  rm(data)
  gc()
  sigma_ml <- t(centered) %*% centered / n
  sigma_ml_inv <- solve(sigma_ml)
  sigma <- sigma_ml * (n / (n - p - 2))
  sigma_inv <- solve(sigma)

  ones <- rep.int(1, p)
  mug <-
    as.numeric(ones %*% sigma_ml_inv %*% mu / as.numeric(ones %*% sigma_ml_inv %*% ones))
  lambda <-
    as.numeric((p + 2) / as.numeric(t(mu - mug * ones) %*% sigma_inv %*% (mu - mug * ones)))

  sigma_mat <- (1 + 1 / (n + lambda)) * sigma + (lambda / (n * (n + 1 + lambda))) * (ones %*% t(ones)) /
    as.numeric(ones %*% sigma_inv %*% ones)

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  list(sigma_mat, NA)
}
