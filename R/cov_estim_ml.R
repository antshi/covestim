#' Maximum-Likelihood Covariance Estimation
#'
#' Computes the Maximum-Likelihood (ML) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#'
#' @details The Maximum-Likelihood estimator of the covariance matrix for a data matrix X
#' is computed with the following formula:
#' \deqn{\hat{\Sigma}=\frac{1}{n} \left(X - \widehat{{\mu}} {1} \right)' \left({X} -  \widehat{{\mu}}{1}\right)}
#' where \eqn{\mu=\bar{x}_{j}=\frac{1}{n}\sum_{i=1}^{n}x_{ij}} for (for \eqn{i=1,\ldots, n}
#' and \eqn{j=1,\ldots,p}) is the sample mean vector and \eqn{1} is an 1xp vector of ones.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_ml <- cov_estim_ml(sp_rets)[[1]]
#'
#' @export cov_estim_ml
#'
cov_estim_ml <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  centered <- apply(data, 2, function(x) {
    x - mean(x)
  })
  sigma_mat <- t(centered) %*% centered / n

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  list(sigma_mat, NA)
}
