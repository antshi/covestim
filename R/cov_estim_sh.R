#' Stein-Haff Covariance Estimation
#'
#' Computes the Stein-Haff (SH) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#' @details The estimation procedure follows \insertCite{stein1975estimation;textual}{covestim} and
#' \insertCite{haff1991variational;textual}{covestim}.
#' Originally found under \insertCite{stcovpackage;textual}{covestim}.
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_sh <- cov_estim_sh(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_sh
#'
cov_estim_sh <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  centered <- apply(data, 2, function(x) {
    x - mean(x)
  })
  sigma_sample <- t(centered) %*% centered / (n - 1)
  eigen_tmp <- eigen(sigma_sample)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  rm(data, centered, sigma_sample, eigen_tmp)
  gc()
  k <- min(n, p)
  eigenval_subset <- eigenval[1:k]
  eigenval_sumdiffs <-
    rowSums(1 / (
      outer(eigenval_subset, eigenval_subset, "-") + diag(Inf, length(eigenval_subset))
    ))
  eigenval_stein <-
    eigenval_subset * n / (abs(n - p) + 1 + 2 * eigenval_subset * eigenval_sumdiffs)
  if (n < p) {
    eigenval_stein <- c(eigenval_stein, rep(0, p - n))
  }
  eigenval_haff <- 1 / stats::isoreg(1 / eigenval_stein[1:k])$yf
  if (n < p) {
    eigenval_haff <- c(eigenval_haff, rep(0, p - n))
  }

  sigma_mat <- eigenvec %*% diag(eigenval_haff) %*% t(eigenvec)

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  list(sigma_mat, NA)
}
