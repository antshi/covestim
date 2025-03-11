#' Eigenvalue Clipping Covariance Estimation (Marcenko-Pastur)
#'
#' Computes the Eigenvalue Clipping (EVC) estimator of the covariance matrix with the Marcenko-Pastur (MP) edge.
#'
#' @param data an nxp data matrix
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#' @details  The eigenvalue clipping covariance matrix estimator is computed with the following formula:
#' \deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
#' where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and
#' \eqn{\hat{\Lambda}} is a diagonal matrix with the "clipped" sample eigenvalues.
#' The clipping procedure follows \insertCite{laloux1999;textual}{covestim}.
#' In particular, when assuming i.i.d returns, the eigenvalues of the sample correlation matrix
#' are distributed according to a Marcenko-Pastur distribution \insertCite{marvcenko1967distribution}{covestim} with
#' \deqn{\lambda_{min, max}=(1\mp\sqrt{p/n})^2}
#' as the smallest and largest eigenvalues of a random correlation matrix.
#' Therefore, only eigenvalues which lie outside this interval can bring useful information.
#' In this eigenvalue clipping procedure the sample eigenvalues bigger that \eqn{\lambda_{max}} are kept and
#' the rest are substituted with their average as in \insertCite{bouchaudpotters2009;textual}{covestim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_evc_mp <- cov_estim_evc_mp(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_evc_mp
#'
cov_estim_evc_mp <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  vola_mat <- diag(apply(data, 2, stats::sd, na.rm = TRUE))
  corr_mat <- stats::cor(data)
  eigen_tmp <- eigen(corr_mat)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  rm(data, corr_mat, eigen_tmp)
  gc()
  meanvar <- mean(eigenval)
  lambdamax <- meanvar * ((1 + sqrt(p / n))^2)
  index <- which(eigenval < lambdamax)
  eigenval_aver <- mean(eigenval[index])
  eigenval_rmt <- eigenval
  eigenval_rmt[index] <- eigenval_aver
  cut_edge <- 1 - index[1] / p

  corr_mat_rmt <- eigenvec %*% diag(eigenval_rmt) %*% t(eigenvec)

  sigma_mat <- as.matrix(vola_mat %*% corr_mat_rmt %*% vola_mat)

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  list(sigma_mat, cut_edge)
}

#'  Eigenvalue Clipping Covariance Estimation (Bouchaud-Potters)
#'
#' Computes the Eigenvalue Clipping (EVC) estimator of the covariance matrix with the Bouchaud-Potters (BP) technique.
#'
#' @param data an nxp data matrix
#' @param cut_edge a double, indicating the proportion for the applied eigenvalue clipping.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the user-supplied cut edge.
#' }
#' @details The eigenvalue clipping covariance matrix estimator is computed with the following formula:
#' \deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
#' where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and
#' \eqn{\hat{\Lambda}} is a diagonal matrix with the clipped sample eigenvalues.
#' The clipping procedure follows \insertCite{bouchaudpotters2009;textual}{covestim}.
#' In particular, the user-defined cutting edge \eqn{s} gives the proportion for the applied eigenvalue clipping
#' so that the \eqn{(1-s)\times p} largest eigenvalues are kept
#' and the remaining \eqn{s\times p} eigenvalues are substituted by their average.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[, -1]
#' sigma_evc_bp <- cov_estim_evc_bp(sp_rets, cut_edge = 0.3)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited
#'
#' @export cov_estim_evc_bp
#'
cov_estim_evc_bp <- function(data, cut_edge) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  p <- dim(data)[2]
  vola_mat <- diag(apply(data, 2, stats::sd, na.rm = TRUE))
  corr_mat <- stats::cor(data)
  eigen_tmp <- eigen(corr_mat)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  rm(data, corr_mat, eigen_tmp)

  keep_edge <- (1 - cut_edge) * p
  if (keep_edge == 0) {
    eigenval_rmt <- rep.int(1, p)
  } else if (keep_edge == p) {
    eigenval_rmt <- eigenval
  } else {
    eigenval_rmt <- eigenval
    eigenval_rmt[(keep_edge + 1):p] <-
      mean(eigenval[(keep_edge + 1):p])
  }
  corr_mat_rmt <- eigenvec %*% diag(eigenval_rmt) %*% t(eigenvec)
  sigma_mat <- as.matrix(vola_mat %*% corr_mat_rmt %*% vola_mat)

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  list(sigma_mat, cut_edge)
}
