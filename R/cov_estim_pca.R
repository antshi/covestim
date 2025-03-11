#' Principal Component Analysis Covariance Estimation
#'
#' Computes a Principal Component Analysis (PCA) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix
#' @param number_pc an integer, indicating the number of principal components. Default value is NULL and
#' the number of principal components is set according to the Marcenko-Pastur edge
#' as in \insertCite{marvcenko1967distribution;textual}{covestim}.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the number of principal components.
#' }
#' @examples
#' data(sp200)
#' sp_rets <- sp200[1:100, -1]
#' # user-defined number of factors
#' sigma_pca <- cov_estim_pca(sp_rets, number_pc = 2)[[1]]
#' # number of factors, defined with MP cut-edge
#' results_pca_mp <- cov_estim_pca(sp_rets)
#' sigma_pca_mp <- results_pca_mp[[1]]
#' number_pc <- results_pca_mp[[2]]
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{johnson2002applied}{covestim}
#'
#' \insertRef{fan2016overview}{covestim}
#'
#' @export cov_estim_pca
#'
cov_estim_pca <- function(data, number_pc = NULL) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  centered <- apply(data, 2, function(x) {
    x - mean(x)
  })
  rm(data)
  gc()
  if (is.null(number_pc)) {
    corr_mat <- stats::cor(centered)
    cor_eigen_tmp <- eigen(corr_mat)
    cor_eigenval <- cor_eigen_tmp$values
    lambdamax <- ((1 + sqrt(p / n))^2)
    index <- which(cor_eigenval >= lambdamax)
    if (length(index) != 0) {
      number_pc <- index[length(index)]
    } else {
      number_pc <- 0
    }
  }

  if (number_pc != 0) {
    centered_t <- t(centered)
    eigenvec <- eigen(t(centered_t) %*% centered_t)$vectors
    fac_pca <- sqrt(n) * eigenvec[, 1:number_pc]
    lam_pca <- centered_t %*% fac_pca / n
    uhat <- centered_t - lam_pca %*% t(fac_pca)
    low_rank <- lam_pca %*% t(lam_pca)
    su_pca <- uhat %*% t(uhat) / n
    su_diag <- diag(diag(su_pca))
    sigma_mat <- low_rank + su_diag
  } else {
    sigma_mat <- diag(1, p)
  }

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  list(sigma_mat, number_pc)
}
