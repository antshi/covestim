% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_help_funcs.R
\name{make_posdef}
\alias{make_posdef}
\title{Matrix Positive-Definiteness}
\usage{
make_posdef(mat, tol = 1e-08)
}
\arguments{
\item{mat}{a pxp matrix.}

\item{tol}{a double, the tolerance for the relative positiveness of eigenvalues compared to the largest eigenvalue.
Default value is tol=1e-8. If NULL, a machine tolerance approximation is applied.}
}
\value{
a positive-definite (to the tolerance of tol) matrix.
}
\description{
Produces a positive-definite matrix (to a certain tolerance).
Originally found under \insertCite{corpcorpackage;textual}{covestim}.
}
\examples{
data(rets_m)
ml_sigma <- cov_estim_wrapper2(rets_m, "ML")
is_posdef(ml_sigma)
ml_sigma_posdef <- make_posdef(ml_sigma)
is_posdef(ml_sigma_posdef)

}
\references{
\insertRef{corpcorpackage}{covestim}
}
