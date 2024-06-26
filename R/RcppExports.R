# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' MIG kernel density estimator
#'
#' Given a data matrix over a half-space defined by \code{beta},
#' compute the log density taking in turn an observation in  \code{newdata}
#' as location vector and computing the kernel density estimate.
#' @inheritParams dmig
#' @param newdata matrix of new observations at which to evaluated the kernel density
#' @return the value of the likelihood cross-validation criterion
#' @keywords internal
mig_kdens_arma <- function(x, newdata, Omega, beta, logd) {
    .Call(`_mig_mig_kdens_arma`, x, newdata, Omega, beta, logd)
}

tnorm_kdens_arma <- function(x, newdata, Omega, beta, logd) {
    .Call(`_mig_tnorm_kdens_arma`, x, newdata, Omega, beta, logd)
}

#' Likelihood cross-validation for kernel density estimation with MIG
#'
#' Given a data matrix over a half-space defined by \code{beta},
#' compute the log density using leave-one-out cross validation,
#' taking in turn an observation as location vector and computing the
#' density of the resulting mixture.
#' @inheritParams dmig
#' @return the value of the likelihood cross-validation criterion
#' @export
mig_lcv <- function(x, beta, Omega) {
    .Call(`_mig_mig_lcv`, x, beta, Omega)
}

#' Robust likelihood cross-validation for kernel density estimation
#'
#' Given a data matrix over a half-space defined by \code{beta},
#' compute the log density using leave-one-out cross validation,
#' taking in turn an observation as location vector and computing the
#' density of the resulting mixture.
#' @inheritParams dmig
#' @param xsamp matrix of points at which to evaluate the integral
#' @param dxsamp density of points
#' @return the value of the likelihood cross-validation criterion
#' @export
mig_rlcv <- function(x, beta, Omega, xsamp, dxsamp) {
    .Call(`_mig_mig_rlcv`, x, beta, Omega, xsamp, dxsamp)
}

mig_lscv <- function(x, beta, Omega, xsamp, dxsamp) {
    .Call(`_mig_mig_lscv`, x, beta, Omega, xsamp, dxsamp)
}

