#' Simulate elliptical vector subject to a linear constraint
#'
#' Simulate multivariate Student-t \eqn{\boldsymbol{x}}
#' with location vector \code{mu}, scale matrix \code{sigma} and  \code{df} (integer) degrees of freedom
#' subject to the linear constraint \eqn{\boldsymbol{\beta}^\top\boldsymbol{x} > 0}.
#' Negative degrees of freedom or values larger than 1000 imply Gaussian vectors are generated instead.
#' @param n number of simulations
#' @param beta \code{d} vector of linear constraints
#' @param mu location vector
#' @param sigma scale matrix
#' @param delta buffer; default to zero
#' @param df degrees of freedom argument
#' @return an \code{n} by \code{d} matrix of random vectors
#' @export
#' @keywords internal
rtellipt <- function(n, beta, mu, sigma, df, delta = 0){
   d <- length(beta)
   Amat <- diag(d)
   Amat[1,] <- beta
   if(missing(df)){
      df <- 1001 # Normal
   }
   n <- as.integer(n)
   stopifnot(n > 0,
             length(df) == 1L,
             isSymmetric(sigma),
             ncol(sigma) == d,
             length(mu) == d)
   if(!isTRUE(df > 0 & df < 1000)){
      samp <- TruncatedNormal::rtmvnorm(
         n = n,
         mu = c(Amat %*% mu),
         sigma = Amat %*% sigma %*% t(Amat),
         lb = c(delta, rep(-Inf, d-1)))
   } else{
      df <- as.integer(df)
      samp <- TruncatedNormal::rtmvt(
         n = n,
         mu = c(Amat %*% mu),
         df = df,
         sigma = Amat %*% sigma %*% t(Amat),
         lb = c(delta, rep(-Inf, d-1)))
   }
   cbind(c(solve(Amat)[1,] %*% t(samp)), samp[,-1])
}

#' Density of elliptical vectors subject to a linear constraint
#'
#' Compute the density of  multivariate Student-t or Gaussian \eqn{\boldsymbol{x}}
#' with location vector \code{mu}, scale matrix \code{sigma} and  \code{df} (integer) degrees of freedom
#' subject to the linear constraint \eqn{\boldsymbol{\beta}^\top\boldsymbol{x} > \delta}.
#' Negative degrees of freedom or values larger than 1000 imply Gaussian vectors are generated instead.
#' @param n number of simulations
#' @param beta \code{d} vector of linear constraints
#' @param mu location vector
#' @param sigma scale matrix
#' @param delta buffer; default to zero
#' @param df degrees of freedom argument
#' @param log logical; if \code{TRUE}, return the log density
#' @return an \code{n} by \code{d} matrix of random vectors
#' @export
#' @keywords internal
dtellipt <- function(x, beta, mu, sigma, df, delta = 0, log = FALSE){
   d <- length(beta)
   betamu <- sum(beta*mu)
   sig_quad <- c(t(beta) %*% sigma %*% beta)
   if(missing(df)){
      df <- 1001 # Normal
   }
   x <- matrix(x, ncol = d)
   n <- nrow(x)
   stopifnot(length(df) == 1L,
             isSymmetric(sigma),
             ncol(sigma) == d,
             length(mu) == d)
   logd <- rep(-Inf, times = n)
   valid <- x %*% beta > delta
   if(!isTRUE(df > 0 & df < 1000)){
      if(sum(valid) > 0){
      logd[valid] <-  TruncatedNormal::dtmvnorm(
         x = x[valid,, drop = FALSE], mu = mu, sigma = sigma, log = TRUE) -
   TruncatedNormal::ptmvnorm(
      q = -delta,
      mu = -betamu, # symmetry of elliptical distribution about the mean
      sigma = sig_quad,
      log = TRUE)
      }
   } else{
      if(sum(valid) > 0){
      logd[valid] <-  TruncatedNormal::dtmvt(
         x = x[valid,,drop = FALSE], mu = mu, sigma = sigma, df = df,  log = TRUE) -
         TruncatedNormal::ptmvt(
            q = -delta,
            mu = -betamu, # symmetry of elliptical distribution about the mean
            sigma = sig_quad,
            df = df,
            log = TRUE)
      }
   }
   if(isTRUE(log)){
      return(logd)
   } else{
      return(exp(logd))
   }

}
