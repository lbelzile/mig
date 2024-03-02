#' Simulate elliptical vector subject to a linear constraint
#'
#' Simulate multivariate Student-t \eqn{\boldsymbol{x}}
#' with location vector \code{mu}, scale matrix \code{sigma} and  \code{df} (integer) degrees of freedom
#' subject to the linear constraint \eqn{\boldsymbol{\beta}^\top\boldsymbol{x} > 0}.
#' Negative degrees of freedom or values larger than 1000 imply Gaussian vectors are generated instead
#' @param n number of simulations
#' @param beta \code{d} vector of linear constraints
#' @param mu location vector
#' @param sigma scale matrix
#' @param df degrees of freedom argument
#' @return an \code{n} by \code{d} matrix of random vectors
#' @export
#' @keywords internal
simEllipt <- function(n, beta, mu, sigma, df){
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
         lb = c(0, rep(-Inf, d-1)))
   } else{
      df <- as.integer(df)
      samp <- TruncatedNormal::rtmvt(
         n = n,
         mu = c(Amat %*% mu),
         df = df,
         sigma = Amat %*% sigma %*% t(Amat),
         lb = c(0, rep(-Inf, d-1)))
   }
   cbind(c(solve(Amat)[1,] %*% t(samp)), samp[,-1])
}

