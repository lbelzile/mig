#' Gradient of the MIG log likelihood with respect to data
#'
#' This function returns the gradient vector of the log likelihood with respect to the
#' argument \code{x}.
#' @inheritParams dmig
#' @keywords internal
#' @export
mig_loglik_grad <- function(x, xi, Omega, beta){
   d <- length(beta)
   x <- matrix(x, ncol = d)
   stopifnot(length(xi) == d,
             nrow(Omega) == ncol(Omega),
             ncol(Omega) == d)
   Omegainv <- solve(Omega)
   if(nrow(x) == 1L){
   x <- c(x)
   betax <- sum(beta*x)
   cprod <- c(Omegainv %*% (x - xi))
   -(d/2+1) * beta / betax - cprod / betax +
      sum((x - xi) * cprod) * beta / (2 * betax^2)
   } else{
      invbetax <- 1/c(x %*% beta)
      cprodmat <- sweep(x, 2, xi) %*% Omegainv
      -(d/2+1) * tcrossprod(invbetax, beta) - sweep(cprodmat, 1, invbetax, "*") +
         tcrossprod(rowSums(sweep(x, 2, xi) * cprodmat)*invbetax^2/2, beta)
   }
}


#' Hessian of the MIG log likelihood with respect to data
#'
#' This function returns the hessian, i.e., the matrix of
#' second derivatives of the log likelihood with respect to the
#' argument \code{x}.
#' @inheritParams dmig
#' @return a \code{d} by \code{d} matrix of second derivatives if \code{x} has length \code{d},
#' else an \code{n} by \code{d} by \code{d} array if \code{x} is an \code{n} by \code{d} matrix
#' @keywords internal
#' @export
mig_loglik_hessian <- function(x, beta, xi, Omega){
   d <- length(beta)
   x <- matrix(x, ncol = d)
   Omegainv <- solve(Omega)
   if(nrow(x) == 1L){
      x <- as.vector(x)
      cprod <- c(Omegainv %*% (x-xi))
      betax <- sum(beta*x)
      ((d/2+1)*tcrossprod(beta)/betax^2 - Omegainv/betax +
            tcrossprod(cprod, beta) / betax^2 +
            tcrossprod(beta, cprod) / betax^2 -
            sum((x - xi) * cprod) * tcrossprod(beta)/(betax^3))
   } else{
      cprodmat <- sweep(x, 2, xi) %*% Omegainv
      crossprod <- outer(cprodmat, beta, FUN = "*")
      invbetax <- 1/c(x %*% beta)
      outer(invbetax^2, (d/2+1)*tcrossprod(beta), FUN = "*") -
      outer(invbetax, Omegainv, FUN = "*") +
      sweep(aperm(crossprod, c(1,3,2)) + crossprod, 1, invbetax^2, FUN = "*") -
      outer(rowSums(cprodmat * sweep(x, 2, xi)) * invbetax^3, tcrossprod(beta), FUN = "*")
   }
}


#' Laplacian of the MIG log likelihood with respect to the data
#'
#' Computes the sum of second derivatives of the multivariate
#' inverse Gaussian likelihood with respect to the data argument
#' \code{x}. The function is vectorized for more efficiency.
#'
#' @inheritParams dmig
#' @return an \code{n} vector
#' @keywords internal
#' @export
mig_loglik_laplacian <- function(x, beta, xi, Omega){
   d <- length(beta)
   x <- matrix(x, ncol = d)
   Omegainv <- solve(Omega);
   if(nrow(x) == 1L){
      x <- c(x)
      cprod <- c(Omegainv %*% (x-xi))
      betax <- sum(beta*x)
      sum(((d/2+1)*beta^2/betax^2 - diag(Omegainv)/betax +
              2*cprod*beta / betax^2 -
              sum((x - xi) * cprod) * beta^2/(betax^3)))
   } else{
      cprodmat <- sweep(x, 2, xi) %*% Omegainv
      betax <- c(x %*% beta)
      (d/2+1)*sum(beta^2)/betax^2 - sum(diag(Omegainv))/betax +
         2 * c(cprodmat %*% beta)/betax^2 -
         rowSums(cprodmat * sweep(x, 2, xi)) * sum(beta^2)/betax^3
   }
}
#' Laplacian of the MIG density with respect to the data
#'
#' Computes the sum of second derivatives of the multivariate
#' inverse Gaussian density with respect to the data argument
#' \code{x}. The function is vectorized for more efficiency.
#'
#' @inheritParams dmig
#' @return an \code{n} vector
#' @keywords internal
#' @export
dmig_laplacian <- function(x, xi, Omega, beta, scale = TRUE){
      laplacian <- rowSums(mig_loglik_grad(x = x, xi = xi, Omega = Omega, beta = beta)^2) +
         mig_loglik_laplacian(x = x, xi = xi, Omega = Omega, beta = beta)
   if(isTRUE(scale)){
      laplacian <- laplacian * dmig(x = x, xi = xi, Omega = Omega, beta = beta, log = FALSE)
   }
   return(laplacian)
}

#
# set.seed(1234)
# d <- 5
# beta <- rexp(d)
# xi <- rexp(d)
# Omega <- matrix(0.5, d, d) + diag(d)
# samp <- rmig(n = 10, xi = xi, beta = beta, Omega = Omega)
# test <- apply(samp, 1, function(x){
#    sum(diag(numDeriv::hessian(func = dmig, x = x, xi = xi, Omega = Omega,
#                               beta = beta, log = FALSE)))})
# test2 <- dmig_laplacian(x = samp, xi = xi, Omega = Omega, beta = beta)
#
# test4 <- array(dim = c(nrow(samp), d, d))
# for(i in 1:nrow(samp)){
#  test4[i,,] <- mig_loglik_hessian(x = samp[i,], xi = xi, Omega = Omega, beta = beta)
# }
# #    numDeriv::hessian(func = dmig, x = x, xi = xi, Omega = Omega, beta = beta,
# #                               log = TRUE)
# test5 <- mig_loglik_hessian(x = samp, beta = beta, xi = xi, Omega = Omega)
# max(abs(test4/test5 - 1))
