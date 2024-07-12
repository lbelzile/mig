#' Gradient of the MIG log likelihood with respect to data
#'
#' This function returns the gradient vector of the log likelihood with respect to the
#' argument \code{x}.
#' @inheritParams dmig
#' @keywords internal
#' @return an \code{n} by \code{d} matrix of first derivatives for the gradient, observation by observation, or a \code{d} vector if \code{x} is a vector.
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

mvnorm_loglik_grad <- function(x, mu, Q){
   d <- ncol(Q)
   x <- matrix(x, ncol = d)
   -scale(x, center = mu, scale = FALSE) %*% Q
}

mvnorm_loglik_hessian <- function(Q){
   -Q
}

mvnorm_loglik_laplacian <- function(Q){
   -sum(diag(Q))
}

dtnorm_laplacian <- function(x, mu, sigma, beta, delta = 0, scale = TRUE){
   Q <- solve(sigma)
   laplacian <- rowSums(mvnorm_loglik_grad(x = x, mu = mu, Q = Q)^2) +
      mvnorm_loglik_laplacian(Q = Q)
   if(isTRUE(scale)){
      laplacian <- laplacian * dtellipt(x = x, beta = beta, mu = mu, sigma = sigma, df = 0, delta = delta, log = FALSE)
   }
   return(laplacian)
}

#' Log of sum of terms
#'
#' Computes the log of a sum of positive components, given on the log scale (\code{lx}), avoiding numerical overflow.
#' @param lx vector or matrix of log of terms to add
#' @param loff log of offset
#' @keywords internal
#' @author Marius Hofert, Martin Maechler (package \code{copula})
#' @return a vector of log sums
#' @export
.lsum <- function (lx, loff) {
   d <- dim(lx)
   rx <- length(d)
   if (mis.off <- missing(loff))
      loff <- {
         if (rx <= 1L)
            max(lx)
         else if (rx == 2L)
            apply(lx, 2L, max)
      }
   if (rx <= 1L) {
      if (is.finite(loff))
         loff + log(sum(exp(lx - loff)))
      else if (mis.off || is.na(loff) || loff == max(lx))
         loff
      else stop("'loff  is infinite but not == max(.)")
   }
   else if (rx == 2L) {
      if (any(x.off <- !is.finite(loff))) {
         if (mis.off || isTRUE(all.equal(loff, apply(lx,
                                                      2L, max)))) {
            if (all(x.off))
               return(loff)
            r <- loff
            iok <- which(!x.off)
            l.of <- loff[iok]
            r[iok] <- l.of + log(colSums(exp(lx[, iok, drop = FALSE] -
                                                rep(l.of, each = d[1]))))
            r
         }
         else stop("'loff' has non-finite values but differs from default max(.)")
      }
      else loff + log(colSums(exp(lx - rep(loff, each = d[1]))))
   }
   else stop("not yet implemented for arrays of rank >= 3")
}


#' Magnetic storms
#'
#' Absolute magnitude of 373 geomagnetic storms lasting more than 48h with absolute magnitude (dst) larger than 100 in 1957-2014.
#'
#' @source Aki Vehtari
#' @references World Data Center for Geomagnetism, Kyoto, M. Nose, T. Iyemori, M. Sugiura, T. Kamei (2015), \emph{Geomagnetic Dst index}, doi:10.17593/14515-74000.
#' @docType data
#' @note For a detailed article presenting the derivation of the Dst index, see \code{http://wdc.kugi.kyoto-u.ac.jp/dstdir/dst2/onDstindex.html}
#' @format a vector of size 373
#' @name geomagnetic
NULL
