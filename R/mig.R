#' Multivariate inverse Gaussian distribution
#'
#' The density of the MIG model is
#' \deqn{f(\boldsymbol{x}+\boldsymbol{a}) =(2\pi)^{-d/2}\boldsymbol{\beta}^{\top}\boldsymbol{\xi}|\boldsymbol{\Omega}|^{-1/2}(\boldsymbol{\beta}^{\top}\boldsymbol{x})^{-(1+d/2)}\exp\left\{-\frac{(\boldsymbol{x}-\boldsymbol{\xi})^{\top}\boldsymbol{\Omega}^{-1}(\boldsymbol{x}-\boldsymbol{\xi})}{2\boldsymbol{\beta}^{\top}\boldsymbol{x}}\right\}}
#' for points in the \code{d}-dimensional half-space \eqn{\{\boldsymbol{x} \in \mathbb{R}^d: \boldsymbol{\beta}^{\top}(\boldsymbol{x}-\boldsymbol{a}) \geq 0\}}
#'
#' Observations are generated using the representation as the first hitting time of a hyperplane of a correlated Brownian motion.
#' @importFrom stats cov nlm optim rnorm
#' @useDynLib mig, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @rdname mig
#' @param x \code{n} by \code{d} matrix of quantiles
#' @param log logical; if \code{TRUE}, return the log density
#' @return for \code{dmig}, the (log)-density
#'
#' @examples
#' # Density evaluation
#' x <- rbind(c(1, 2), c(2,3), c(0,-1))
#' beta <- c(1, 0)
#' xi <- c(1, 1)
#' Omega <- matrix(c(2, -1, -1, 2), nrow = 2, ncol = 2)
#' dmig(x, xi = xi, Omega = Omega, beta = beta)
#' @export
dmig <- function(x, xi, Omega, beta, shift, log = TRUE){
   # Cast vectors to matrix
   d <- length(beta)
   x <- matrix(x, ncol = d)
   if(!missing(shift)){
      stopifnot(length(shift) == d)
      x <- scale(x, center = shift, scale = FALSE)
      xi <- xi - shift
   }
   Omega <- as.matrix(Omega)
   stopifnot(length(xi) == d,
             nrow(Omega) == ncol(Omega),
             isSymmetric(Omega),
             ncol(Omega) == d,
             sum(beta*xi) > 0,
             isTRUE(all(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values > 0)))
   # Out of the halfspace if x %*% beta <= 0
   xtbeta <- as.numeric(x %*% beta)
   logdens <- rep(-Inf, length.out = nrow(x))
   insideDomain <- xtbeta > 0
   Lm <- backsolve(chol(Omega), diag(d))
   logdens[insideDomain] <-
      -d/2*log(2*pi) + sum(log(diag(Lm))) + log(sum(beta*xi)) -
      (d/2+1)*log(xtbeta[insideDomain]) -
      0.5/xtbeta[insideDomain] * apply(
         x[insideDomain,, drop = FALSE],
         MARGIN = 1,
         FUN = function(xv){
            crossprod(crossprod(Lm, xv - xi))
         }
      )
   if(isTRUE(log)){
      logdens
   } else{
      exp(logdens)
   }
}


#' Simulation of multivariate inverse Gaussian vectors
#'
#' @param n number of observations
#' @param xi \code{d} vector of location parameters \eqn{\boldsymbol{\xi}}, giving the expected value
#' @param Omega \code{d} by \code{d} positive definite scale matrix \eqn{\boldsymbol{\Omega}}
#' @param beta \code{d} vector \eqn{\boldsymbol{\beta}} defining the half-space through \eqn{\boldsymbol{\beta}^{\top}\boldsymbol{\xi}>0}
#' @param shift \code{d} translation for the half-space \eqn{\boldsymbol{a}}
#' @param timeinc time increment for multivariate simulation algorithm based on the hitting time of Brownian motion, default to \code{1e-3}.
#' @param method string; one of inverse system (\code{invsim}, default), Brownian motion (\code{bm})
#' @return for \code{rmig}, an \code{n} vector if \code{d=1} (univariate) or an \code{n} by \code{d} matrix if \code{d > 1}
#' @author Frederic Ouimet (\code{bm}), Leo Belzile (\code{invsim})
#' @rdname mig
#' @export
#' @examples
#' # Random number generation
#' d <- 5
#' beta <- runif(d)
#' xi <- rexp(d)
#' Omega <- matrix(0.5, d, d) + diag(d)
#' samp <- rmig(n = 1000, beta = beta, xi = xi, Omega = Omega)
#' mle <- mig_mle(samp, beta = beta)
rmig <- function(n, xi, Omega, beta, shift, method = c("invsim", "bm"),
                 timeinc = 1e-3){
   method <- match.arg(method)
   n <- as.integer(n[1])
   d <- length(xi)
   Omega <- as.matrix(Omega)
   if(missing(shift)){
      shift <- rep(0, d)
   }

   # Simulate the MIG vectors and store them in output_matrix
   if (d >= 2) { # when d >= 2, use the Brownian representation hitting a linear surface for the first time
   stopifnot(n >= 1,
             length(beta) == d,
             nrow(Omega) == ncol(Omega),
             isSymmetric(Omega),
             sum(beta*xi) > 0,
             ncol(Omega) == d,
             isTRUE(all(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values > 0)))
      if(method == "bm"){
   covBrownian <- t(chol(timeinc * Omega))
   output_matrix <- matrix(nrow = n, ncol = d)

      for (i in seq_len(n)) {
         # Initialize Brownian path and time
         X <- xi; t <- 0;

         # Simulate Brownian path until the linear time barrier is hit
         while (beta %*% X > t) {
            X <- X + covBrownian %*% rnorm(d)
            t <- t + timeinc
         }

         # Store the first hitting time
         output_matrix[i,] <- as.numeric(X) + shift
      }
      } else if(method == "invsim"){
         d <- length(beta)
         beta <- as.numeric(beta)
         # Project onto orthogonal complement of vector
         Mbeta <- (diag(d) - tcrossprod(beta)/(sum(beta^2)))
         # Matrix is rank-deficient: compute eigendecomposition
         # Shed matrix to remove the eigenvector corresponding to the 0 eigenvalue
         Q2 <- t(eigen(Mbeta, symmetric = TRUE)$vectors[,-d, drop = FALSE])
         # all.equal(rep(0, d-1), c(Q2 %*% beta)) # check orthogonality
         # all.equal(diag(d-1), tcrossprod(Q2)) # check basis is orthonormal
         Qmat <- rbind(beta, Q2)
         covmat <- solve(Q2 %*% solve(Omega) %*% t(Q2))

         mu <- sum(beta*xi)
         omega <- sum(beta * c(Omega %*% beta))
         Z1 <- statmod::rinvgauss(n = n, mean = mu, shape = mu^2 / omega)
         Z2 <- sweep(as.matrix(TruncatedNormal::rtmvnorm(n = n, mu = rep(0, d-1), sigma = covmat)), 1, sqrt(Z1), "*")
         Z2 <- sweep(Z2, 2, c(Q2 %*% xi), "+") + tcrossprod(Z1 - mu, c(Q2 %*% c(Omega %*% beta)/omega))
         svdQ <- svd(Qmat)
         Qinv <- svdQ$v %*% diag(1/svdQ$d) %*% t(svdQ$u)
         return(sweep(t(Qinv %*% t(cbind(Z1, Z2))), MARGIN = 2, STATS = shift, FUN = "+"))
      }
   } else { # when d = 1, use the statmod function rinvgauss
      output_matrix <- shift + statmod::rinvgauss(n = n, mean = xi, shape = xi^2 / Omega)
   }
   return(output_matrix)
}



#' Maximum likelihood estimation of multivariate inverse Gaussian vectors
#'
#' The maximum likelihood estimators are obtained for fixed shift vector \eqn{\boldsymbol{a}}
#' and half-space vector \eqn{\boldsymbol{\beta}}.
#' @inheritParams dmig
#' @return a list with components:
#' \itemize{
#' \item \code{xi}: MLE of the expectation or location vector
#' \item \code{Omega}: MLE of the scale matrix
#' }
#' @export
mig_mle <- function(x, beta, shift){
   d <- length(beta)
   x <- matrix(x, ncol = d)
   n <- nrow(x)
   if(!missing(shift)){
      stopifnot(length(shift) == d)
      x <- scale(x, center = shift, scale = FALSE)
   } else{
      shift <- rep(0, d)
   }
   xi_hat <- colMeans(x)
   Omega_hat <- matrix(0, ncol = d, nrow = d)
   for(i in seq_len(n)){
      Omega_hat <- Omega_hat + tcrossprod(x[i,] - xi_hat)/sum(beta*x[i,])
   }
   list(xi = xi_hat + shift,
        Omega = Omega_hat/n)
}


#' Method of moments estimators for multivariate inverse Gaussian vectors
#'
#' These estimators are based on the sample mean and covariance.
#' @inheritParams dmig
#' @return a list with components:
#' \itemize{
#' \item \code{xi}: MOM estimator of the expectation or location vector
#' \item \code{Omega}: MOM estimator of the scale matrix
#' }
#' @export
mig_mom <- function(x, beta, shift){
   d <- length(beta)
   x <- matrix(x, ncol = d)
   n <- nrow(x)
   if(!missing(shift)){
      stopifnot(length(shift) == d)
      x <- scale(x, center = shift, scale = FALSE)
   } else{
      shift <- rep(0, d)
   }
   xi_hat <- colMeans(x)
   Omega_hat <- cov(x)/sum(beta*xi_hat)
   list(xi = xi_hat + shift,
        Omega = Omega_hat)
}
