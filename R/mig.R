#' Multivariate inverse Gaussian distribution
#'
#' The density of the MIG model is
#' \deqn{f(\boldsymbol{x}+\boldsymbol{a}) =(2\pi)^{-d/2}\boldsymbol{\beta}^{\top}\boldsymbol{\xi}|\boldsymbol{\Omega}|^{-1/2}(\boldsymbol{\beta}^{\top}\boldsymbol{x})^{-(1+d/2)}\exp\left\{-\frac{(\boldsymbol{x}-\boldsymbol{\xi})^{\top}\boldsymbol{\Omega}^{-1}(\boldsymbol{x}-\boldsymbol{\xi})}{2\boldsymbol{\beta}^{\top}\boldsymbol{x}}\right\}}
#' for points in the \code{d}-dimensional half-space \eqn{\{\boldsymbol{x} \in \mathbb{R}^d: \boldsymbol{\beta}^{\top}(\boldsymbol{x}-\boldsymbol{a}) \geq 0\}}
#'
#' Observations are generated using the representation as the first hitting time of a hyperplane of a correlated Brownian motion.
#' @importFrom stats cov nlm optim rnorm sd runif
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
dmig <- function(x, xi, Omega, beta, shift, log = FALSE){
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
#' d <- 5L
#' beta <- runif(d)
#' xi <- rexp(d)
#' Omega <- matrix(0.5, d, d) + diag(d)
#' samp <- rmig(n = 1000, beta = beta, xi = xi, Omega = Omega)
#' mle <- fit_mig(samp, beta = beta, method = "mle")
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
         # Q2 <- t(svd(Mbeta)$u[,-d, drop = FALSE])
         # all.equal(rep(0, d-1), c(Q2 %*% beta)) # check orthogonality
         # all.equal(diag(d-1), tcrossprod(Q2)) # check basis is orthonormal
         Qmat <- rbind(beta, Q2)
         covmat <- solve(Q2 %*% solve(Omega) %*% t(Q2))

         mu <- sum(beta*xi)
         omega <- sum(beta * c(Omega %*% beta))
         Z1 <- statmod::rinvgauss(n = n, mean = mu, shape = mu^2 / omega)
         Z2 <- sweep(as.matrix(TruncatedNormal::rtmvnorm(
            n = n, mu = rep(0, d-1), sigma = covmat)), 1, sqrt(Z1), "*")
         Z2 <- sweep(Z2, 2, c(Q2 %*% xi), "+") +
            tcrossprod(Z1 - mu, c(Q2 %*% c(Omega %*% beta)/omega))
         svdQ <- svd(Qmat)
         Qinv <- svdQ$v %*% diag(1/svdQ$d) %*% t(svdQ$u)
         return(sweep(t(Qinv %*% t(cbind(Z1, Z2))),
                      MARGIN = 2, STATS = shift, FUN = "+"))
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
#' @keywords internal
.mig_mle <- function(x, beta, shift){
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
#' \item \code{xi}: MOM estimate of the expectation or location vector
#' \item \code{Omega}: MOM estimate of the scale matrix
#' }
#' @export
#' @keywords internal
.mig_mom <- function(x, beta, shift){
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

#' Fit multivariate inverse Gaussian distribution
#' @inheritParams .mig_mle
#' @param method string, one of \code{mle} for maximum likelihood estimation, or \code{mom} for method of moments.
#' @return a list with components:
#' \itemize{
#' \item \code{xi}: estimate of the expectation or location vector
#' \item \code{Omega}: estimate of the scale matrix
#' }
#' @export
fit_mig <- function(x, beta, method = c("mle","mom"), shift){
   method <- match.arg(method)
   if(method == "mle"){
      .mig_mle(x = x, beta = beta, shift = shift)
   } else if(method == "mom"){
      .mig_mom(x = x, beta = beta, shift = shift)
   }
}


rtinvgauss <- function(n, lower = 0, upper = Inf, xi, Omega){
   n <- as.integer(n)
   stopifnot(n >= 1, length(upper) == 1L, isTRUE(upper > 0),
             length(lower) == 1L, isTRUE(lower >= 0), lower <= upper)
   if(lower == upper){ return(rep(lower, n))}
   lo <- ifelse(lower == 0, 0, statmod::pinvgauss(lower, mean = xi, shape = xi^2/Omega))
   up <- ifelse(is.finite(upper), statmod::pinvgauss(upper, mean = xi, shape = xi^2/Omega), 1)
   statmod::qinvgauss(lo + (up-lo)*runif(n),
                      mean = xi, shape = xi^2/Omega)
}


#' Distribution function of multivariate inverse Gaussian vectors
#'
#' @param q \code{n} by \code{d} matrix of quantiles
#' @param log logical; if \code{TRUE}, returns log probabilities
#' @param B number of Monte Carlo replications for the SOV estimator
#' @inheritParams dmig
#' @return an \code{n} vector of (log) probabilities
#' @author Leo Belzile
#' @rdname mig
#' @export
#' @examples
#' set.seed(1234)
#' d <- 2L
#' beta <- runif(d)
#' Omega <- rWishart(n = 1, df = 2*d, Sigma = matrix(0.5, d, d) + diag(d))[,,1]
#' xi <- rexp(d)
#' q <- mig::rmig(n = 10, beta = beta, Omega = Omega, xi = xi)
#' pmig(q, xi = xi, beta = beta, Omega = Omega)
pmig <- function(q, xi, Omega, beta, log = FALSE, method = c("sov", "mc"), B = 1e4L){
   d <- length(xi)
   method <- match.arg(method)
   log <- isTRUE(log)
   if(d == 1L){
      stopifnot(length(Omega) == 1L)
      return(statmod::pinvgauss(as.vector(q),
                                mean = xi, shape = xi^2/Omega, log = log))
   }
   if((d > 2) & (method == "sov")){
      method <- "mc"
      # warning("Separation-of-variables estimator not implemented beyond bivariate case.")
   }
   if(method == "mc"){
      n_batch <- ceiling(B / 1e5)
      mc <- matrix(nrow = nrow(q), ncol = n_batch)
      for(i in seq_len(n_batch)){
         mcsamp <- mig::rmig(n = min(B, 1e5), beta = beta, Omega = Omega, xi = xi)
         mc[,i] <- as.numeric(apply(q, 1, function(x){
            mean(apply(sweep(mcsamp, 2, x, FUN = "-"), 1, max) < 0)}))
      }
      est <- rowMeans(mc)
      if(n_batch >= 5){
        attr(est, "sd") <- apply(mc, 1, sd)/sqrt(n_batch)
      }
      return(est)
   }
   beta <- as.numeric(beta)
   Mbeta <- (diag(d) - tcrossprod(beta)/(sum(beta^2)))
   Q2 <- t(eigen(Mbeta, symmetric = TRUE)$vectors[,-d, drop = FALSE])
   Q2 <- matrix(c(-beta[2], beta[1])/sqrt(sum(beta^2)), nrow = 1) # only for d=2
   Qmat <- rbind(beta, Q2)
   covmat <- solve(Q2 %*% solve(Omega) %*% t(Q2))
   # Mean and scale of the univariate inverse Gaussian
   xi_r <- sum(beta*xi)
   omega_r <- c(t(beta) %*% Omega %*% beta)
   mu_zp1 <- c(Q2 %*% xi)
   mu_zp2 <- c(Q2 %*% Omega %*% beta / omega_r)
   q <- matrix(q, ncol = d)
   n <- nrow(q)
   B <- as.integer(B)
   # TODO reordering based on partial information
   # about bounds to reduce the error
   L <- t(chol(covmat))
   # Container for log probability
   lprobs <- rep(-Inf, n)
   negBetas <- isTRUE(any(beta < 0))
   for(i in seq_len(n)){
      rmax <- Inf
      rmin <- 0
      tq <- sum(beta * q[i,]) # as.numeric(q %*% beta)
      if(!negBetas){ # all positive, so q necessarily in the half space
         rmax <- max(0, tq)
      }
      if(isTRUE(all(beta <= 0))){
         rmin <- max(0, tq)
      }
      if(rmax > 0){
         if(d > 2L){
            Z <- matrix(0, nrow = d - 2, ncol = B)
         }
         # create array for variables
         R <- rtinvgauss(n = B, lower = rmin, upper = rmax, xi = xi_r, Omega = omega_r)
         p <- rep(0, B)
         tv1 <- Qmat %*% c(q[i,1], -beta[1]/beta[2]*q[i,1])
         tv2 <- Qmat %*% c(-beta[2]/beta[1]*q[i,2], q[i,2])
         if(isTRUE(all(beta > 0))){
            zmin <- Q2[2]/beta[2]*R + tv1[2]
            zmax <- Q2[1]/beta[1]*R + tv2[2]
         } else if((beta[1] < 0) & (beta[2] > 0)){
            zmax <- Inf
            zmin <- pmax(Q2[1]/beta[1]*R + tv2[2], Q2[2]/beta[2]*R + tv1[2])
         } else if((beta[1] > 0) & (beta[2] < 0)){
            zmin <- -Inf
            zmax <- pmin(Q2[1]/beta[1]*R + tv2[2], Q2[2]/beta[2]*R + tv1[2])
         } else if(isTRUE(all(beta < 0))){
            zmax <- Q2[2]/beta[2]*R + tv1[2]
            zmin <- Q2[1]/beta[1]*R + tv2[2]
         } else if(beta[1] == 0 & beta[2] > 0){
            zmax <- Inf
            zmin <- -q[i,1]
            rmin = 0
            rmax = tq
         } else if(beta[2] == 0 & beta[1] > 0){
            zmin <- -Inf
            zmax <- q[i,2]
            rmin = 0
            rmax = tq
         }  else if(beta[1] == 0 & beta[2] < 0){
            zmax <- q[i,1]
            zmin <- -Inf
            rmin <- tq
            rmax <- Inf
         } else if(beta[2] == 0 & beta[1] < 0){
            zmax <- Inf
            zmin <- -q[i,2]
            rmin <- tq
            rmax <- Inf
         } else{
            stop("Invalid input for \"beta\"")
         }
         upp <- zmax
         low <- zmin
         # compute likelihood ratio for R
         # In d >= 2, we would need to solve linear programs for each point
         for(k in seq_len(d-1)){
            #    capture.output(low_bd <- limSolve::linp(
            #       E = Q[1:k,],
            #       F = c(r[b], Z[b]),
            #       G = rbind(beta, -diag(d)),
            #       H = c(0, -q[i,]),
            #       Cost = Q2[k,],
            #       ispos = FALSE, verbose = FALSE),  file = nullfile())
            #    if(isTRUE(low_bd$IsError)){
            #       low <- -Inf
            #    } else{
            #       low <- low_bd$solutionNorm
            #    }
            #    capture.output(upp_bd <- limSolve::linp(G = rbind(beta, -diag(d)),
            #                             H = c(0, -q[i,]),
            #                             Cost = -Q2[k,],
            #                             ispos = FALSE, verbose = FALSE),  file = nullfile())
            #    if(isTRUE(upp_bd$IsError)){
            #       upp <- Inf
            #    } else{
            #       upp <- -upp_bd$solutionNorm
            #    }
            # print(integrate(function(r){
            #    exp(TruncatedNormal::lnNpr((low - low/rmax*r)/sqrt(r)/L[k,k], (upp - upp/rmax*r)/sqrt(r)/L[k,k])) *
            #       statmod::dinvgauss(x = r, mean = xi_r, shape = xi_r^2/omega_r)},  lower = 0, upper = rmax))
            # compute matrix multiplication L*Z
            if(k > 1){
               col <- c(L[k,seq_len(k)] %*% Z[seq_len(k),])
               # compute limits of truncation
               tu <- ((upp - (mu_zp1[k] + mu_zp2[k]*(R-xi_r)))/sqrt(R) - col)/L[k,k]
               tl <- ((low - (mu_zp1[k] + mu_zp2[k]*(R-xi_r)))/sqrt(R) - col)/L[k,k]
            } else{
               tu <- (upp - (mu_zp1[k] + mu_zp2[k]*(R-xi_r)))/(sqrt(R)*L[k,k])
               tl <- (low - (mu_zp1[k] + mu_zp2[k]*(R-xi_r)))/(sqrt(R)*L[k,k])
            }
            #simulate N(0, 1) conditional on [tl, tu]
            if(k < (d-1)){
               # final Z(d) need not be simulated
               Z[k,] <- TruncatedNormal::trandn(tl, tu);
            }
            # update likelihood ratio
            # browser()
            p <- p + TruncatedNormal::lnNpr(tl, tu)
         }
         # now switch back from logarithmic scale
         lprobs[i] <- mig::.lsum(p) +
            log(statmod::pinvgauss(rmax, mean = xi_r, shape = xi_r^2/omega_r) -
                   statmod::pinvgauss(rmin, mean = xi_r, shape = xi_r^2/omega_r))
      }
   }
   lprobs <- lprobs - log(B)
   if(log){
      return(as.numeric(lprobs))
   } else{
      return(exp(as.numeric(lprobs)))
   }
}
