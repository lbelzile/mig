
#' Optimal scale matrix for kernel density estimation
#'
#' Given an \code{n} sample from a multivariate distribution on the half-space defined by
#' \eqn{\{\boldsymbol{x} \in \mathbb{R}^d: \boldsymbol{\beta}^\top\boldsymbol{x}>0\}},
#' the function computes the bandwidth (\code{type="isotropic"}) or scale
#' matrix that minimizes the asymptotic mean integrated squared error away from the boundary.
#' The latter depend on the true unknown density, which is replaced by the kernel density or
#'  a MIG distribution evaluated at the maximum likelihood estimator. The integral or the integrated
#'  squared error are obtained by Monte Carlo integration with \code{N} simulations
#'
#' @param x an \code{n} by \code{d} matrix of observations
#' @param beta \code{d} vector defining the half-space
#' @param shift location vector for translating the half-space. If missing, defaults to zero
#' @param type string indicating whether to compute an isotropic model or estimate the optimal scale matrix via optimization
#' @param family distribution for smoothing, either \code{mig} for multivariate inverse Gaussian, \code{tnorm} for truncated normal on the half-space and \code{hsgauss} for the Gaussian smoothing after suitable transformation.
#' @param method estimation criterion, either \code{amise} for the expression that minimizes the asymptotic integrated squared error, \code{lcv} for likelihood (leave-one-out) cross-validation, \code{lscv} for least-square cross-validation or \code{rlcv} for robust cross validation of Wu (2019)
#' @param N integer number of simulations for Monte Carlo integration
#' @param approx string; distribution to approximate the true density function \eqn{f(x)}; either \code{kernel} for the kernel estimator evaluated at the sample points (except for \code{method="amise"}, which isn't supported),  \code{mig} for multivariate inverse Gaussian with the method of moments or \code{tnorm} for the multivariate truncated Gaussian evaluated by maximum likelihood.
#' @param transformation string for optional scaling of the data before computing the bandwidth. Either standardization to unit variance \code{scaling}, spherical transformation to unit variance and zero correlation (\code{spherical}), or \code{none} (default).
#' @param pointwise if \code{NULL}, evaluates the mean integrated squared error, otherwise a \code{d} vector to evaluate the bandwidth or scale pointwise
#' @param maxiter integer; max number of iterations in the call to \code{optim}.
#' @param buffer double indicating the buffer from the half-space
#' @param ... additional parameters, currently ignored
#' @return a \code{d} by \code{d} scale matrix
#' @references
#' Wu, X. (2019). Robust likelihood cross-validation for kernel density estimation. \emph{Journal of Business & Economic Statistics}, 37(\bold{4}), 761–770. \doi{10.1080/07350015.2018.1424633}
#' Bowman, A.W. (1984). An alternative method of cross-validation for the smoothing of density estimates, \emph{Biometrika}, 71(\bold{2}), 353–360. \doi{10.1093/biomet/71.2.353}
#' Rudemo, M. (1982). Empirical choice of histograms and kernel density estimators. \emph{Scandinavian Journal of Statistics}, 9(\bold{2}), 65–78. http://www.jstor.org/stable/4615859
#' @export
kdens_bandwidth<- function(
      x,
      beta,
      shift,
      family = c("mig","hsgauss","tnorm"),
      method = c("amise", "lcv", "lscv", "rlcv"),
      type = c("isotropic", "diag", "full"),
      approx = c("kernel", "mig", "tnorm"),
      transformation = c("none", "scaling", "spherical"),
      N = 1e4L,
      buffer = 0,
      pointwise = NULL,
      maxiter = 2e3L,
      ...){
   maxiter <- max(as.integer(maxiter), 1e4L, na.rm = TRUE)
   transformation <- match.arg(transformation)
   approx <- match.arg(approx)
   mckern <- isTRUE(approx == "kernel")
   if(transformation %in% c("scaling","spherical")){
      if(transformation == "scaling"){
       covM <- diag(apply(x, 2, sd))
      } else{
       covM <- cov(x)
      }
      Tm <- chol(covM)
      Tminv <- solve(Tm)
      bw <- kdens_bandwidth(
         x = x %*% Tminv,
         beta = c(Tm %*% beta),
         shift = shift,
         method = method,
         type = type,
         approx = approx,
         transformation = 'none',
         N = N,
         buffer = buffer,
         pointwise = pointwise)
      return(t(Tm)%*% bw %*% Tm)
   }
   buffer <- pmax(buffer[1], 0)
   stopifnot(buffer >= 0, length(buffer) == 1, is.finite(buffer))
   method <- match.arg(method)
   approx <- match.arg(approx)
   type <- match.arg(type)
   # if(method == "amise" & approx == "empirical"){
   #    warning("Method not supported for method \"amise\".")
   #    approx <- "mig"
   # }
   N <- as.integer(N)
   stopifnot(isTRUE(is.matrix(x)))
   d <- length(beta)
   x <- as.matrix(x, ncol = d)
   n <- nrow(x)
   if(!missing(shift)){
      stopifnot(length(shift) == d)
      x <- scale(x, center = shift, scale = FALSE)
   } else{
      shift <- rep(0, d)
   }
   # Numerical optimization with Cholesky root
   chol2cov <- function(pars, d){
      stopifnot(length(pars) == d*(d+1)/2)
      chol <- matrix(0, nrow = d, ncol = d)
      diag(chol) <- abs(pars[1:d])
      chol[lower.tri(chol, diag = FALSE)] <- pars[-(1:d)]
      chol <- t(chol)
      crossprod(chol)
   }
   cov2chol <- function(cov){
     L <- t(chol(cov))
     c(abs(diag(L)), L[lower.tri(L, diag = FALSE)])
   }
   if(type == "full"){
     start <- cov2chol(n^(-1/(d+2)) * cov(x))
   } else if(type == "diag"){
      start <- log(n^(-1/(d+2)) * diag(cov(x)))
   }
   if(method == "amise"){
      if(family != "mig"){
       stop("Invalid option for the chosen \"family\".")
      }
      # Compute moment estimator
      if(approx == "mig"){
         estim <- .mig_mom(x = x, beta = beta)
      } else { # if approximation is Gaussian or truncated Gaussian
         estim <- mle_truncgauss(x = x, beta = beta)
      }
      # Simulate observations from MLE to get Monte Carlo sample
      # to approximate the integral of the MISE
      if(!is.null(pointwise)){
         stopifnot(length(pointwise) == d)
         samp <- matrix(pointwise, ncol = d)
      } else{
         if(approx == "mig"){
           samp <- rmig(n = N, xi = estim$xi,
                        Omega = estim$Omega, beta = beta)
         } else if(approx == "tnorm"){
           samp <- rtellipt(n = N, beta = beta, mu = estim$loc,
                            sigma = estim$scale, delta = buffer)
         } else{
           # approx must be 'kernel'
         stop("Invalid method") # this should not get executed
         }
      }
      xbeta <- c(samp %*% beta)
      logxbeta <- log(xbeta)
      # Integral blows up near the boundary
      # samp <- samp[xbeta > buffer,,drop = FALSE]
      # xbeta <- xbeta[xbeta > buffer]
      int1 <- exp(-(d/2)*log(4*pi) + .lsum((-d/2)*logxbeta) - log(N))
      if(type == "isotropic"){
         if(approx == "mig"){
            hopt <-  (d / n * int1 /
            mean(xbeta^2 * dmig_laplacian(x = samp, xi = estim$xi, Omega = estim$Omega, beta = beta, scale = FALSE)^2 *
                dmig(x = samp,  xi = estim$xi, Omega = estim$Omega, beta = beta, log = FALSE)))^(1/(d+4))
         } else{
            hopt <-  (d / n * int1 /
                         mean(xbeta^2 * dtnorm_laplacian(
                            x = samp, mu = estim$mu, sigma = estim$sigma, beta = beta, delta = buffer, scale = FALSE)^2 *
                                 dtellipt(x = samp,  mu = estim$mu, sigma = estim$sigma, beta = beta, delta = buffer, log = FALSE)))^(1/(d+4))
         }
      return(diag(rep(hopt, d)))
      } else {
          if(approx == "mig"){
         gradients <- mig_loglik_grad(x = samp, beta = beta, xi = estim$xi, Omega = estim$Omega)
         hessians <- mig_loglik_hessian(x = samp, beta = beta, xi = estim$xi, Omega = estim$Omega)
         logdens <- dmig(x = samp, xi = estim$xi, Omega = estim$Omega, beta = beta, log = TRUE)
         for(i in seq_len(nrow(gradients))){
            hessians[i,,] <- hessians[i,,] + tcrossprod(gradients[i,])
         }
         } else{
            Q <- solve(estim$sigma)
            gradients <- mvnorm_loglik_grad(x = samp, mu = estim$mu, Q = Q)
            hessian <- mvnorm_loglik_hessian(Q = Q)
            hessians <- array(dim = c(dim(gradients), d))
            for(i in seq_len(nrow(gradients))){
              hessians[i,,] <- hessian + tcrossprod(gradients[i,])
            }
            logdens <- dtellipt(x = samp, beta = beta, mu = estim$mu, sigma = estim$sigma, delta = buffer, log = TRUE)
         }
         logscalefact <- logdens + 2*logxbeta
         if(type == "full"){
       optfun <- function(pars, ...){
          Hmat <- chol2cov(pars, d = d)
          exp(-0.5*determinant(Hmat)$modulus - log(n)) * int1 +
             mean(exp(2*log(abs(apply(hessians, 1, function(x){sum(Hmat * x)}))) + logscalefact))/4

       }
         } else if(type == "diag"){
            optfun <- function(pars, ...){
               Hmat <- diag(exp(pars))
               exp(-0.5*determinant(Hmat)$modulus - log(n)) * int1 +
                  mean(exp(2*log(abs(apply(hessians, 1, function(x){sum(Hmat * x)}))) + logscalefact))/4

            }
         }
       # Test different optimization routines
       if (!requireNamespace("minqa", quietly = TRUE)) {
           optH <- optim(par = start,
                        fn = optfun,
                        method = "BFGS",
                        control = list(maxit = maxiter))
       } else{
       optH <- minqa::newuoa(par = start,
                             fn = optfun,
                             control = list(maxfun = maxiter))
       }
         if(type == "full"){
            return(chol2cov(optH$par, d = d))
         } else{
            return(diag(exp(optH$par)))
         }
      }
   } else if(method %in% c("lcv", "lscv", "rlcv")){
      if(family == "hsgauss"){
         Mbeta <- (diag(d) - tcrossprod(beta)/(sum(beta^2)))
         if(d > 2){
            Q2 <- t(eigen(Mbeta, symmetric = TRUE)$vectors[,-d, drop = FALSE])
         } else if(d == 2){
            Q2 <- matrix(c(-beta[2], beta[1])/sqrt(sum(beta^2)), nrow = 1) # only for d=2
         }
         # Standardize first component
         Qmat <- rbind(beta/sqrt(sum(beta^2)), Q2) # better to standardize
         # isTRUE(all.equal(diag(d), tcrossprod(Qmat)))
         map <- function(x){
            tx <- t(tcrossprod(Qmat, x))
            tx[,1] <- log(tx[,1])
            return(tx)
         }
         tx <- map(x)
         # Jacobian of transformation to Rd, akin to a weighting scheme
         logweights <- -tx[,1]
      }
      if(method == "lcv"){
        optfun2 <- function(pars, x, xsamp, dxsamp, ...){
         if(length(pars) == 1L){
            Hmat <- exp(pars[1])*diag(d)
         } else if(length(pars) == d*(d+1)/2){
            Hmat <- chol2cov(pars = pars, d = d)
         } else if(length(pars) == d){
            Hmat <- diag(exp(pars))
         }
         # -mean(sapply(seq_len(n), function(i){
         #    .lsum(dmig(x[-i, ,drop = FALSE],
         #                       xi = as.numeric(x[i,]),
         #                       Omega = Hmat, beta = beta, log = TRUE))
         # }) - log(n-1))
         if(family == "mig"){
            -mig_lcv(x = x, beta = beta, Omega = Hmat)
         } else if(family == "tnorm"){
            -tnorm_lcv(x = x, beta = beta, Omega = Hmat)
         } else if(family == "hsgauss"){
            -gauss_lcv(x = tx, Sigma = Hmat, logweights = logweights)
         }
      }
      } else { # method = least squares or robust
         covX <- cov(x)
         n <- nrow(x)
         if(approx == "mig"){
            estim <- .mig_mom(x = x, beta = beta)
         } else if(approx == "tnorm"){
            estim <- mle_truncgauss(x)
         }
         if(approx == "mig"){
            samp <- rmig(n = N, xi = estim$xi, Omega = estim$Omega, beta = beta)
            dsamp <- dmig(x = samp, xi = estim$xi, Omega = estim$Omega,
                          beta = beta, log = FALSE)
         } else if(approx == "tnorm"){
            samp <- rtellipt(n = N, beta = beta, mu = estim$loc,
                             sigma = estim$scale, delta = buffer)
            dsamp <- dtellipt(x = samp, beta = beta, mu = estim$loc,
                              sigma = estim$scale, delta = buffer, log = FALSE)
         } else{
            # Use the kernel estimator with the sample to approximate the integral
            xsamp <- matrix(0, nrow = 1, ncol = d)
            dxsamp <- 0
         }
         if(method == "rlcv"){
            a <- an(x = x)
         optfun2 <- function(pars, x, xsamp, dxsamp, ...){
            if(length(pars) == 1L){
               Hmat <- exp(pars[1])*diag(d)
            } else if(length(pars) == d*(d+1)/2){
               Hmat <- chol2cov(pars = pars, d = d)
            } else if(length(pars) == d){
              Hmat <- diag(exp(pars))
            }
            if(family == "mig"){
               - mig_rlcv(x = x, beta = beta, Omega = Hmat, xsamp = xsamp, an = a,
                              dxsamp = dxsamp, mckern = mckern)
            } else if(family == "tnorm"){
               - tnorm_rlcv(x = x, Omega = Hmat, beta = beta, an = a,
                            xsamp = xsamp, dxsamp = dxsamp, mckern = mckern)
            } else if(family == "hsgauss"){
               - gauss_rlcv(x = tx, Sigma = Hmat, logweights = logweights,
                            an = a, xsamp = xsamp, dxsamp = dxsamp, mckern = mckern)
            }
            }
         } else if(method == "lscv"){
         optfun2 <- function(pars, x, xsamp, dxsamp, ...){
            if(length(pars) == 1L){
               Hmat <- exp(pars[1])*diag(d)
            } else if(length(pars) == d*(d+1)/2){
               Hmat <- chol2cov(pars = pars, d = d)
            } else if(length(pars) == d){
               Hmat <- diag(exp(pars))
            }
            # fx <- mig_kdens_arma(x = x,
            #                      newdata = xsamp,
            #                      Omega = Hmat,
            #                      beta = beta,
            #                      logd = FALSE)
            # bias <- 0.5*mean(fx^2/dxsamp)
            #
            # checkR <- mean(exp(sapply(seq_len(n), function(i){
            #    .lsum(dmig(x[-i, ,drop = FALSE],
            #               xi = as.numeric(x[i,]),
            #               Omega = Hmat, beta = beta, log = TRUE))
            # }) - log(n-1))) - bias
            if(family == "mig"){
            -mig_lscv(x = x, beta = beta, Omega = Hmat, xsamp = xsamp,
                             dxsamp = dxsamp, mckern = mckern)
            } else if(family == "tnorm"){
               -tnorm_lscv(x = x, beta = beta, Omega = Hmat, xsamp = xsamp,
                         dxsamp = dxsamp, mckern = mckern)
            } else if(family == "hsgauss"){
               -gauss_lscv(x = x, Sigma = Hmat, logweights = logweights, xsamp = xsamp,
                           dxsamp = dxsamp, mckern = mckern)
            }
         }
      }
      }
   if(type == "isotropic"){
      convergence <- FALSE
      opt <- suppressWarnings(
         optim(fn = optfun2,
               par = (-1/(d+2))*(log(n)+log(d+2)),
               control = list(maxit = maxiter),
               method = "Brent",
               x = x,
               xsamp = samp,
               dxsamp = dsamp,
               lower = -20, upper = 5))
      if(isTRUE(opt$convergence == 0)){
            convergence <- TRUE
      }
      if(!convergence){
         warning("Optimization did not converge")
         return(diag(NA, nrow = d, ncol = d))
      } else{
         return(diag(d) * exp(opt$par))
      }
      } else {
         if (!requireNamespace("minqa", quietly = TRUE)) {
            optH <- optim(par = start,
                        fn = optfun2,
                        x = x,
                        xsamp = samp,
                        dxsamp = dsamp,
                        method = "Nelder-Mead",
                        control = list(maxit = maxiter))
        } else{
         optH <- minqa::newuoa(par = start,
                               fn = optfun2,
                               x = x,
                               xsamp = samp,
                               dxsamp = dsamp,
                               control = list(maxfun = maxiter))
         optH$convergence <- optH$ierr
        }
         if(optH$convergence != 0){
            warning("Optimization of bandwidth matrix did not converge.")
         }
         if(type == "full"){
           return(chol2cov(optH$par, d = d))
         } else if(type == "diag"){
            return(diag(exp(optH$par)))
         }
      }
   }
}

#' Multivariate inverse Gaussian kernel density estimator
#'
#' Given a matrix of new observations, compute the density of the multivariate
#' inverse Gaussian mixture defined by assigning equal weight to each component where
#' \eqn{\boldsymbol{\xi}} is the location parameter.
#' @param newdata matrix of new observations at which to evaluated the kernel density
#' @inheritParams dmig
#' @export
#' @return value of the (log)-density at \code{newdata}
mig_kdens <- function(x, newdata, Omega, beta, log = FALSE){
   # d <- length(beta)
   # newdata <- as.matrix(newdata, ncol = d)
   # x <- matrix(x, ncol = d)
   # valid <- newdata %*% beta > 0
   # logd <- rep(-Inf, nrow(newdata))
   # logd[valid] <- apply(newdata[valid,], 1, function(xi){
   #    .lsum(dmig(x = x, beta = beta, Omega = Omega, xi = xi, log = TRUE))}) - log(nrow(x))
      return(mig_kdens_arma(x = x, newdata = newdata, Omega = Omega, beta = beta, logd = log))
}



#' Truncated Gaussian kernel density estimator
#'
#' Given a data matrix over a half-space defined by \code{beta},
#' compute the log density of the asymmetric truncated Gaussian kernel density estimator,
#' taking in turn an observation as location vector.
#' @inheritParams dmig
#' @param newdata matrix of new observations at which to evaluated the kernel density
#' @param Sigma scale matrix
#' @return a vector containing the value of the kernel density at each of the \code{newdata} points
#' @export
tnorm_kdens <- function(x, newdata, Sigma, beta, log = TRUE, ...){
   tnorm_kdens_arma(x = x, newdata = newdata, Omega = Sigma, beta = beta, logd = log)
}

#' Gaussian kernel density estimator on half-space
#'
#' Given a data matrix over a half-space defined by \code{beta}, compute an homeomorphism to
#' \eqn{\mathbb{R}^d} and perform kernel smoothing based on a Gaussian kernel density estimator,
#' taking each turn an observation as location vector.
#' @inheritParams dmig
#' @param newdata matrix of new observations at which to evaluated the kernel density
#' @param Sigma scale matrix
#' @return a vector containing the value of the kernel density at each of the \code{newdata} points
#' @export
hsgauss_kdens <- function(x, newdata, Sigma, beta, log = TRUE, ...){
   beta <- as.numeric(beta)
   d <- length(beta)
   Mbeta <- (diag(d) - tcrossprod(beta)/(sum(beta^2)))
   if(d > 2){
      Q2 <- t(eigen(Mbeta, symmetric = TRUE)$vectors[,-d, drop = FALSE])
   } else if(d == 2){
      Q2 <- matrix(c(-beta[2], beta[1])/sqrt(sum(beta^2)), nrow = 1) # only for d=2
   }
   # Standardize first component
   # Standardize first component
   Qmat <- rbind(beta/sqrt(sum(beta^2)), Q2) # better to standardize
   # isTRUE(all.equal(diag(d), tcrossprod(Qmat)))
   map <- function(x){
      tx <- t(tcrossprod(Qmat, x))
      tx[,1] <- suppressWarnings(log(tx[,1]))
      return(tx)
   }
   tx <- map(x)
   # Jacobian of transformation to Rd, akin to a weighting scheme
   logjac <- -tx[,1] # det(Qmat) = 1
   tnd <- map(newdata)
   admis <- is.finite(tnd[,1])
   logd <- rep(-Inf, nrow(tnd))
   logd[admis] <- gauss_kdens_arma(x = tx, newdata = tnd[admis,,drop=FALSE],
                                   Sigma = Sigma, logweights = logjac, logd = TRUE)
   if(isTRUE(log)){
      return(logd)
   } else{
      return(exp(logd))
   }
}


#' Normal bandwidth rule
#'
#' Given an \code{n} by \code{d} matrix of observations, compute the
#' bandwidth according to Scott's rule.
#'
#' @param x \code{n} by \code{d} matrix of observations
#' @return a \code{d} by \code{d} diagonal bandwidth matrix
#' @keywords internal
#' @export
normalrule_bandwidth <- function(x){
   sigmas <- apply(x, 2, sd)
   n <- nrow(x)
   d <- ncol(x)
   (4/(d+2)/n)^(1/(d+4)) * diag(sigmas)
}

#' Threshold selection for bandwidth
#'
#' Automated thresholds selection for the robust likelihood cross validation.
#' The cutoff is based on the covariance matrix of the sample data.
#' @param x matrix of observations
#' @return cutoff for robust selection
an <- function(x){
   n <- NROW(x)
   d <- NCOL(x)
   sigma <- cov(x)
 exp(-0.5*as.numeric(determinant(sigma)$modulus) - (d/2)*log(2*pi) + lgamma(d/2) + (1-d/2)*log(log(n)) - log(n))
}
