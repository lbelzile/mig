
#' Optimal scale matrix for MIG kernel density estimation
#'
#' Given an \code{n} sample from a multivariate
#' inverse Gaussian distribution on the halfspace defined by
#' \eqn{\{\boldsymbol{x} \in \mathbb{R}^d: \boldsymbol{\beta}^\top\boldsymbol{x}>0\}},
#' the function computes the bandwidth (\code{type="isotropic"}) or scale
#' matrix that minimizes the asymptotic mean integrated squared error away from the boundary.
#' The latter depend on the true unknown density, which is replaced using as plug-in
#'  a MIG distribution evaluated at the maximum likelihood estimator. The integral or the integrated
#'  squared error are obtained by Monte Carlo integration with \code{N} simulations
#'
#' @param x an \code{n} by \code{d} matrix of observations
#' @param beta \code{d} vector defining the halfspace
#' @param shift location vector for translating the halfspace. If missing, defaults to zero
#' @param type string indicating whether to compute an isotropic model or estimate the optimal scale matrix via optimization
#' @param method estimation criterion, either \code{amise} for the expression that minimizes the asymptotic integrated squared error or \code{lcv} for likelihood (leave-one-out) cross-validation
#' @param N integer number of simulations to evaluate the integrals of the MISE by Monte Carlo
#' @param approx string; distribution to approximate the true density function \eqn{f(x)}; either \code{mig} for multivariate inverse Gaussian, or \code{tnorm} for truncated Gaussian.
#' @param transformation string for optional scaling of the data before computing the bandwidth. Either standardization to unit variance \code{scaling}, spherical transformation to unit variance and zero correlation (\code{spherical}), or \code{none} (default).
#' @param pointwise if \code{NULL}, evaluates the mean integrated squared error, otherwise a \code{d} vector to evaluate the bandwidth or scale pointwise
#' @return a \code{d} by \code{d} scale matrix
#' @export
mig_kdens_bandwidth<- function(
      x,
      beta,
      shift,
      method = c("amise", "lcv"),
      type = c("isotropic", "full"),
      approx = c("mig", "tnorm"),
      transformation = c("none","scaling","spherical"),
      N = 1e4L,
      buffer = 0.25,
      pointwise = NULL){
   transformation <- match.arg(transformation)
   if(transformation %in% c("scaling","spherical")){
      if(transformation == "scaling"){
       covM <- diag(apply(x, 2, sd))
      } else{
       covM <- cov(x)
      }
      Tm <- chol(covM)
      Tminv <- solve(Tm)
      bw <- mig_kdens_bandwidth(
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
      diag(chol) <- exp(pars[1:d])
      chol[lower.tri(chol, diag = FALSE)] <- pars[-(1:d)]
      chol <- t(chol)
      crossprod(chol)
   }
   cov2chol <- function(cov){
     L <- t(chol(cov))
     c(log(diag(L)), L[lower.tri(L, diag = FALSE)])
   }
   start <- cov2chol(n^(-1/(d+2)) * cov(x))
   if(method == "amise"){
      # Compute maximum likelihood estimate
      if(approx == "mig"){
         estim <- .mig_mom(x = x, beta = beta)
      } else{
         estim <- list(mu = colMeans(x), sigma = cov(x))
      }
      # Simulate observations from MLE to get Monte Carlo sample
      # to approximate the integral of the MISE
      if(!is.null(pointwise)){
         stopifnot(length(pointwise) == d)
         samp <- matrix(pointwise, ncol = d)
      } else{
         if(approx == "mig"){
           samp <- rmig(n = N, xi = estim$xi, Omega = estim$Omega, beta = beta)
         } else{
           samp <- rtellipt(n = N, beta = beta, mu = estim$mu,
                            sigma = estim$sigma, delta = buffer)
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
      } else if(type == "full"){
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
       optfun <- function(pars){
          Hmat <- chol2cov(pars, d = d)
          exp(-0.5*determinant(Hmat)$modulus - log(n)) * int1 +
             mean(exp(2*log(abs(apply(hessians, 1, function(x){sum(Hmat * x)}))) + logscalefact))/4

       }
        optH <- optim(par = start,
                     fn = optfun,
                     method = "BFGS",
                     control = list(maxit = 1e4L))
       return(chol2cov(optH$par, d = d))
      }
   } else if(method == "lcv"){
      optfun <- function(pars){
         if(length(pars) == 1L){
            Hmat <- exp(pars[1])*diag(d)
         } else if(length(pars) == d*(d+1)/2){
            Hmat <- chol2cov(pars = pars, d = d)
         }
         # -mean(sapply(seq_len(n), function(i){
         #    .lsum(dmig(x[-i, ,drop = FALSE],
         #                       xi = as.numeric(x[i,]),
         #                       Omega = Hmat, beta = beta, log = TRUE))
         # }) - log(n-1))
         -mig_lcv(x = x, beta = beta, Omega = Hmat)
      }
   if(type == "isotropic"){
      opt <- suppressWarnings(nlm(f = optfun,  p = (-1/(d+2))*log(n),
                                  iterlim = 1e4))
      if(inherits(opt, "try-error")){
         convergence <- FALSE
      } else{
         if(isTRUE(opt$code %in% 1:2)){
            convergence <- TRUE
         } else{
            convergence <- FALSE
            warning("Optimization did not converge")
         }
      }
      if(!convergence){
         warning("Optimization did not converge")
         return(diag(NA, nrow = d, ncol = d))
      } else{
         return(diag(d) * exp(opt$estimate))
      }
      } else if(type == "full"){
         optH <- optim(par = start, #rep(0, d*(d+1)/2),
                       fn = optfun,
                       method = "Nelder-Mead",
                       control = list(maxit = 1e4L))
         if(optH$convergence != 0){
            warning("Optimization of bandwidth matrix did not converge.")
         }
         return(chol2cov(optH$par, d = d))
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
   if(log){
      return(mig_kdens_arma(x = x, newdata = newdata, Omega = Omega, beta = beta, logd = log))
   } else{
      return(exp(mig_kdens_arma(x = x, newdata = newdata, Omega = Omega, beta = beta, logd = log)))
   }
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
#' @keywords internal
#' @export
tellipt_kdens <- function(x, newdata, Sigma, beta, log = TRUE, ...){
   #   d <- length(beta)
   #   stopifnot(ncol(x) == d,
   #             ncol(newdata) == d,
   #             ncol(Sigma) == d)
   #   valid <- newdata %*% beta > 0
   #   logd <- rep(-Inf, nrow(newdata))
   #   Sigbeta <- t(beta) %*% Sigma %*% beta
   #   logd[valid] <- apply(newdata[valid,], 1, function(mu){
   #       .lsum(TruncatedNormal::dtmvnorm(
   #         x = x, mu = mu, sigma = Sigma, log = TRUE) -
   #           TruncatedNormal::ptmvnorm(
   #         q = 0,
   #         mu = -sum(beta * mu), # symmetry of elliptical distribution about the mean
   #         sigma = Sigbeta,
   #         log = TRUE)) - log(nrow(x))
   # })
   logd <- tnorm_kdens_arma(x = x, newdata = newdata, Omega = Sigma, beta = beta, logd = log)
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

