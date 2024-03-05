
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
#' @param pointwise if \code{NULL}, evaluates the mean integrated squared error, otherwise a \code{d} vector to evaluate the bandwidth or scale pointwise
#' @return a \code{d} by \code{d} scale matrix
#' @export
mig_kdens_bandwidth<- function(
      x,
      beta,
      shift,
      method = c("amise", "lcv"),
      type = c("isotropic", "full"),
      N = 1e4L,
      # buffer = 0.01,
      pointwise = NULL){
   method <- match.arg(method)
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
   if(method == "amise"){
      # Compute maximum likelihood estimate
      mle <- .mig_mom(x = x, beta = beta)
      # Simulate observations from MLE to get Monte Carlo sample
      # to approximate the integral of the MISE
      if(!is.null(pointwise)){
         stopifnot(length(pointwise) == d)
         samp <- matrix(pointwise, ncol = d)
      } else{
         samp <- rmig(n = N, xi = mle$xi, Omega = mle$Omega, beta = beta)
      }
      xbeta <- c(samp %*% beta)
      logxbeta <- log(xbeta)
      browser()
      # Integral blows up near the boundary
      # samp <- samp[xbeta > buffer,,drop = FALSE]
      # xbeta <- xbeta[xbeta > buffer]
      int1 <- exp(-(d/2)*log(4*pi) + .lsum((-d/2)*logxbeta) - log(N))
      if(type == "isotropic"){
      hopt <-  (d / n * int1 /
      mean(xbeta^2 * dmig_laplacian(x = samp, xi = mle$xi, Omega = mle$Omega, beta = beta, scale = FALSE)^2 *
          dmig(x = samp,  xi = mle$xi, Omega = mle$Omega, beta = beta, log = FALSE)))^(1/(d+4))
      return(diag(rep(hopt, d)))
      } else if(type == "full"){
         gradients <- mig_loglik_grad(x = samp, beta = beta, xi = mle$xi, Omega = mle$Omega)
         hessians <- mig_loglik_hessian(x = samp, beta = beta, xi = mle$xi, Omega = mle$Omega)
         logdens <- dmig(x = samp, xi = mle$xi, Omega = mle$Omega, beta = beta, log = TRUE)
         logscalefact <- logdens + 2*logxbeta
         for(i in seq_len(nrow(gradients))){
            hessians[i,,] <- hessians[i,,] + tcrossprod(gradients[i,])
         }
       optfun <- function(pars){
          Hmat <- chol2cov(pars, d = d)
          exp(-sum(abs(pars[1:d])) - log(n))* int1 +
             mean(exp(2*log(abs(apply(hessians, 1, function(x){sum(Hmat * x)}))) + logscalefact)) - log(4)

       }
       cholOm <- chol(mle$Omega)
       start <- c(log(diag(cholOm)), cholOm[upper.tri(cholOm, diag = FALSE)])
        optH <- optim(par = start,#rep(0, d*(d+1)/2),
                     fn = optfun,
                     method = "BFGS",
                     control = list(maxit = 1e4L))
       return(chol2cov(optH$par, d = d))
      }
   } else if(method == "lcv"){
      optfun <- function(pars){
         if(length(pars) == 1L){
            Hmat <- pars[1]^2*diag(d)
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
      opt <- try(nlm(f = optfun, p = mean(diag(cov(x)))))
      if(!inherits(opt, "try-error")){
         convergence <- FALSE
      } else{
         if(isTRUE(opt$code %in% c(1,2))){
            return(diag(d) * opt$estimate^2)
         } else{
            convergence <- FALSE
         }
      }
      if(!convergence){
         warning("Optimization did not converge")
      }
      } else if(type == "full"){
         optH <- optim(par = rep(0, d*(d+1)/2),
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

#' Kernel density estimator
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
