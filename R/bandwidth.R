
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
#' @param method estimation criterion, either \code{asymptotic} for the expression that minimizes the asymptotic integrated squared error or \code{lcv} for likelihood (leave-one-out) cross-validation
#' @param N integer number of simulations to evaluate the integrals of the MISE by Monte Carlo
#' @param pointwise if \code{NULL}, evaluates the mean integrated squared error, otherwise a \code{d} vector to evaluate the bandwidth or scale pointwise
#' @return a \code{d} by \code{d} scale matrix
#' @export
mig_kdens_bandwidth<- function(
      x,
      beta,
      shift,
      method = c("asymptotic", "lcv"),
      type = c("isotropic", "full"),
      N = 1e3L,
      pointwise = NULL){
   method <- match.arg(method)
   type <- match.arg(type)
   N <- as.integer(N)
    stopifnot("\"x\" is not a matrix." ~ is.matrix(x))
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
      tcrossprod(chol)
   }
   if(method == "asymptotic"){
      # Compute maximum likelihood estimate
      mle <- mig_mle(x = x, beta = beta)
      # Simulate observations from MLE to get Monte Carlo sample
      # to approximate the integral of the MISE
      if(!is.null(pointwise)){
         stopifnot(length(pointwise) == d)
         samp <- matrix(pointwise, ncol = d)
      } else{
         samp <- rmig(n = N, xi = mle$xi, Omega = mle$Omega, beta = beta)
      }
      xbeta <- c(samp %*% beta)
      int1 <- mean((4*pi*xbeta)^(-d/2))
      if(type == "isotropic"){
      hopt <-  (d / n * int1 /
      mean(xbeta^2 * dmig_laplacian(x = samp, xi = mle$xi, Omega = mle$Omega, beta = beta, scale = FALSE)^2 *
          dmig(x = samp,  xi = mle$xi, Omega = mle$Omega, beta = beta, log = FALSE)))^(2/(d+4))
      return(diag(rep(hopt, d)))
      } else if(type == "full"){
         gradients <- mig_loglik_grad(x = samp, beta = beta, xi = mle$xi, Omega = mle$Omega)
         hessians <- mig_loglik_hessian(x = samp, beta = beta, xi = mle$xi, Omega = mle$Omega)
         logdens <- dmig(x = samp, xi = mle$xi, Omega = mle$Omega, beta = beta, log = TRUE)
         logscalefact <- logdens - 2*log(xbeta) - log(4)
         for(i in seq_len(nrow(gradients))){
            hessians[i,,] <- hessians[i,,] + tcrossprod(gradients[i,])
         }

       start <- rep(0, d*(d+1)/2)
       optfun <- function(pars){
          Hmat <- chol2cov(pars, d = d)
          exp(-2*sum(pars[1:d]) - log(n))* int1 +
             mean(exp(2*log(abs(apply(hessians, 1, function(x){sum(Hmat * x)}))) + logscalefact))

       }
       optH <- optim(par = rep(0, d*(d+1)/2),
                     fn = optfun,
                     method = "Nelder-Mead",
                     control = list(maxit = 1e4L))
       return(chol2cov(optH$par, d = d))
      }
   } else if(method == "lcv"){
      objective_function <- function(pars){
         if(length(pars) == 1L){
            Hmat <- pars[1]^2*diag(d)
         } else{
            Hmat <- chol2cov(pars = pars, d = d)
         }
         -mean(sapply(seq_len(n), function(i){
            copula:::lsum(dmig(x[-i, ,drop = FALSE],
                               xi = as.numeric(x[i,]),
                               Omega = Hmat, beta = beta, log = TRUE))
         }) - log(n-1))
      }
   if(type == "isotropic"){
      opt <- try(nlm(f = objective_function, p = mean(diag(cov(x)))))
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
         return(chol2cov(optH$par, d = d))
      }
    }
}
