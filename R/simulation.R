
#' Simulate elliptical vector subject to a linear constraint
#'
#' Simulate multivariate Student-t \eqn{\boldsymbol{x}}
#' with location vector \code{mu}, scale matrix \code{sigma} and  \code{df} (integer) degrees of freedom
#' subject to the linear constraint \eqn{\boldsymbol{\beta}^\top\boldsymbol{x} > 0}.
#' Negative degrees of freedom or values larger than 1000 imply Gaussian vectors are generated instead
#' @param nsim number of simulations
#' @param beta \code{d} vector of linear constraints
#' @param mu location vector
#' @param sigma scale matrix
#' @param df degrees of freedom argument
#' @return an \code{nsim} by \code{d} matrix of random vectors
simEllipt <- function(nsim, beta, mu, sigma, df){
   d <- length(beta)
   Amat <- diag(d)
   Amat[1,] <- beta
   df <- as.integer(df)
   nsim <- as.integer(nsim)
   stopifnot(nsim > 0)
   if(df < 0 & df > 1000){
      samp <- TruncatedNormal::rtmvnorm(
         n = nsim,
         mu = c(Amat %*% mu),
         sigma = Amat %*% sigma %*% t(Amat),
         lb = c(0, rep(-Inf, d-1)))
   } else{
      samp <- TruncatedNormal::rtmvt(
         n = nsim,
         mu = c(Amat %*% mu),
         df = df,
         sigma = Amat %*% sigma %*% t(Amat),
         lb = c(0, rep(-Inf, d-1)))
   }
   cbind(c(solve(Amat)[1,] %*% t(samp)), samp[,-1])
}



# Use accept-reject to simulate vectors
# The distribution is marginal Gumbel, but the penultimate approximation is
# quite heavy tailed
# d <- 5
# beta <- rexp(d)
# xi <- rexp(d)
# df0 <- 10
# Omega <- matrix(0.5, d, d) + diag(d)
# accept_reject_mig <- function(n, xi, beta, Omega, df0 = 10){
#    d <- length(beta)
#    Amat <- diag(d)
#    Amat[1,] <- beta
#    logdensratio <- function(x, xi, Omega, beta, mu, df, sc, ...){
#       # Compute negative log density ratio
#       # over all potential values of x, with a scale multiple of Omega
#       - dmig(x = x, xi = xi, beta = beta, Omega = Omega, log = TRUE) +
#          TruncatedNormal::dtmvt(x, mu = mu, sigma = sc*Omega, df = df0, log = TRUE)
#    }
#    # Find the mode (multivariate Student is mu, but MIG has mean at mu, and a different mode...)
#    mode <- optim(par = xi,
#                  fn = dmig,
#                  gr = mig_loglik_grad,
#                  method = "BFGS",
#                  xi = xi,
#                  Omega = Omega,
#                  beta = beta,
#                  control = list(fnscale = -1))$par
#    gradlogdensratio <- function(x, xi, Omega, beta, mu, df, sc, ...){
#       d <- length(beta)
#       cprod <- c(solve(Omega) %*% (x - xi))
#       quad <- sum((x-xi) * cprod)
#       cprod2 <- c(solve(Omega) %*% (x - mu))
#       quad2 <- sum((x - mu) * cprod2)
#       betatx <- sum(beta*x)
#       (d/2+1) * beta/betatx + (2 * cprod * betatx - quad * beta)/(2*betatx^2) +
#        - (df + d)/(df * sc) * cprod2 / (1+quad2/(df*sc))
#    }
#    # xtest <- rexp(d)
#    # isTRUE(all.equal(
#    #    gradlogdensratio(x = xtest, xi = xi, Omega = Omega, beta = beta, mu = mode, df = df0, sc = 1.5),
#    #    numDeriv::grad(logdensratio, x = xtest, xi = xi, Omega = Omega, beta = beta, mu = mode, df = df0, sc = 1.5)
#    # ))
#    quad_betaOmega <- c(t(beta) %*% Omega %*% beta)
#    nmode <- sum(beta*mode)
#    opt_scale <- optim(
#       par = df0/((df0-2)*sum(beta*xi)),
#       method = "Brent", lower = 0.01, upper = 10,
#       control = list(fnscale = -1),
#       fn = function(scale){
#          alabama::auglag(par = xi,
#                          hin = function(x, beta, ...){sum(beta*x) },
#                fn = logdensratio,
#                gr = gradlogdensratio,
#                beta = beta, xi = xi, Omega = Omega, df = df0, mu = mode,
#                # control = list(reltol = 1e-4, abstol = 1e-3),
#                sc = scale,
#                control.outer = list(trace = FALSE))$value -
#             TruncatedNormal::ptmvt(
#                q = 0, mu = -nmode,
#                sigma = scale * quad_betaOmega,
#                df = df0,
#                log = TRUE)
#       })
#    # Compute optimal scaling of scale matrix
#    optScale <- opt_scale$par * Omega
#    # Compute minimum ratio of density
#    logC <- -opt_scale$value
#    stopifnot(logC > 0)
#    out <- matrix(nrow = 0, ncol = d)
#    nleft <- n
#    # Compute the normalizing constant for the proposal (do once)
#    # Suppress warnings because of disabled scrambling
   # normCst <- TruncatedNormal::ptmvt(
   #    q = 0,
   #    mu = -nmode, # symmetry of elliptical distribution about the mean
   #    sigma = opt_scale$par * quad_betaOmega,
   #    df = df0,
   #    log = TRUE)
   # ## Check acceptance rate via brute-force Monte Carlo
   # # mu <- replicate(n = 100, expr = {
   # #    mean((TruncatedNormal::rtmvt(
   # #    n = 1e5,
   # #    mu = mode,
   # #    df = df0,
   # #    sigma = optScale) %*% beta) >0) })
#    while(nleft > 0){
#    # Simulate proposal values from truncated elliptical Student-t
#    simCandidate <- simEllipt(
#       nsim = pmin(5e4L, ceiling(ifelse(nleft < 100, 1.5, 1.2)*nleft*exp(logC))),
#       beta = beta,
#       df = df0,
#       mu = mode,
#       sigma = optScale)
#    # Log of acceptance ratio
#    logR <- dmig(x = simCandidate, xi = xi, Omega = Omega, beta = beta, log = TRUE) -
#       TruncatedNormal::dtmvt(x = simCandidate, mu = mode, sigma = optScale,
#                              df = df0, log = TRUE) + normCst
#    stopifnot(max(logR) < logC)
#    mean(logR >= (-rexp(nrow(simCandidate)) + max(logR)))
#    1/exp(logC)
#    # Append results to existing matrix
#    out <- rbind(out, simCandidate[logR >= -rexp(nrow(simCandidate)) + logC,])
#    nleft <- n - nrow(out)
#    }
#    if(nleft < 0){
#       out <- out[seq_len(n),]
#    }
#    # Check scale matrix of sample
#    # mle <- mig_mle(out, beta = beta)
#    # mle$xi; mle$Omega
#   return(out)
# }

