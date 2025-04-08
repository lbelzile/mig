library(tinytest)
set.seed(2024)
d <- 5L
n <- 100L
beta <- rexp(d)
xi <- rexp(d)
Omega <- matrix(0.5, d, d) + diag(d)
x <- rmig(n = n, beta = beta, xi = xi, Omega = Omega)

lcv <- sapply(1:n, function(i) {
  mig_kdens_arma(
    x = x[-i, ],
    newdata = x[i, , drop = FALSE],
    Omega = Omega,
    beta = beta,
    logd = TRUE
  )
})
lcv2 <- as.numeric(mig_loo(x = x, Omega = Omega, beta = beta))
expect_equal(lcv, lcv2)

lcv <- sapply(1:n, function(i) {
  tnorm_kdens_arma(
    x = x[-i, ],
    newdata = x[i, , drop = FALSE],
    Omega = Omega,
    beta = beta,
    logd = TRUE
  )
})
lcv2 <- as.numeric(tnorm_loo(x = x, Omega = Omega, beta = beta))
expect_equal(lcv, lcv2)


Mbeta <- (diag(d) - tcrossprod(beta) / (sum(beta^2)))
if (d > 2) {
  Q2 <- t(eigen(Mbeta, symmetric = TRUE)$vectors[, -d, drop = FALSE])
} else if (d == 2) {
  Q2 <- matrix(c(-beta[2], beta[1]) / sqrt(sum(beta^2)), nrow = 1) # only for d=2
}
# Standardize first component
Qmat <- rbind(beta / sqrt(sum(beta^2)), Q2) # better to standardize
# isTRUE(all.equal(diag(d), tcrossprod(Qmat)))
map <- function(x) {
  tx <- t(tcrossprod(Qmat, x))
  tx[, 1] <- exp(tx[, 1])
  return(tx)
}
tx <- map(x)
# Jacobian of transformation to Rd, akin to a weighting scheme
logjac <- -log(tx[, 1])
Sigma <- Qmat %*% Omega %*% t(Qmat)
lcv <- sapply(1:n, function(i) {
  gauss_kdens_arma(
    x = tx[-i, ],
    newdata = tx[i, , drop = FALSE],
    Sigma = Sigma,
    logweights = logjac[-i],
    logd = TRUE
  )
})
lcv2 <- as.numeric(gauss_loo(x = tx, Sigma = Sigma, logweights = logjac))
expect_equal(lcv, lcv2)
