
library(tinytest)
d <- rpois(n = 1, lambda = 4) + 1
beta <- rexp(d)
xi <- rexp(d)
sigma <- rWishart(n = 1, df = d + 1, Sigma = diag(d) +
                     matrix(runif(1), nrow = d, ncol = d))[,,1]
sim1 <- mig::rtellipt(10, df = 2, beta = beta, mu = xi, sigma = sigma)
sim2 <- mig::rtellipt(10, beta = beta, mu = xi, sigma = sigma) # Gaussian default
sim3 <- mig::rtellipt(10, df = Inf, beta = beta, mu = xi, sigma = sigma) # Gaussian default
expect_true(isTRUE(all(sim1 %*% beta > 0)))
expect_true(isTRUE(all(sim2 %*% beta > 0)))
expect_true(isTRUE(all(sim3 %*% beta > 0)))


set.seed(1234)
samp1 <- rmig(n = 10, xi = xi, Omega = sigma, beta = beta)
set.seed(1234)
samp2 <- rmig(n = 10, xi = xi, Omega = sigma, beta = beta, shift = 1:d)
samp2
expect_equal((samp2 - samp1)[1,], 1:d)
mle1 <- mig::fit_mig(x = samp1, beta = beta, method = "mle")
mle2 <- mig::fit_mig(x = samp2, beta = beta, shift = 1:d, method = "mle")
mle1$xi - mle2$xi
expect_equal(mle1$Omega, mle2$Omega)


set.seed(1234)
d <- 2L
Omega <- rWishart(n = 1, df = d + 3, matrix(0.5, d, d) + diag(rep(0.5, d)))[,,1]
beta <- rexp(2)
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- as.numeric(pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5))
expect_equal(mc_est, p_est, tolerance = 2e-3)

set.seed(1234)
beta <- -rexp(2)
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5)
expect_equal(mc_est, p_est, tolerance = 3e-3)

set.seed(1234)
beta <- c(-rexp(1), rexp(1))
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5)
expect_equal(mc_est, p_est, tolerance = 3e-3)

set.seed(1234)
beta <- c(rexp(1), -rexp(1))
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5)
expect_equal(mc_est, p_est, tolerance = 3e-3)

set.seed(1234)
beta <- c(rexp(1), 0)
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5)
expect_equal(mc_est, p_est, tolerance = 3e-3)

set.seed(1234)
beta <- c(-rexp(1), 0)
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 2e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 2e5)
expect_equal(mc_est, p_est, tolerance = 2e-3)

set.seed(1234)
beta <- c(0, -rexp(1))
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5)
expect_equal(mc_est, p_est, tolerance = 2e-3)

set.seed(1234)
beta <- c(0, -rexp(1))
xi <- rexp(2) * sign(beta)
q <- rbind(xi - 1, xi, xi - c(0,-1), xi + c(2,0), xi + 1, xi + 5)
p_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, B = 1e5)
mc_est <- pmig(q, xi = xi, beta = beta, Omega = Omega, method = "mc", B = 1e5)
expect_equal(mc_est, p_est, tolerance = 2e-3)
