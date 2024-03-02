
library(tinytest)
d <- rpois(n = 1, lambda = 4) + 1
beta <- rexp(d)
xi <- rexp(d)
sigma <- rWishart(n = 1, df = d + 1, Sigma = diag(d) +
                     matrix(runif(1), nrow = d, ncol = d))[,,1]
sim1 <- mig::simEllipt(10, df = 2, beta = beta, mu = xi, sigma = sigma)
sim2 <- mig::simEllipt(10, beta = beta, mu = xi, sigma = sigma) # Gaussian default
sim3 <- mig::simEllipt(10, df = Inf, beta = beta, mu = xi, sigma = sigma) # Gaussian default
expect_true(isTRUE(all(sim1 %*% beta > 0)))
expect_true(isTRUE(all(sim2 %*% beta > 0)))
expect_true(isTRUE(all(sim3 %*% beta > 0)))


set.seed(1234)
samp1 <- rmig(n = 10, xi = xi, Omega = sigma, beta = beta)
set.seed(1234)
samp2 <- rmig(n = 10, xi = xi, Omega = sigma, beta = beta, shift = 1:d)
samp2
expect_equal((samp2 - samp1)[1,], 1:d)
mle1 <- mig_mle(x = samp1, beta = beta)
mle2 <- mig_mle(x = samp2, beta = beta, shift = 1:d)
mle1$xi - mle2$xi
expect_equal(mle1$Omega, mle2$Omega)
