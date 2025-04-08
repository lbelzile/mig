# Tests for leave-one-out cross validation
d <- 5
beta <- rexp(d)
xi <- rexp(d)
Omega <- diag(d) + matrix(1, d, d)
x <- mig::rmig(n = 100, xi = xi, Omega = Omega, beta = beta)
lcv_optfun <- function(x, beta, Omega) {
  n <- nrow(x)
  mean(sapply(seq_len(n), function(i) {
    mig::.lsum(dmig(
      x[-i, , drop = FALSE],
      xi = as.numeric(x[i, ]),
      Omega = Omega,
      beta = beta,
      log = TRUE
    ))
  })) -
    log(n - 1)
}
test <- sum(lcv_optfun(x = x, beta = beta, Omega = Omega))
test2 <- mig::mig_lcv(x = x, beta = beta, Omega = Omega)
expect_equal(test, test2)

# Check there are no errors
if (at_home()) {
  i <- 0L
  for (d in 2:3) {
    for (fam in c("mig", "hsgauss", "tnorm")) {
      for (method in c("amise", "lcv", "rlcv", "lscv")) {
        for (transfo in c("none", "scaling", "spherical")) {
          for (type in c("isotropic", "diag", "full")) {
            for (approx in c("kernel", "mig", "tnorm")) {
              error <- FALSE
              if (fam != "mig" & method == "amise") {
                error <- TRUE
              } else if (fam == "mig" & approx == "kernel") {
                error <- TRUE
              }
              beta <- rexp(d)
              xi <- rexp(d)
              Omega <- rWishart(
                n = 1,
                df = d + 4L,
                Sigma = diag(d) + matrix(1, d, d)
              )[,, 1]
              x <- mig::rmig(n = 100, xi = xi, Omega = Omega, beta = beta)
              band <- try(
                expr = mig::kdens_bandwidth(
                  x = x,
                  beta = beta,
                  family = fam,
                  method = method,
                  transformation = transfo,
                  type = type,
                  N = 1e3L,
                  approx = approx
                ),
                silent = TRUE
              )
              i <- i + 1L
              print(i)
              if (isTRUE(error)) {
                expect_inherits(current = band, class = "try-error")
              } else {
                expect_inherits(current = band, class = "matrix")
              }
            }
          }
        }
      }
    }
  }
}
