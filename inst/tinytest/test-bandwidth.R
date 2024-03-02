
d <- 5;
beta <- rexp(d);
xi <- rexp(d);
Omega <- diag(d) + matrix(1, d, d)
x <- mig::rmig(n = 100, xi = xi, Omega = Omega, beta = beta)
lcv_optfun <- function(x, beta, Omega){
   n <- nrow(x)
   mean(sapply(seq_len(n), function(i){
      mig::.lsum(dmig(x[-i, ,drop = FALSE],
                 xi = as.numeric(x[i,]),
                 Omega = Omega, beta = beta, log = TRUE))
   })) - log(n-1)
}
test <- sum(lcv_optfun(x = x, beta = beta, Omega = Omega))
test2 <- mig::mig_lcv(x = x, beta = beta, Omega = Omega)
expect_equal(test, test2)
