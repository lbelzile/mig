
set.seed(1234)
d <- 5
beta <- rexp(d)
xi <- rexp(d)
Omega <- matrix(0.5, d, d) + diag(d)
samp <- rmig(n = 10, xi = xi, beta = beta, Omega = Omega)
test1 <- apply(samp, 1, function(x){
   sum(diag(numDeriv::hessian(func = dmig, x = x, xi = xi, Omega = Omega,
                              beta = beta, log = FALSE)))})
test2 <- dmig_laplacian(x = samp, xi = xi, Omega = Omega, beta = beta)
# Check Laplacian match for the density
tinytest::expect_equal(test2, test1, tolerance = 1e-3)
test3 <- array(dim = c(nrow(samp), d, d))
for(i in 1:nrow(samp)){
 test3[i,,] <- mig_loglik_hessian(x = samp[i,], xi = xi, Omega = Omega, beta = beta)
}
test4 <- mig_loglik_hessian(x = samp, beta = beta, xi = xi, Omega = Omega)
# Check hessians of log likelihood match (correct matrix multiplications)
tinytest::expect_equal(0, max(abs(test3/test4 - 1)))
# Check gradients of log likelihood
expect_equal(max(abs(mig_loglik_grad(x = samp, xi = xi, Omega = Omega, beta = beta) -
           t(apply(samp, 1, function(x){
              numDeriv::grad(func = dmig, x = x, xi = xi, Omega = Omega,
               beta = beta, log = TRUE)})))), 0)
# Check Laplacian
expect_equal(max(abs(mig_loglik_laplacian(x = samp, xi = xi, Omega = Omega, beta = beta) -
                        apply(mig_loglik_hessian(x = samp, xi = xi, Omega = Omega,
                                                 beta = beta), 1, function(y){sum(diag(y))}))), 0)
