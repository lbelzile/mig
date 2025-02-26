// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' MIG kernel density estimator
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the log density taking in turn an observation in  \code{newdata}
//' as location vector and computing the kernel density estimate.
//' @inheritParams dmig
//' @param newdata matrix of new observations at which to evaluated the kernel density
//' @return the value of the likelihood cross-validation criterion
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector mig_kdens_arma(arma::mat x, arma::mat newdata, arma::mat Omega, arma::colvec beta, bool logd){
 arma::uword  n = newdata.n_rows;
 arma::uword  d = beta.n_elem;
 double maxele = 0;
 double lognx = log(x.n_rows);
 double betaxi = 0;
 if((newdata.n_cols != d) | (x.n_cols != d)){
   Rcpp::stop("Invalid data: \"newdata\" and \"x\" must have the same number of columns as the length of \"beta\".");
 }
 if((Omega.n_rows != d) | (Omega.n_cols != d)){
   Rcpp::stop("Invalid dimensions for the scale matrix.");
 }
 arma::colvec betax = x * beta;
 if(!all(betax > 0)){
   Rcpp::stop("Invalid \"x\" matrix input: some entries are not points in the half space.");
 }
 arma::colvec temp_container(x.n_rows);
 arma::mat cholinv = chol(inv_sympd(Omega));
 // Container for log density
 Rcpp::NumericVector logdens(n);
 double logcst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag()));
 for(arma::uword i=0; i < n; i++){
   betaxi = sum(newdata.row(i) * beta);
   if(betaxi > 0){
     for(arma::uword j=0; j < x.n_rows; j++){
       temp_container(j) = -0.5*sum(pow(cholinv * (newdata.row(i) - x.row(j)).t(), 2))/betax(j) -(0.5*d+1.0)*log(betax(j));
     }
     maxele = max(temp_container);
     logdens[i] = logcst + log(betaxi) + maxele + log(sum(exp(temp_container - maxele)))- lognx;
     temp_container.zeros();
   } else{
     logdens[i] = - arma::datum::inf;
   }
 }
 if(!logd){
   logdens = exp(logdens);
 }
 return logdens;
}

//' Truncated Gaussian kernel density estimator
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the log density of the asymmetric truncated Gaussian kernel density estimator,
//' taking in turn an observation as location vector.
//' @inheritParams dmig
//' @param newdata matrix of new observations at which to evaluated the kernel density
//' @return the value of the likelihood cross-validation criterion
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector tnorm_kdens_arma(arma::mat x, arma::mat newdata, arma::mat Omega, arma::colvec beta, bool logd){
 arma::uword  n = newdata.n_rows;
 arma::uword  d = beta.n_elem;
 double maxele = 0;
 double lognx = log(x.n_rows);
 double betaxi = 0;
 arma::vec sdTr = sqrt(trans(beta) * Omega * beta);
 if((newdata.n_cols != d) | (x.n_cols != d)){
   Rcpp::stop("Invalid data: \"newdata\" and \"x\" must have the same number of columns as the length of \"beta\".");
 }
 if((Omega.n_rows != d) | (Omega.n_cols != d)){
   Rcpp::stop("Invalid dimensions for the scale matrix.");
 }
 arma::colvec betax = x * beta;
 if(!all(betax > 0)){
   Rcpp::stop("Invalid \"x\" matrix input: some entries are not points in the half space.");
 }
 arma::colvec temp_container(x.n_rows);
 arma::mat cholinv = chol(inv_sympd(Omega));
 // Container for log density
 Rcpp::NumericVector logdens(n);
 double logcst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag()));
 for(arma::uword i=0; i < n; i++){
   betaxi = sum(newdata.row(i) * beta);
   if(betaxi > 0){
     for(arma::uword j=0; j < x.n_rows; j++){
       temp_container(j) = -0.5*sum(pow(cholinv * (newdata.row(i) - x.row(j)).t(), 2));
     }
     maxele = max(temp_container);
     logdens[i] = logcst + maxele + log(sum(exp(temp_container - maxele))) - lognx - log(arma::normcdf(0.0, -betaxi, sdTr(0)));
     temp_container.zeros();
   } else{
     logdens[i] = - arma::datum::inf;
   }
 }
 if(!logd){
   logdens = exp(logdens);
 }
 return logdens;
}

//' Gaussian kernel density estimator
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the log density of the asymmetric truncated Gaussian kernel density estimator,
//' taking in turn an observation as location vector.
//' @param x matrix of observations
//' @param newdata matrix of new observations at which to evaluated the kernel density
//' @param Sigma covariance matrix
//' @param logd logical; if \code{TRUE}, return the log density
//' @return log density estimator
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector gauss_kdens_arma(arma::mat x, arma::mat newdata, arma::mat Sigma, arma::vec logweights, bool logd){
 arma::uword  d = x.n_cols;
 double maxele = 0;
 double lognx = log(x.n_rows);
 if(newdata.n_cols != d){
   Rcpp::stop("Invalid data: \"newdata\" and \"x\" must have the same number of columns.");
 }
 if((Sigma.n_rows != d) | (Sigma.n_cols != d)){
   Rcpp::stop("Invalid dimensions for the scale matrix.");
 }
 if(logweights.n_elem != x.n_rows){
    Rcpp::stop("Invalid weight vector.");
 }
 arma::colvec temp_container(x.n_rows);
 arma::mat cholinv = chol(inv_sympd(Sigma));
 // Container for log density
 Rcpp::NumericVector logdens(newdata.n_rows);
 double logcst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag())) - lognx;
 for(arma::uword i=0; i < newdata.n_rows; i++){
   for(arma::uword j=0; j < x.n_rows; j++){
     temp_container(j) = -0.5*sum(pow(cholinv * (newdata.row(i) - x.row(j)).t(), 2)) + logweights(j);
   }
   maxele = max(temp_container);
   logdens[i] = logcst + maxele + log(sum(exp(temp_container - maxele)));
   temp_container.zeros();
 }
 if(!logd){
   logdens = exp(logdens);
 }
 return logdens;
}


//' Leave-one-out cross-validation for kernel density estimation with MIG
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the log density using leave-one-out cross validation,
//' taking in turn an observation as location vector and computing the
//' density of the resulting mixture.
//' @inheritParams dmig
//' @return the value of the likelihood cross-validation criterion
//' @export
//' @keywords internal
// [[Rcpp::export]]
arma::colvec mig_loo(arma::mat x, arma::mat Omega, arma::colvec beta) {
 arma::mat cholinv = chol(inv_sympd(Omega));
 arma::uword d = x.n_cols;
 arma::uword n = x.n_rows;
 if(n < 2){
   Rcpp::stop("Invalid \"x\" argument: must have at least two rows.");
 }
 arma::uword k = 0;
 // arma::mat quad = arma::mat(n, n);
 arma::colvec quadv(n*(n-1)/2);
 arma::colvec objective(n);
 //double objective = 0;
 arma::uword im = 0;
 arma::uword jm = 0;
 arma::uword ind = 0;
 //arma::colvec objective(n);
 arma::colvec terms(n-1);
 double maxterms = 0;

 // double term = 0;
 if(beta.n_rows != d){
   Rcpp::stop("Invalid arguments: \"beta\" must be the same length as the number of columns of \"x\".");
 }
 arma::colvec cp = arma::colvec(n);
 k = 0;
 for(arma::uword i = 0; i < n-1; i++){
   for(arma::uword j = i+1; j < n; j++){ // assumes that n>1
     quadv(k) = sum(pow(cholinv * (x.row(i) - x.row(j)).t(), 2));
     k++;
   }
   cp(i) = sum(beta.t() % x.row(i));
 }
 cp(n-1) =  sum(beta.t() % x.row(n-1));
 double cst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag())) - log(n-1.0);
 for(arma::uword i = 0; i < n; i++){
   terms.zeros();
   k = 0;
   for(arma::uword j = 0; j < n; j++){
     if(i != j){
       im = std::min(i,j);
       jm = std::max(i,j);
       ind = im*n-im*(im+1)/2+(jm-im)-1;
       terms(k) = -(0.5*d+1.0)*log(cp(j))-0.5*quadv(ind)/cp(j);
       k++;
     }
   }
   maxterms = max(terms);
   objective(i) += maxterms + log(sum(exp(terms - maxterms))) + cst + log(cp(i));
 }
 return objective;
}



//' Robust likelihood cross-validation for kernel density estimation for MIG
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the log density using leave-one-out cross validation,
//' taking in turn an observation as location vector and computing the
//' density of the resulting mixture.
//' @inheritParams dmig
//' @param xsamp matrix of points at which to evaluate the integral
//' @param dxsamp density of points
//' @return the value of the likelihood cross-validation criterion
//' @export
//' @keywords internal
// [[Rcpp::export]]
double mig_rlcv(arma::mat x, arma::colvec beta, arma::mat Omega, double an,
               arma::mat xsamp, arma::vec dxsamp, bool mckern = true) {
 arma::mat cholinv = chol(inv_sympd(Omega));
 arma::uword d = x.n_cols;
 arma::uword n = x.n_rows;
 if(n < 2){
   Rcpp::stop("Invalid \"x\" argument: must have at least two rows.");
 }
 arma::uword k = 0;
 // arma::mat quad = arma::mat(n, n);
 arma::colvec quadv(n*(n-1)/2);
 double objective = 0;
 arma::uword im = 0;
 arma::uword jm = 0;
 arma::uword ind = 0;
 //arma::colvec objective(n);
 arma::colvec terms(n-1);
 double logobj = 0;
 double maxterms = 0;
 //double logan = -0.5*d*log(2*arma::datum::pi) - 0.5*log_det_sympd(cov(x)) + lgamma(0.5*d) + (1-0.5*d)*log(log(1.0 * n)) - log(1.0 * n);
 // double an = exp(logan);
 double logan = log(an);
 // double term = 0;
 if(beta.n_rows != d){
   Rcpp::stop("Invalid arguments: \"beta\" must be the same length as the number of columns of \"x\".");
 }
 arma::colvec cp = arma::colvec(n);
 k = 0;
 for(arma::uword i = 0; i < n-1; i++){
   for(arma::uword j = i+1; j < n; j++){ // assumes that n>1
     quadv(k) = sum(pow(cholinv * (x.row(i) - x.row(j)).t(), 2));
     k++;
   }
   cp(i) = sum(beta.t() % x.row(i));
 }
 cp(n-1) =  sum(beta.t() % x.row(n-1));
 arma::colvec logcp = log(cp);
 // Dropping constant terms
 double logcst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag())) - log(n-1.0);
 for(arma::uword i = 0; i < n; i++){
   terms.zeros();
   k = 0;
   for(arma::uword j = 0; j < n; j++){
     if(i != j){
       im = std::min(i,j);
       jm = std::max(i,j);
       ind = im*n-im*(im+1)/2+(jm-im)-1;
       terms(k) = -(0.5*d+1.0)*logcp(j)-0.5*quadv(ind)/cp(j);
       k++;
     }
   }
   maxterms = max(terms);
   // Obtain the log penalty
   logobj = maxterms + log(sum(exp(terms - maxterms))) + logcst + logcp(i);
   if(logobj >= logan){
     objective += logobj;
   } else{ // Modify log if x < a
     objective += logan - 1.0 + exp(logobj - logan);
   }
 }
 // Compute bias term
 double bias = 0;
 if(mckern){
   Rcpp::NumericVector fx = mig_kdens_arma(x, x, Omega, beta, false);
   for(arma::uword i = 0; i < xsamp.n_rows; i++){
     if(fx(i) >= an){
       bias += 1;
     } else{
       bias += fx(i)*0.5/an;
     }
   }
 } else{
   Rcpp::NumericVector fx = mig_kdens_arma(x, xsamp, Omega, beta, false);
   for(arma::uword i = 0; i < xsamp.n_rows; i++){
     if(fx(i) >= an){
       bias += fx(i)/dxsamp(i);
     } else{
       bias += fx(i)*fx(i)*0.5/(dxsamp(i)*an);
     }
   }
 }
 bias = bias / ((double) xsamp.n_rows);
 objective = objective/n - bias;
 return objective;
}



//' Least squares cross-validation for MIG density estimation
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the average using leave-one-out cross validation density minus
//' half the squared density.
//' @inheritParams dmig
//' @inheritParams mig_rlcv
//' @return the value of the least square cross-validation criterion
//' @export
//' @keywords internal
// [[Rcpp::export]]
double mig_lscv(arma::mat x, arma::colvec beta, arma::mat Omega,
               arma::mat xsamp, arma::vec dxsamp, bool mckern = true) {
 arma::mat cholinv = chol(inv_sympd(Omega));
 arma::uword d = x.n_cols;
 arma::uword n = x.n_rows;
 if(n < 2){
   Rcpp::stop("Invalid \"x\" argument: must have at least two rows.");
 }
 arma::uword k = 0;
 // arma::mat quad = arma::mat(n, n);
 arma::colvec quadv(n*(n-1)/2);
 double objective = 0;
 arma::uword im = 0;
 arma::uword jm = 0;
 arma::uword ind = 0;
 //arma::colvec objective(n);
 arma::colvec terms(n-1);
 double logobj = 0;
 double maxterms = 0;
 // double term = 0;
 if(beta.n_rows != d){
   Rcpp::stop("Invalid arguments: \"beta\" must be the same length as the number of columns of \"x\".");
 }
 arma::colvec cp = arma::colvec(n);
 k = 0;
 for(arma::uword i = 0; i < n-1; i++){
   for(arma::uword j = i+1; j < n; j++){ // assumes that n>1
     quadv(k) = sum(pow(cholinv * (x.row(i) - x.row(j)).t(), 2));
     k++;
   }
   cp(i) = sum(beta.t() % x.row(i));
 }
 cp(n-1) =  sum(beta.t() % x.row(n-1));
 arma::colvec logcp = log(cp);
 // We are dropping constant terms here
 // Namely -0.5*d*log(2*pi) + log(cp(i)) +
 double logcst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag())) - log(n-1.0);
 for(arma::uword i = 0; i < n; i++){
   terms.zeros();
   k = 0;
   for(arma::uword j = 0; j < n; j++){
     if(i != j){
       im = std::min(i,j);
       jm = std::max(i,j);
       ind = im*n-im*(im+1)/2+(jm-im)-1;
       terms(k) = -(0.5*d+1.0)*logcp(j)-0.5*quadv(ind)/cp(j);
       k++;
     }
   }
   maxterms = max(terms);
   // Obtain the log penalty
   logobj = maxterms + log(sum(exp(terms - maxterms)));
   objective += exp(logobj + logcst + logcp(i));
 }
 // Compute squared density term using Monte Carlo draws xsamp
 // with density dxsamp
 double quadT = 0;
 if(mckern){
   for(arma::uword i = 0; i < xsamp.n_rows; i++){
     Rcpp::NumericVector fx = mig_kdens_arma(x, x, Omega, beta, false);
     quadT += fx(i);
   }
 } else{
   Rcpp::NumericVector fx = mig_kdens_arma(x, xsamp, Omega, beta, false);
   for(arma::uword i = 0; i < xsamp.n_rows; i++){
     quadT += fx(i)*fx(i)/dxsamp(i);
   }
 }
 quadT = quadT * 0.5 / xsamp.n_rows;
 objective = objective/n - quadT;
 return objective;
}



//' Leave-one-out cross-validation for Gaussian kernel density estimation
//'
//' Given a data matrix, compute the log density using leave-one-out
//' cross validation, taking in turn an observation as location vector
//' and computing the density of the resulting mixture.
//' @param x \code{n} by \code{d} matrix of observations
//' @param Sigma smoothing positive-definite matrix
//' @param logweights vector of log weights
//' @return a vector of values for the weighted leave-one-out likelihood
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::colvec gauss_loo(arma::mat x, arma::mat Sigma, arma::vec logweights) {
 arma::mat cholinv = chol(inv_sympd(Sigma));
 arma::uword d = x.n_cols;
 arma::uword n = x.n_rows;
 if(n < 2){
   Rcpp::stop("Invalid \"x\" argument: must have at least two rows.");
 }
 if(logweights.n_elem != n){
    Rcpp::stop("Invalid \"logweights\" argument: must have the same length as the number of rows of \"x\".");
 }
 arma::uword k = 0;
 // arma::mat quad = arma::mat(n, n);
 arma::colvec quadv(n*(n-1)/2);
 arma::colvec objective(n);
 // double objective = 0;
 arma::uword im = 0;
 arma::uword jm = 0;
 arma::uword ind = 0;
 //arma::colvec objective(n);
 arma::colvec terms(n-1);
 double maxterms = 0;
 double cst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag())) - log(n-1.0);
 // double term = 0;
 k = 0;
 for(arma::uword i = 0; i < n-1; i++){
   for(arma::uword j = i+1; j < n; j++){ // assumes that n>1
     quadv(k) = sum(pow(cholinv * (x.row(i) - x.row(j)).t(), 2));
     k++;
   }
 }
 for(arma::uword i = 0; i < n; i++){
   terms.zeros();
   k = 0;
   for(arma::uword j = 0; j < n; j++){
     if(i != j){
       im = std::min(i,j);
       jm = std::max(i,j);
       ind = im*n-im*(im+1)/2+(jm-im)-1;
       terms(k) = -0.5*quadv(ind) + logweights(j);
       k++;
     }
   }
   maxterms = max(terms);
   objective(i) += maxterms + log(sum(exp(terms - maxterms))) + cst;
 }
 return objective;
}


//' Leave-one-out cross-validation for truncated Gaussian kernel density estimation
//'
//' Given a data matrix, compute the log density using leave-one-out
//' cross validation, taking in turn an observation as location vector
//' and computing the density of the resulting mixture.
//' @param x \code{n} by \code{d} matrix of observations
//' @param Sigma smoothing positive-definite matrix
//' @param beta vector of constraints for the half-space
//' @return a vector of values for the weighted leave-one-out likelihood
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::colvec tnorm_loo(arma::mat x, arma::mat Omega, arma::vec beta) {
    arma::mat cholinv = chol(inv_sympd(Omega));
    arma::uword d = x.n_cols;
    arma::uword n = x.n_rows;
    if(n < 2){
       Rcpp::stop("Invalid \"x\" argument: must have at least two rows.");
    }
    if(Omega.n_cols != Omega.n_rows){
     Rcpp::stop("Invalid covariance matrix.");
    }
    if(Omega.n_cols != d){
     Rcpp::stop("Invalid dimension for covariance matrix.");
    }
    if(beta.n_elem != d){
      Rcpp::stop("Invalid \"beta\" vector.");
    }
    arma::uword k = 0;
    // arma::mat quad = arma::mat(n, n);
    arma::colvec quadv(n*(n-1)/2);
    arma::colvec objective(n);
    for(arma::uword i = 0; i < n-1; i++){
       for(arma::uword j = i+1; j < n; j++){ // assumes that n>1
          quadv(k) = sum(pow(cholinv * (x.row(i) - x.row(j)).t(), 2));
          k++;
       }
    }
    // double objective = 0;
    arma::uword im = 0;
    arma::uword jm = 0;
    arma::uword ind = 0;
    //arma::colvec objective(n);
    arma::colvec terms(n-1);
    double maxterms = 0;
    double betaxi = 0;
    arma::vec sdTr = sqrt(trans(beta) * Omega * beta);
    double cst = -0.5*d*log(2*arma::datum::pi) + sum(log(cholinv.diag())) - log(n-1.0);
    for(arma::uword i = 0; i < n; i++){
       betaxi = sum(x.row(i) * beta);
       if(betaxi > 0){
       terms.zeros();
       k = 0;
       for(arma::uword j = 0; j < n; j++){
          if(i != j){
             im = std::min(i,j);
             jm = std::max(i,j);
             ind = im*n-im*(im+1)/2+(jm-im)-1;
             terms(k) = -0.5*quadv(ind);
             k++;
          }
       }
       maxterms = max(terms);
       objective(i) += maxterms + log(sum(exp(terms - maxterms))) + cst - log(arma::normcdf(0.0, -betaxi, sdTr(0)));
       } else{
        objective(i) = -arma::datum::inf;
       }
    }
    return objective;
 }


//' Likelihood cross-validation for MIG density estimation
//'
//' Likelihood cross-validation criterion function.
//' @param x \code{n} by \code{d} matrix of observations
//' @inheritParams dmig
//' @return LCV criterion value
//' @keywords internal
//' @export
// [[Rcpp::export]]
double mig_lcv(arma::mat x, arma::mat Omega, arma::vec beta){
   double objective = mean(mig_loo(x, Omega, beta));
   return objective;
}

//' Likelihood cross-validation for Gaussian kernel density estimation
//'
//' Likelihood cross-validation criterion function.
//' @param x \code{n} by \code{d} matrix of observations
//' @param Sigma smoothing positive-definite matrix
//' @param logweights log vector of weights
//' @return LCV criterion value
//' @export
//' @keywords internal
// [[Rcpp::export]]
double gauss_lcv(arma::mat x, arma::mat Sigma, arma::vec logweights){
   double objective = mean(gauss_loo(x, Sigma, logweights));
   return objective;
}

//' Likelihood cross-validation for truncated normal kernel density estimation
//'
//' Likelihood cross-validation criterion function.
//' @param x \code{n} by \code{d} matrix of observations
//' @param Omega smoothing positive-definite matrix
//' @param beta vector of constraints for the half-space
//' @return LCV criterion value
//' @export
//' @keywords internal
// [[Rcpp::export]]
double tnorm_lcv(arma::mat x, arma::mat Omega, arma::vec beta){
   double objective = mean(tnorm_loo(x, Omega, beta));
   return objective;
}

//' Least squares cross-validation for Gaussian kernel density estimation
//'
//' Least squares cross-validation for weighted Gaussian samples.
//' @param xsamp \code{n} by \code{d} random sample for Monte Carlo estimation of bias
//' @param dxsamp \code{n} vector of density for the points from \code{xsamp}
//' @param mckern logical; if \code{TRUE}, uses the kernel as sampler for Monte Carlo estimation
//' @inheritParams gauss_lcv
//' @return least square criterion value
//' @export
//' @keywords internal
// [[Rcpp::export]]
double gauss_lscv(arma::mat x, arma::mat Sigma, arma::vec logweights,
                  arma::mat xsamp, arma::vec dxsamp, bool mckern = true){
   double objective = 2*mean(exp(gauss_loo(x, Sigma, logweights)));
   if(mckern){
      objective -= mean(gauss_kdens_arma(x, x, Sigma, logweights, false));
   } else{
      arma::vec f_x = gauss_kdens_arma(x, xsamp, Sigma, logweights, false);
      objective -= mean(f_x%f_x/dxsamp);
   }
   return objective;
}

//' Least squares cross-validation for truncated Gaussian kernel density estimation
//'
//' Least squares cross-validation for truncated Gaussian samples.
//' @inheritParams gauss_lscv
//' @inheritParams tnorm_lcv
//' @return least square criterion value
//' @export
//' @keywords internal
// [[Rcpp::export]]
double tnorm_lscv(arma::mat x, arma::mat Omega, arma::vec beta,
                  arma::mat xsamp, arma::vec dxsamp, bool mckern = true){
   double objective = 2*mean(exp(tnorm_loo(x, Omega, beta)));
   if(mckern){
      objective -= mean(tnorm_kdens_arma(x, x, Omega, beta, false));
   } else{
      arma::vec f_x = tnorm_kdens_arma(x, xsamp, Omega, beta, false);
      objective -= mean(f_x%f_x/dxsamp);
   }
   return objective;
}

//' Robust likelihood cross-validation for Gaussian kernel density estimation
//'
//' Robust likelihood cross-validation criterion function of Wu.
//' @inheritParams gauss_lscv
//' @inheritParams gauss_lcv
//' @param an threshold for linear approximation
//' @return RLCV criterion value
//' @export
//' @keywords internal
// [[Rcpp::export]]
double gauss_rlcv(arma::mat x, arma::mat Sigma, arma::vec logweights, double an,
                  arma::mat xsamp, arma::vec dxsamp, bool mckern = true){
  arma::uword n = x.n_rows;
  arma::vec logf_ix = gauss_loo(x, Sigma, logweights);
  double logan = log(an);
  double objective = 0;
  int k = 0;
  for(arma::uword i = 0; i < n; i++){
     if(logf_ix(i) >= logan){
       objective += logf_ix(i);
     } else{
        k++;
       objective += exp(logf_ix(i) - logan);
     }
  }
  objective += (logan - 1.0)*k;
  k = 0;
  double bias = 0;
  if(mckern){
     arma::vec f_x = gauss_kdens_arma(x, x, Sigma, logweights, false);
     for(arma::uword i = 0; i < n; i++){
        if(f_x(i) >= an){
           k++;
        } else{
           bias += f_x(i);
        }
     }
     bias = (bias*0.5/an + (double) k) / ((double) n);
  } else{
     for(arma::uword i = 0; i < xsamp.n_rows; i++){
        arma::vec f_x = gauss_kdens_arma(x, xsamp, Sigma, logweights, false);
        if(f_x(i) >= an){
           bias += f_x(i)/dxsamp(i);
        } else{
           bias += f_x(i)*f_x(i)*0.5/(dxsamp(i)*an);
        }
     }
     bias = bias / ((double) xsamp.n_rows);
  }
  objective = objective / ((double) n) - bias;
  return objective;
}


//' Likelihood cross-validation for truncated normal kernel density estimation
//'
//' Robust likelihood cross-validation criterion function.
//' @inheritParams gauss_lscv
//' @inheritParams tnorm_lcv
//' @param an threshold for linear approximation
//' @return RLCV criterion value
//' @export
//' @keywords internal
// [[Rcpp::export]]
double tnorm_rlcv(arma::mat x, arma::mat Omega, arma::vec beta, double an,
                  arma::mat xsamp, arma::vec dxsamp, bool mckern = true){
   arma::uword n = x.n_rows;
   arma::vec logf_ix = tnorm_loo(x, Omega, beta);
   double logan = log(an);
   double objective = 0;
   for(arma::uword i = 0; i < n; i++){
      if(logf_ix(i) >= logan){
         objective += logf_ix(i);
      } else{
         objective += logan - 1 + exp(logf_ix(i) - logan);
      }
   }
   double bias = 0;
   if(mckern){
      arma::vec f_x = tnorm_kdens_arma(x, x, Omega, beta, false);
      for(arma::uword i = 0; i < n; i++){
         if(f_x(i) >= an){
            bias += 1;
         } else{
            bias += f_x(i)*0.5/an;
         }
      }
      bias = bias / ((double) n);
   } else{
      if(xsamp.n_rows != dxsamp.n_rows){
         Rcpp::stop("Invalid Monte Carlo sample or density.");
      }
      arma::vec f_x = tnorm_kdens_arma(x, xsamp, Omega, beta, false);
      for(arma::uword i = 0; i < xsamp.n_rows; i++){
         if(f_x(i) >= an){
            bias += f_x(i)/dxsamp(i);
         } else{
            bias += f_x(i)*f_x(i)*0.5/(dxsamp(i)*an);
         }
      }
      bias = bias / ((double) xsamp.n_rows);
   }
   objective = objective / n - bias;
   return objective;
}



//' Likelihood cross-validation for truncated normal kernel density estimation
 //'
 //' Robust likelihood cross-validation criterion function.
 //' @inheritParams gauss_lscv
 //' @inheritParams tnorm_lcv
 //' @param an threshold for linear approximation
 //' @return RLCV criterion value
 //' @export
 //' @keywords internal
 // [[Rcpp::export]]
double mig2_rlcv(arma::mat x, arma::mat Omega, arma::vec beta, double an,
                   arma::mat xsamp, arma::vec dxsamp, bool mckern = true){
    arma::uword n = x.n_rows;
    arma::vec logf_ix = mig_loo(x, Omega, beta);
    double logan = log(an);
    double objective = 0;
    for(arma::uword i = 0; i < n; i++){
       if(logf_ix(i) >= logan){
          objective += logf_ix(i);
       } else{
          objective += logan - 1 + exp(logf_ix(i) - logan);
       }
    }
    double bias = 0;
    if(mckern){
       arma::vec f_x = mig_kdens_arma(x, x, Omega, beta, false);
       for(arma::uword i = 0; i < n; i++){
          if(f_x(i) >= an){
             bias += 1;
          } else{
             bias += f_x(i)*0.5/an;
          }
       }
       bias = bias / ((double) n);
    } else{
       if(xsamp.n_rows != dxsamp.n_rows){
          Rcpp::stop("Invalid Monte Carlo sample or density.");
       }
       arma::vec f_x = mig_kdens_arma(x, xsamp, Omega, beta, false);
       for(arma::uword i = 0; i < xsamp.n_rows; i++){
          if(f_x(i) >= an){
             bias += f_x(i)/dxsamp(i);
          } else{
             bias += f_x(i)*f_x(i)*0.5/(dxsamp(i)*an);
          }
       }
       bias = bias / ((double) xsamp.n_rows);
    }
    objective = objective / n - bias;
    return objective;
 }


//' Least squares cross-validation for MIG kernel density estimation
 //'
 //' Least squares cross-validation for multivariate inverse Gaussian samples.
 //' @inheritParams gauss_lscv
 //' @inheritParams tnorm_lcv
 //' @return least square criterion value
 //' @export
 //' @keywords internal
 // [[Rcpp::export]]
double mig2_lscv(arma::mat x, arma::mat Omega, arma::vec beta,
                   arma::mat xsamp, arma::vec dxsamp, bool mckern = true){
    double objective = 2*mean(exp(mig_loo(x, Omega, beta)));
    if(mckern){
       objective -= mean(mig_kdens_arma(x, x, Omega, beta, false));
    } else{
       arma::vec f_x = mig_kdens_arma(x, xsamp, Omega, beta, false);
       objective -= mean(f_x%f_x/dxsamp);
    }
    return objective;
 }


//' Likelihood cross-validation for MIG kernel density estimation
 //'
 //' Likelihood cross-validation criterion function.
 //' @param x \code{n} by \code{d} matrix of observations
 //' @param Omega smoothing positive-definite matrix
 //' @param beta vector of constraints for the half-space
 //' @return LCV criterion value
 //' @export
 //' @keywords internal
 // [[Rcpp::export]]
 double mig2_lcv(arma::mat x, arma::mat Omega, arma::vec beta){
    double objective = mean(mig_loo(x, Omega, beta));
    return objective;
 }
