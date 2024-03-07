// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' Likelihood cross-validation for kernel density estimation with MIG
//'
//' Given a data matrix over a half-space defined by \code{beta},
//' compute the log density using leave-one-out cross validation,
//' taking in turn an observation as location vector and computing the
//' density of the resulting mixture.
//' @inheritParams dmig
//' @return the value of the likelihood cross-validation criterion
//' @export
// [[Rcpp::export]]
double mig_lcv(arma::mat x, arma::rowvec beta, arma::mat Omega) {
  arma::mat cholinv = chol(Omega.i());
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
  double maxterms = 0;
  // double term = 0;
  if(beta.n_cols != d){
   Rcpp::stop("Invalid arguments: \"beta\" must be the same length as the number of columns of \"x\".");
  }
  arma::colvec cp = arma::colvec(n);
  k = 0;
  for(arma::uword i = 0; i < n-1; i++){
   for(arma::uword j = i+1; j < n; j++){ // assumes that n>1
    quadv(k) = sum(pow(cholinv * (x.row(i) - x.row(j)).t(), 2));
      k++;
    }
   cp(i) = sum(beta % x.row(i));
  }
  cp(n-1) =  sum(beta % x.row(n-1));
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
     objective += maxterms + log(sum(exp(terms - maxterms)));
  }
  objective = objective/n -0.5*d*log(2*M_PI) + sum(log(cp))/n + sum(log(cholinv.diag())) - log(n-1.0);
 return objective;
    // Rcpp::List::create(
    // Rcpp::Named("quadv") = quadv,
    // Rcpp::Named("quad") = quad,
    // Rcpp::Named("total") = term,
    // Rcpp::Named("objective") = objective);
}


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
 arma::mat cholinv = chol(Omega.i());
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
    arma::mat cholinv = chol(Omega.i());
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
    return logdens;
 }
