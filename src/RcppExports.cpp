// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mig_kdens_arma
Rcpp::NumericVector mig_kdens_arma(arma::mat x, arma::mat newdata, arma::mat Omega, arma::colvec beta, bool logd);
RcppExport SEXP _mig_mig_kdens_arma(SEXP xSEXP, SEXP newdataSEXP, SEXP OmegaSEXP, SEXP betaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type newdata(newdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(mig_kdens_arma(x, newdata, Omega, beta, logd));
    return rcpp_result_gen;
END_RCPP
}
// tnorm_kdens_arma
Rcpp::NumericVector tnorm_kdens_arma(arma::mat x, arma::mat newdata, arma::mat Omega, arma::colvec beta, bool logd);
RcppExport SEXP _mig_tnorm_kdens_arma(SEXP xSEXP, SEXP newdataSEXP, SEXP OmegaSEXP, SEXP betaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type newdata(newdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(tnorm_kdens_arma(x, newdata, Omega, beta, logd));
    return rcpp_result_gen;
END_RCPP
}
// gauss_kdens_arma
Rcpp::NumericVector gauss_kdens_arma(arma::mat x, arma::mat newdata, arma::mat Sigma, arma::vec logweights, bool logd);
RcppExport SEXP _mig_gauss_kdens_arma(SEXP xSEXP, SEXP newdataSEXP, SEXP SigmaSEXP, SEXP logweightsSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type newdata(newdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logweights(logweightsSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_kdens_arma(x, newdata, Sigma, logweights, logd));
    return rcpp_result_gen;
END_RCPP
}
// mig_loo
arma::colvec mig_loo(arma::mat x, arma::colvec beta, arma::mat Omega);
RcppExport SEXP _mig_mig_loo(SEXP xSEXP, SEXP betaSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    rcpp_result_gen = Rcpp::wrap(mig_loo(x, beta, Omega));
    return rcpp_result_gen;
END_RCPP
}
// mig_rlcv
double mig_rlcv(arma::mat x, arma::colvec beta, arma::mat Omega, double an, arma::mat xsamp, arma::vec dxsamp, bool mckern);
RcppExport SEXP _mig_mig_rlcv(SEXP xSEXP, SEXP betaSEXP, SEXP OmegaSEXP, SEXP anSEXP, SEXP xsampSEXP, SEXP dxsampSEXP, SEXP mckernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xsamp(xsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dxsamp(dxsampSEXP);
    Rcpp::traits::input_parameter< bool >::type mckern(mckernSEXP);
    rcpp_result_gen = Rcpp::wrap(mig_rlcv(x, beta, Omega, an, xsamp, dxsamp, mckern));
    return rcpp_result_gen;
END_RCPP
}
// mig_lscv
double mig_lscv(arma::mat x, arma::colvec beta, arma::mat Omega, arma::mat xsamp, arma::vec dxsamp, bool mckern);
RcppExport SEXP _mig_mig_lscv(SEXP xSEXP, SEXP betaSEXP, SEXP OmegaSEXP, SEXP xsampSEXP, SEXP dxsampSEXP, SEXP mckernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xsamp(xsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dxsamp(dxsampSEXP);
    Rcpp::traits::input_parameter< bool >::type mckern(mckernSEXP);
    rcpp_result_gen = Rcpp::wrap(mig_lscv(x, beta, Omega, xsamp, dxsamp, mckern));
    return rcpp_result_gen;
END_RCPP
}
// gauss_loo
arma::colvec gauss_loo(arma::mat x, arma::mat Sigma, arma::vec logweights);
RcppExport SEXP _mig_gauss_loo(SEXP xSEXP, SEXP SigmaSEXP, SEXP logweightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logweights(logweightsSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_loo(x, Sigma, logweights));
    return rcpp_result_gen;
END_RCPP
}
// tnorm_loo
arma::colvec tnorm_loo(arma::mat x, arma::mat Omega, arma::vec beta);
RcppExport SEXP _mig_tnorm_loo(SEXP xSEXP, SEXP OmegaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(tnorm_loo(x, Omega, beta));
    return rcpp_result_gen;
END_RCPP
}
// gauss_lcv
double gauss_lcv(arma::mat x, arma::mat Sigma, arma::vec logweights);
RcppExport SEXP _mig_gauss_lcv(SEXP xSEXP, SEXP SigmaSEXP, SEXP logweightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logweights(logweightsSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_lcv(x, Sigma, logweights));
    return rcpp_result_gen;
END_RCPP
}
// tnorm_lcv
double tnorm_lcv(arma::mat x, arma::mat Omega, arma::vec beta);
RcppExport SEXP _mig_tnorm_lcv(SEXP xSEXP, SEXP OmegaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(tnorm_lcv(x, Omega, beta));
    return rcpp_result_gen;
END_RCPP
}
// gauss_lscv
double gauss_lscv(arma::mat x, arma::mat Sigma, arma::vec logweights, arma::mat xsamp, arma::vec dxsamp, bool mckern);
RcppExport SEXP _mig_gauss_lscv(SEXP xSEXP, SEXP SigmaSEXP, SEXP logweightsSEXP, SEXP xsampSEXP, SEXP dxsampSEXP, SEXP mckernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logweights(logweightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xsamp(xsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dxsamp(dxsampSEXP);
    Rcpp::traits::input_parameter< bool >::type mckern(mckernSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_lscv(x, Sigma, logweights, xsamp, dxsamp, mckern));
    return rcpp_result_gen;
END_RCPP
}
// tnorm_lscv
double tnorm_lscv(arma::mat x, arma::mat Omega, arma::vec beta, arma::mat xsamp, arma::vec dxsamp, bool mckern);
RcppExport SEXP _mig_tnorm_lscv(SEXP xSEXP, SEXP OmegaSEXP, SEXP betaSEXP, SEXP xsampSEXP, SEXP dxsampSEXP, SEXP mckernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xsamp(xsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dxsamp(dxsampSEXP);
    Rcpp::traits::input_parameter< bool >::type mckern(mckernSEXP);
    rcpp_result_gen = Rcpp::wrap(tnorm_lscv(x, Omega, beta, xsamp, dxsamp, mckern));
    return rcpp_result_gen;
END_RCPP
}
// gauss_rlcv
double gauss_rlcv(arma::mat x, arma::mat Sigma, arma::vec logweights, double an, arma::mat xsamp, arma::vec dxsamp, bool mckern);
RcppExport SEXP _mig_gauss_rlcv(SEXP xSEXP, SEXP SigmaSEXP, SEXP logweightsSEXP, SEXP anSEXP, SEXP xsampSEXP, SEXP dxsampSEXP, SEXP mckernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logweights(logweightsSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xsamp(xsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dxsamp(dxsampSEXP);
    Rcpp::traits::input_parameter< bool >::type mckern(mckernSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_rlcv(x, Sigma, logweights, an, xsamp, dxsamp, mckern));
    return rcpp_result_gen;
END_RCPP
}
// tnorm_rlcv
double tnorm_rlcv(arma::mat x, arma::mat Omega, arma::vec beta, double an, arma::mat xsamp, arma::vec dxsamp, bool mckern);
RcppExport SEXP _mig_tnorm_rlcv(SEXP xSEXP, SEXP OmegaSEXP, SEXP betaSEXP, SEXP anSEXP, SEXP xsampSEXP, SEXP dxsampSEXP, SEXP mckernSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xsamp(xsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dxsamp(dxsampSEXP);
    Rcpp::traits::input_parameter< bool >::type mckern(mckernSEXP);
    rcpp_result_gen = Rcpp::wrap(tnorm_rlcv(x, Omega, beta, an, xsamp, dxsamp, mckern));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mig_mig_kdens_arma", (DL_FUNC) &_mig_mig_kdens_arma, 5},
    {"_mig_tnorm_kdens_arma", (DL_FUNC) &_mig_tnorm_kdens_arma, 5},
    {"_mig_gauss_kdens_arma", (DL_FUNC) &_mig_gauss_kdens_arma, 5},
    {"_mig_mig_loo", (DL_FUNC) &_mig_mig_loo, 3},
    {"_mig_mig_rlcv", (DL_FUNC) &_mig_mig_rlcv, 7},
    {"_mig_mig_lscv", (DL_FUNC) &_mig_mig_lscv, 6},
    {"_mig_gauss_loo", (DL_FUNC) &_mig_gauss_loo, 3},
    {"_mig_tnorm_loo", (DL_FUNC) &_mig_tnorm_loo, 3},
    {"_mig_gauss_lcv", (DL_FUNC) &_mig_gauss_lcv, 3},
    {"_mig_tnorm_lcv", (DL_FUNC) &_mig_tnorm_lcv, 3},
    {"_mig_gauss_lscv", (DL_FUNC) &_mig_gauss_lscv, 6},
    {"_mig_tnorm_lscv", (DL_FUNC) &_mig_tnorm_lscv, 6},
    {"_mig_gauss_rlcv", (DL_FUNC) &_mig_gauss_rlcv, 7},
    {"_mig_tnorm_rlcv", (DL_FUNC) &_mig_tnorm_rlcv, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_mig(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
