// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Cpporder
NumericVector Cpporder(NumericVector y, NumericVector x);
RcppExport SEXP _Rgof_Cpporder(SEXP ySEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Cpporder(y, x));
    return rcpp_result_gen;
END_RCPP
}
// TS_cont
Rcpp::NumericVector TS_cont(Rcpp::NumericVector x, Rcpp::NumericVector Fx, Rcpp::NumericVector param, Rcpp::Function qnull);
RcppExport SEXP _Rgof_TS_cont(SEXP xSEXP, SEXP FxSEXP, SEXP paramSEXP, SEXP qnullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fx(FxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type qnull(qnullSEXP);
    rcpp_result_gen = Rcpp::wrap(TS_cont(x, Fx, param, qnull));
    return rcpp_result_gen;
END_RCPP
}
// TS_disc
NumericVector TS_disc(IntegerVector x, NumericVector Fx, NumericVector vals);
RcppExport SEXP _Rgof_TS_disc(SEXP xSEXP, SEXP FxSEXP, SEXP valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Fx(FxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vals(valsSEXP);
    rcpp_result_gen = Rcpp::wrap(TS_disc(x, Fx, vals));
    return rcpp_result_gen;
END_RCPP
}
// TSw_cont
Rcpp::NumericVector TSw_cont(Rcpp::NumericVector x, Rcpp::NumericVector Fx, Rcpp::NumericVector w);
RcppExport SEXP _Rgof_TSw_cont(SEXP xSEXP, SEXP FxSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fx(FxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(TSw_cont(x, Fx, w));
    return rcpp_result_gen;
END_RCPP
}
// TSw_disc
NumericVector TSw_disc(IntegerVector x, NumericVector Fx, NumericVector vals, NumericVector w);
RcppExport SEXP _Rgof_TSw_disc(SEXP xSEXP, SEXP FxSEXP, SEXP valsSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Fx(FxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(TSw_disc(x, Fx, vals, w));
    return rcpp_result_gen;
END_RCPP
}
// bincounter
Rcpp::IntegerVector bincounter(Rcpp::NumericVector x, Rcpp::NumericVector bins);
RcppExport SEXP _Rgof_bincounter(SEXP xSEXP, SEXP binsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bins(binsSEXP);
    rcpp_result_gen = Rcpp::wrap(bincounter(x, bins));
    return rcpp_result_gen;
END_RCPP
}
// gof_cont
Rcpp::NumericMatrix gof_cont(Rcpp::NumericVector x, Rcpp::Function pnull, Rcpp::Function rnull, Rcpp::Function qnull, Rcpp::Function w, Rcpp::Function phat, Rcpp::Function TS, int typeTS, Rcpp::List TSextra, int B);
RcppExport SEXP _Rgof_gof_cont(SEXP xSEXP, SEXP pnullSEXP, SEXP rnullSEXP, SEXP qnullSEXP, SEXP wSEXP, SEXP phatSEXP, SEXP TSSEXP, SEXP typeTSSEXP, SEXP TSextraSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type qnull(qnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type w(wSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type TS(TSSEXP);
    Rcpp::traits::input_parameter< int >::type typeTS(typeTSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type TSextra(TSextraSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(gof_cont(x, pnull, rnull, qnull, w, phat, TS, typeTS, TSextra, B));
    return rcpp_result_gen;
END_RCPP
}
// gof_disc
NumericMatrix gof_disc(Rcpp::IntegerVector x, Rcpp::Function pnull, Rcpp::Function rnull, Rcpp::NumericVector vals, Rcpp::Function phat, Rcpp::Function TS, int typeTS, Rcpp::List TSextra, double rate, int B);
RcppExport SEXP _Rgof_gof_disc(SEXP xSEXP, SEXP pnullSEXP, SEXP rnullSEXP, SEXP valsSEXP, SEXP phatSEXP, SEXP TSSEXP, SEXP typeTSSEXP, SEXP TSextraSEXP, SEXP rateSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type TS(TSSEXP);
    Rcpp::traits::input_parameter< int >::type typeTS(typeTSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type TSextra(TSextraSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(gof_disc(x, pnull, rnull, vals, phat, TS, typeTS, TSextra, rate, B));
    return rcpp_result_gen;
END_RCPP
}
// power_cont
Rcpp::NumericMatrix power_cont(Rcpp::Function pnull, Rcpp::Function rnull, Rcpp::Function qnull, Rcpp::Function ralt, Rcpp::NumericVector param_alt, Rcpp::Function w, Rcpp::Function phat, Rcpp::Function TS, int typeTS, Rcpp::List TSextra, Rcpp::IntegerVector B, const double alpha);
RcppExport SEXP _Rgof_power_cont(SEXP pnullSEXP, SEXP rnullSEXP, SEXP qnullSEXP, SEXP raltSEXP, SEXP param_altSEXP, SEXP wSEXP, SEXP phatSEXP, SEXP TSSEXP, SEXP typeTSSEXP, SEXP TSextraSEXP, SEXP BSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type qnull(qnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type ralt(raltSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param_alt(param_altSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type w(wSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type TS(TSSEXP);
    Rcpp::traits::input_parameter< int >::type typeTS(typeTSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type TSextra(TSextraSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(power_cont(pnull, rnull, qnull, ralt, param_alt, w, phat, TS, typeTS, TSextra, B, alpha));
    return rcpp_result_gen;
END_RCPP
}
// power_disc
Rcpp::NumericMatrix power_disc(Rcpp::Function pnull, Rcpp::Function rnull, Rcpp::NumericVector vals, Rcpp::Function ralt, Rcpp::NumericVector param_alt, Rcpp::Function phat, Rcpp::Function TS, int typeTS, Rcpp::List TSextra, double rate, Rcpp::IntegerVector B, const double alpha);
RcppExport SEXP _Rgof_power_disc(SEXP pnullSEXP, SEXP rnullSEXP, SEXP valsSEXP, SEXP raltSEXP, SEXP param_altSEXP, SEXP phatSEXP, SEXP TSSEXP, SEXP typeTSSEXP, SEXP TSextraSEXP, SEXP rateSEXP, SEXP BSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type ralt(raltSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param_alt(param_altSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type TS(TSSEXP);
    Rcpp::traits::input_parameter< int >::type typeTS(typeTSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type TSextra(TSextraSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(power_disc(pnull, rnull, vals, ralt, param_alt, phat, TS, typeTS, TSextra, rate, B, alpha));
    return rcpp_result_gen;
END_RCPP
}
// wbincounter
Rcpp::NumericVector wbincounter(Rcpp::NumericVector x, Rcpp::NumericVector bins, Rcpp::NumericVector w);
RcppExport SEXP _Rgof_wbincounter(SEXP xSEXP, SEXP binsSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(wbincounter(x, bins, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rgof_Cpporder", (DL_FUNC) &_Rgof_Cpporder, 2},
    {"_Rgof_TS_cont", (DL_FUNC) &_Rgof_TS_cont, 4},
    {"_Rgof_TS_disc", (DL_FUNC) &_Rgof_TS_disc, 3},
    {"_Rgof_TSw_cont", (DL_FUNC) &_Rgof_TSw_cont, 3},
    {"_Rgof_TSw_disc", (DL_FUNC) &_Rgof_TSw_disc, 4},
    {"_Rgof_bincounter", (DL_FUNC) &_Rgof_bincounter, 2},
    {"_Rgof_gof_cont", (DL_FUNC) &_Rgof_gof_cont, 10},
    {"_Rgof_gof_disc", (DL_FUNC) &_Rgof_gof_disc, 10},
    {"_Rgof_power_cont", (DL_FUNC) &_Rgof_power_cont, 12},
    {"_Rgof_power_disc", (DL_FUNC) &_Rgof_power_disc, 12},
    {"_Rgof_wbincounter", (DL_FUNC) &_Rgof_wbincounter, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rgof(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
