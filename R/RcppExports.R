# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' sort vector y by values in vector x
#' 
#' @param y numeric vector
#' @param x numeric vector
#' @keywords internal
#' @return numeric vector
Cpporder <- function(y, x) {
    .Call(`_Rgof_Cpporder`, y, x)
}

#' Find test statistics for continuous data
#' 
#' @param x A numeric vector.
#' @param pnull cdf.
#' @param param parameters for pnull  in case of parameter estimation.
#' @param qnull An R function, the quantile function under the null hypothesis.
#' @return A numeric vector with test statistics
TS_cont <- function(x, pnull, param, qnull) {
    .Call(`_Rgof_TS_cont`, x, pnull, param, qnull)
}

#' Find test statistics for discrete data
#' 
#' @param x An integer vector.
#' @param pnull cdf.
#' @param param parameters for pnull in case of parameter estimation.
#' @param vals A numeric vector with the values of the discrete rv.
#' @return A vector with test statistics
TS_disc <- function(x, pnull, param, vals) {
    .Call(`_Rgof_TS_disc`, x, pnull, param, vals)
}

#' Find test statistics for continuous data with weights
#' 
#' @param x A numeric vector.
#' @param pnull cdf.
#' @param param parameters for pnull in case of parameter estimation.
#' @param w numeric vector of weights
#' @keywords internal
#' @return A numeric vector with test statistics
TSw_cont <- function(x, pnull, param, w) {
    .Call(`_Rgof_TSw_cont`, x, pnull, param, w)
}

#' Find test statistics for discrete data
#' 
#' @param x An integer vector.
#' @param pnull cdf.
#' @param param parameters for pnull in case of parameter estimation.
#' @param vals A numeric vector with the values of the discrete rv.
#' @param w weights 
#' @keywords internal
#' @return A vector with test statistics
TSw_disc <- function(x, pnull, param, vals, w) {
    .Call(`_Rgof_TSw_disc`, x, pnull, param, vals, w)
}

#' count events in bins. Useful for power calculations. Replaces hist command from R.
#' 
#' @param x numeric vector
#' @param bins numeric vector
#' @keywords internal
#' @return Integer vector of counts
bincounter <- function(x, bins) {
    .Call(`_Rgof_bincounter`, x, bins)
}

#' This function calculates the test statistics for  data
#' @param  dta data set as a list
#' @param  TS routine
#' @param  typeTS format of TS
#' @param  TSextra list passed to TS function
#' @keywords internal
#' @return A vector of numbers
calcTS <- function(dta, TS, typeTS, TSextra) {
    .Call(`_Rgof_calcTS`, dta, TS, typeTS, TSextra)
}

#' helper functions to do p value adjustment
#' @param dta data set
#' @param rnull R function (generate data under null hypothesis)
#' @param vals values of discrete random variable
#' @param TS function to calculate test statistics
#' @param typeTS integer indicating type of test statistic
#' @param TSextra list to pass to TS
#' @param B  =1000 Number of simulation runs
#' @keywords internal
#' @return A matrix of powers
gof_adj_C1 <- function(dta, rnull, vals, TS, typeTS, TSextra, B = 1000L) {
    .Call(`_Rgof_gof_adj_C1`, dta, rnull, vals, TS, typeTS, TSextra, B)
}

#' helper functions to do p value adjustment
#' @param dta data set
#' @param rnull R function (generate data under null hypothesis)
#' @param vals values of discrete random variable
#' @param TS function to calculate test statistics
#' @param typeTS integer indicating type of test statistic
#' @param TSextra list to pass to TS
#' @param A matrix of test statistic values 
#' @param B  =1000 Number of simulation runs
#' @keywords internal
#' @return A matrix of powers
gof_adj_C2 <- function(dta, rnull, vals, TS, typeTS, TSextra, A, B = 1000L) {
    .Call(`_Rgof_gof_adj_C2`, dta, rnull, vals, TS, typeTS, TSextra, A, B)
}

#' find power of gof tests for continuous data
#' 
#' @param rnull R function (generate data under null hypothesis)
#' @param vals values of discrete random variable
#' @param ralt  R function to generate data under alternative
#' @param param_alt parameters of ralt
#' @param TS function to calculate test statistics
#' @param typeTS integer indicating type of test statistic
#' @param TSextra list to pass to TS
#' @param B  =1000 Number of simulation runs
#' @keywords internal
#' @return A matrix of powers
gof_power_C <- function(rnull, vals, ralt, param_alt, TS, typeTS, TSextra, B = 1000L) {
    .Call(`_Rgof_gof_power_C`, rnull, vals, ralt, param_alt, TS, typeTS, TSextra, B)
}

#' run gof tests for continuous data
#' 
#' @param dta A numeric vector of data
#' @param rnull R function (generate data under null hypothesis)
#' @param TS function that calculates test statistics
#' @param typeTS integer indicating type of test statistic
#' @param TSextra list to pass to TS
#' @param B (=5000) Number of simulation runs 
#' @keywords internal
#' @return A matrix of numbers
gof_test_C <- function(dta, rnull, TS, typeTS, TSextra, B = 5000L) {
    .Call(`_Rgof_gof_test_C`, dta, rnull, TS, typeTS, TSextra, B)
}

#' a local function needed for the vignette
#' 
#' @param x An integer vector.
#' @param pnull cdf.
#' @param param parameters for pnull in case of parameter estimation.
#' @param vals A numeric vector with the values of the discrete rv.
#' @return A vector with test statistics
#' @export
newTSdisc <- function(x, pnull, param, vals) {
    .Call(`_Rgof_newTSdisc`, x, pnull, param, vals)
}

#' Find counts or sum of weights in bins. Useful for power calculations. Replaces hist command from R.
#' 
#' @param x numeric vector
#' @param bins numeric vector
#' @param w numeric vector of weights 
#' @keywords internal
#' @return sum of weights in bins
wbincounter <- function(x, bins, w) {
    .Call(`_Rgof_wbincounter`, x, bins, w)
}

