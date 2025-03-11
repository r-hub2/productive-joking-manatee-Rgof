#ifndef TS_C_H
#define TS_C_H

#include <Rcpp.h>
Rcpp::NumericVector ts_C(
       int typeTS,
       Rcpp::NumericVector x,
       Rcpp::Function TS, 
       Rcpp::Function pnull, 
       Rcpp::NumericVector param, 
       Rcpp::Function w, 
       Rcpp::Function qnull,
       Rcpp::List TSextra
        );

#endif

#ifndef TS_D_H
#define TS_D_H

#include <Rcpp.h>
Rcpp::NumericVector ts_D(
       int typeTS,
       Rcpp::IntegerVector x,
       Rcpp::Function TS, 
       Rcpp::Function pnull, 
       Rcpp::NumericVector param,        
       Rcpp::NumericVector vals,
       Rcpp::List TSextra
        );

#endif
