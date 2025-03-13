#ifndef CALCTS_H
#define CALCTS_H

#include <Rcpp.h>
Rcpp::NumericVector calcTS(
     Rcpp::List dta,
     Rcpp::Function pnull,
     Rcpp::NumericVector param,
     Rcpp::Function TS,      
     int typeTS,
     Rcpp::List TSextra
        );

#endif

