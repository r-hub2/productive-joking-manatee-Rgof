#include <Rcpp.h>
using namespace Rcpp;
//' This function calculates the test statistics for continuous data
//' @param  typeTS format of TS
//' @param  x continuous data set
//' @param  TS routine
//' @param  pnull  cdf under the null hypothesis
//' @param  param estimated parameters
//' @param  w routine to calculate weights
//' @param  qnull quantile function
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector ts_C(
       int typeTS,
       Rcpp::NumericVector x,     
       Rcpp::Function TS, 
       Rcpp::Function pnull,
       Rcpp::NumericVector param,
       Rcpp::Function w,
       Rcpp::Function qnull,
       Rcpp::List TSextra) {
  NumericVector TS_data, wx;
  if(typeTS==1) TS_data=TS(x, pnull, param, qnull);
  if(typeTS==2) {
    if(std::abs(param(0)+99)<0.01) wx=w(x);
    else wx=w(x, param);
    TS_data=TS(x, pnull, param, wx);
  }  
  if(typeTS==3) TS_data=TS(x, pnull, param);
  if(typeTS==4) TS_data=TS(x, pnull, param, TSextra);
  return  TS_data;
}

//' This function calculates the test statistics for discrete data
//' @param  typeTS format of TS
//' @param  x discrete data set (counts)
//' @param  TS routine 
//' @param  pnull  cdf under the null hypothesis
//' @param  param estimated parameters
//' @param  vals values of discrete RV
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector ts_D(
     int typeTS,
     Rcpp::IntegerVector x,     
     Rcpp::Function TS, 
     Rcpp::Function pnull,
     Rcpp::NumericVector param,
     Rcpp::NumericVector vals, 
     Rcpp::List TSextra) {
   NumericVector TS_data;
   if(typeTS<=6) TS_data=TS(x, pnull, param, vals); 
   if(typeTS==7) TS_data=TS(x, pnull, param, vals, TSextra);
   return  TS_data;
 }
 
 
