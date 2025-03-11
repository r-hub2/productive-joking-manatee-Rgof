#include <Rcpp.h>
using namespace Rcpp;
#include "teststatistics.h"
//' run gof tests for continuous data
//' 
//' @param x A numeric vector of data
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param qnull R function (quantiles under null hypothesis)
//' @param w function to calculate weights, returns -99 if no weights
//' @param phat  function to set or estimate parameters of pnull 
//' @param TS function that calculates test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B (=5000) Number of simulation runs 
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function rnull, 
        Rcpp::Function qnull,
        Rcpp::Function w,
        Rcpp::Function phat, 
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        int B=5000) {

  int n=x.size(), i, j;
  NumericVector Fx(n), wx(n);
  NumericVector p=phat(x), psim=phat(x);
/* Check the number of arguments of various functions and run them accordingly */  
  if(std::abs(p(0)+99)<0.01) Fx=pnull(x);
  else Fx=pnull(x, p);
/* Find test statistics for the data */  
  NumericVector TS_data=ts_C(typeTS, x, TS, pnull, p, w, qnull, TSextra);
    int const nummethods=TS_data.size();
  Rcpp::CharacterVector methods=TS_data.names();
  NumericVector xsim(n), TS_sim(nummethods), pvals(nummethods);
  NumericMatrix out(2, nummethods);
  colnames(out) = methods;
  for(i=0;i<B;++i) {
    if(std::abs(p(0)+99)<0.01) xsim=rnull();
    else xsim=rnull(p);
    psim=phat(xsim);
    NumericVector TS_sim=ts_C(typeTS, xsim, TS, pnull, psim, w, qnull, TSextra);
    for(j=0;j<nummethods;++j) {
      if(TS_data(j)<TS_sim(j)) pvals(j)=pvals(j)+1;
    }
  }
  for(j=0;j<nummethods;++j) {
      out(0, j)=TS_data(j);
      out(1, j)=pvals(j)/B;
  }
  return out;
}
