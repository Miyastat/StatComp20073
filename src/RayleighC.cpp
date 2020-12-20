#include <Rcpp.h>
using namespace Rcpp;

//' @title A Rayleigh sampler using Rcpp
//' @description A Rayleigh sampler using Rcpp
//' @param sigma parameter of rayleigh density
//' @param n size of Rayleigh random number
//' @return a Rayleigh random sample
//' @export
// [[Rcpp::export]]
NumericVector RayleighC(int n, double sigma) {
  NumericVector y(n);
  int k = 0,j = 0;
  double u=0, x=0;
  while (k < n) {
    u = runif(1)[0];
    j = j + 1;
    x = rchisq(1,4)[0];
    if (exp(-pow(x,2)/(2*pow(sigma,2))+x/2-pow(sigma,2)/8) > u) {
      k++;
      y[k] = x;
    }
  }
  return y;
}
