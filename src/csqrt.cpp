#include <Rcpp.h>
using namespace Rcpp;

//' Take the square root of a number
//'
//' @param x A number.
//' @return sqrt(x)
//' @export
// [[Rcpp::export]]
NumericVector csqrt(NumericVector x) {
  return sqrt(x);
}
