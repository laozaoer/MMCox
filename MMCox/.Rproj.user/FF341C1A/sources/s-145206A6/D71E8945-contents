#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::mat add1(const arma::mat&x, const arma::uword&e) {
  arma::mat y=x;
  y.insert_rows(e-1,2*y.row(e-1));
  return y;
}