#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Simple matrix multiplication
//'
//' @param A Matrix
//' @param B Matrix
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' A <- matrix(1:9, 3, 3)
//' B <- matrix(11:19, 3, 3)
//' matrix_mult_cpp(A, B)
// [[Rcpp::export]]
arma::mat matrix_mult_cpp(arma::mat A, arma::mat B) {
  return A * B;
}
