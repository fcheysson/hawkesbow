#include "utils.hpp"

// // Sinus cardinal: sinc(x) = sin(x) / x
// //' @export
// // [[Rcpp::export]]
// arma::vec sinc( arma::vec x ) {
//     arma::vec y(x.n_elem);
//
//     // Iterators
//     arma::vec::iterator it_x = x.begin();
//     arma::vec::iterator it_x_end = x.end();
//     arma::vec::iterator it_y = y.begin();
//     Rcpp::Rcout << "just checking" << std::endl;
//
//     // Loop on x
//     for (; it_x != it_x_end; ++it_x, ++it_y) {
//         if (*it_x == 0.0)
//             *it_y = 1.0;
//         else
//             *it_y = sin(*it_x) / *it_x;
//     }
//     return y;
// };

arma::uword modulus( arma::sword a, arma::sword b ) {
  return (arma::uword)( (a%b+b)%b );
}

arma::uword modulus( arma::sword a, arma::uword b ) {
  arma::sword c = (arma::sword)b;
  return (arma::uword)( (a%c+c)%c );
}
