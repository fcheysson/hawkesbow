#pragma once
#include <RcppArmadillo.h>

#if !defined(MYLIB_UTILS)
#define MYLIB_UTILS

const arma::cx_double i(0.0, 1.0);
const double inv_pi = 1.0 / arma::datum::pi;

// arma::vec sinc( arma::vec x );
// sinc is now a function in armadillo

arma::uword modulus( arma::sword a, arma::sword b );
arma::uword modulus( arma::sword a, arma::uword b );

#endif
