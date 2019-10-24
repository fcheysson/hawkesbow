#pragma once
#include <RcppArmadillo.h>

#if !defined(MYLIB_UTILS)
#define MYLIB_UTILS

const arma::cx_double i(0.0, 1.0);
const double inv_pi = 1.0 / arma::datum::pi;

arma::vec _sinc( arma::vec x );

arma::uword modulus( arma::sword a, arma::sword b );
arma::uword modulus( arma::sword a, arma::uword b );

#endif
