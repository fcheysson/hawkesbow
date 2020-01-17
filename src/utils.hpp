#pragma once
#include <RcppArmadillo.h>

#if !defined(MYLIB_UTILS)
#define MYLIB_UTILS

const arma::cx_double i(0.0, 1.0);
const double inv_pi = 1.0 / arma::datum::pi;

arma::vec _sinc( arma::vec x );

// Powers of 10
double quick_pow10(int n);
double quick_negpow10(int n);

// Padé approximants for the exponential integral of imaginary argument
// cf. https://en.wikipedia.org/wiki/Exponential_integral
// E1(x) = i * (- 1/2 * pi + Si(x)) - Ci(x),    x > 0
double padef( double x );
double padeg( double x );
double Ci( double x );
double Si( double x );

// arma::uword modulus( arma::sword a, arma::sword b );
// arma::uword modulus( arma::sword a, arma::uword b );

#endif
