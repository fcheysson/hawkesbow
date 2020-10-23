#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo, BH)]]
using namespace Rcpp;

#if !defined(MYLIB_UTILS)
#define MYLIB_UTILS

const arma::cx_double i(0.0, 1.0);
const double inv_pi = 1.0 / arma::datum::pi;
arma::cx_double pow_i(int x);
double pow_m1(int x);

arma::vec _sinc( arma::vec x );

// Contour integration for the incomplete gamma function with imaginary argument
// Can be checked using function gamma from R::base, Gradshteyn & Ryzhik, 2007,
// and Barakat, 1960
double integral_midpoint(double(*f)(double x), double a, double b, int n);
double integral_midpoint(double(*f)(double x, double nu), double a, double b, int n, double nu);
double integral_simpson(double(*f)(double x), double a, double b, int n);
double integral_simpson(double(*f)(double x, double nu), double a, double b, int n, double nu);
double quadrant_real(double x, double nu);
double quadrant_imag(double x, double nu);
arma::cx_double contour_quadrant(double nu);
arma::cx_double inc_gamma_imag(double nu);

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
arma::cx_double E1_imaginary( double x );

#endif
