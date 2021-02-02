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
double integral_simpson(double(*f)(double x, double nu, double r), double a, double b, int n, double nu, double r);
double quadrant_real(double x, double nu, double r);
double quadrant_imag(double x, double nu, double r);
arma::cx_double contour_quadrant(double nu, double r);

//' Incomplete gamma function of imaginary argument
//'
//' Calculates the value of
//' \deqn{\Gamma_1(\nu, r) = \int_r^\infty y^{\nu-1} e^{-iy} \mathrm{d}y}
//' for \eqn{0 < \nu < 1} through the following relation (obtained by contour integration)
//' \deqn{\int_r^\infty y^{\nu-1} e^{-iy} \mathrm{d}y =
//' e^{-i\frac{\pi}{2}\nu} \int_r^\infty x^{\nu-1} e^{-x} \mathrm{d}x -
//' e^{-i\frac{\pi}{2}(\nu-1)} r^\nu \int_0^{\pi/2}
//' e^{i\theta \nu}e^{-re^{i\theta}}\mathrm{d}\theta.}
//' The first integral is calculated using function "tgamma" from the library
//' "boost::math", while the second is approximated via Simpson's rule.
//'
//' @param nu A number between 0 and 1 (strictly)
//' @param r A non-negative number
//'
//' @return The incomplete gamma function of imaginary argument (see Details)
//' @export
//'
//' @examples
//' inc_gamma_imag(0.5, 1.0)
// [[Rcpp::export]]
arma::cx_double inc_gamma_imag(double nu, double r);

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

//' Exponential integral of imaginary argument
//'
//' Calculates the value of
//' \deqn{E_1(ix) = \int_1^\infty \frac{e^{-ixt}}{t} \mathrm{d}t}
//' using its relation to the trigonometric integrals
//' (cf. \url{https://en.wikipedia.org/wiki/Exponential_integral#Exponential_integral_of_imaginary_argument}):
//' \deqn{E_1(ix) = i \left[ -\frac{1}{2} \pi + Si(x) \right] - Ci(x)}
//' and their Pad\'e approximants
//' (cf. \url{https://en.wikipedia.org/wiki/Trigonometric_integral#Efficient_evaluation})
//'
//' @param x A non-negative number
//'
//' @return The exponential integral of argument ix
//' @export
//'
//' @examples
//' E1_imaginary(1.0)
// [[Rcpp::export]]
arma::cx_double E1_imaginary( double x );

#endif