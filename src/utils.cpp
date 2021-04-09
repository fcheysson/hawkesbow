#define BOOST_DISABLE_ASSERTS

#include "utils.h"
#include <string>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>

// Powers of i
arma::cx_double pow_i(int x) {
    int rem = x % 4;
    switch(rem) {
    case 1:
        return i;
    case 2:
        return -1;
    case 3:
        return -i;
    case 0:
        return 1;
    }
    return 0;
}

// Powers of -1
double pow_m1(int x) {
    int rem = x % 2;
    switch(rem) {
    case 1:
        return -1;
    case 0:
        return 1;
    }
    return 0;
}

// Sinus cardinal: sinc(x) = sin(x) / x
arma::vec _sinc( arma::vec x ) {
    arma::vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::vec::iterator it_y = y.begin();

    // Loop on x
    for (; it_x != it_x_end; ++it_x, ++it_y) {
        if (*it_x == 0.0)
            *it_y = 1.0;
        else
            *it_y = sin(*it_x) / *it_x;
    }
    return y;
}

arma::cx_vec Etheta_imaginary( double theta, arma::vec x ) {

    arma::cx_vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::cx_vec::iterator it_y = y.begin();

    if (std::fmod(theta, 1.0) == 0) { // Theta is integer

        int itheta = std::floor(theta);

        if (itheta == 1) {

            // Loop on x
            for (; it_x != it_x_end; ++it_x, ++it_y) {
                *it_y = - i * *it_x * exp(i * *it_x) * E1_imaginary(*it_x);
            }

        } else {

            // Term 1 with sum from k = 1 up to theta - 1
            arma::mat xiprod = arma::cumprod(arma::kron(arma::ones<arma::rowvec>(itheta - 1), x), 1);
            arma::vec summands_den = arma::cumprod(theta - arma::regspace(1, itheta - 1));
            arma::cx_vec summands_num = arma::cumprod(-i * arma::ones<arma::cx_vec>(itheta - 1));
            arma::cx_vec term1 = xiprod * (summands_num / summands_den);

            // Get last elements that appear in the sum
            arma::vec last_xiprod = xiprod.col(itheta - 2);
            arma::cx_double last_num = summands_num.back();
            double last_den = summands_den.back();

            // Term 2 with exponential integral
            arma::cx_vec term2(x.n_elem);
            arma::cx_vec::iterator it_term2 = term2.begin();
            arma::vec::iterator it_last = last_xiprod.begin();

            // Loop on xi
            for (; it_x != it_x_end; ++it_x, ++it_term2, ++it_last) {
                *it_term2 = - i * *it_x * last_num * *it_last * exp(i * *it_x) * E1_imaginary(*it_x) / last_den;
            }

            y = term1 + term2;

        }

    } else { // Theta is double, non integer

        int itheta = std::floor(theta);

        if (itheta == 0) {

            // Loop on xi
            for (; it_x != it_x_end; ++it_x, ++it_y) {
                *it_y = - i * *it_x * exp(i * *it_x) * exp((theta-1.0)*log(*it_x)) * inc_gamma_imag(*it_x, 1.0-theta);
            }

        } else {

            // Term 1 with sum from k = 1 up to itheta
            arma::mat xiprod = arma::cumprod(arma::kron(arma::ones<arma::rowvec>(itheta), x), 1);
            arma::vec summands_den = arma::cumprod(theta - arma::regspace(1, itheta));
            arma::cx_vec summands_num = arma::cumprod(-i * arma::ones<arma::cx_vec>(itheta));
            arma::cx_vec term1 = xiprod * (summands_num / summands_den);

            // Get last elements that appear in the denominator
            double last_den = summands_den.back();

            // Term 2 with exponential integral
            arma::cx_vec term2(x.n_elem);
            arma::cx_vec::iterator it_term2 = term2.begin();

            // Loop on xi
            for (; it_x != it_x_end; ++it_x, ++it_term2) {
                *it_term2 = pow_m1(itheta + 1) * pow_i(itheta + 1) * exp(i * *it_x) * exp(theta*log(*it_x)) * inc_gamma_imag(*it_x, 1-theta+itheta) / last_den;
            }

            y = term1 + term2;

        }

    }

    return y;

}

// // Taylor approximants for the incomplete gamma function with imaginary argument
// double Ci( double x, double alpha ) {
//
//     if (x < 0.0)
//         Rcpp::stop("ERROR in Ci: 'x' cannot be negative.");
//
//     if (x == 0.0)
//         return 0.0;
//
//     if (alpha <= 0.0 || alpha >= 1.0) {
//         Rcpp::stop("ERROR in Ci: 'alpha' must be between 0 and 1 strictly.");
//     }
//
//     int n = std::ceil(7.0 + 1.36 * x); // Change constant term 'a' to get an approximation to order 10^{-a-1}
//     double x2 = x * x;
//     double xalpha = exp(alpha * log(x));
//
//     arma::vec evens = 2.0 * arma::regspace(1, n);
//     arma::vec odds = 2.0 * arma::regspace(1, n) - 1.0;
//     arma::vec summands = arma::cumprod(-x2 / (evens % odds)) / (evens + alpha);
//
//     return xalpha * (1.0 / alpha + arma::accu(summands));
// }
//
// double Si( double x, double alpha ) {
//
//     if (x < 0.0)
//         Rcpp::stop("ERROR in Si: 'x' cannot be negative.");
//
//     if (x == 0.0)
//         return 0.0;
//
//     if (alpha <= 0.0 || alpha >= 1.0) {
//         Rcpp::stop("ERROR in Si: 'alpha' must be between 0 and 1 strictly.");
//     }
//
//     int n = std::ceil(7.0 + 1.36 * x); // Change constant term 'a' to get an approximation to order 10^{-a-1}
//     double x2 = x * x;
//     double xalpha = exp(alpha * log(x));
//
//     arma::vec evens = 2.0 * arma::regspace(1, n);
//     arma::vec odds = 2.0 * arma::regspace(1, n) + 1.0;
//     arma::vec summands = arma::cumprod(-x2 / (evens % odds)) / (odds + alpha);
//
//     return x * xalpha * (1.0 / (1.0 + alpha) + arma::accu(summands));
// }
//
// double taylorf( double x, double alpha ) {
//
//     if (x < 20.0)
//         Rcpp::stop("ERROR in taylorf: 'x' must be above 20 for correct approximation.");
//
//     if (alpha <= 0.0 || alpha >= 1.0) {
//         Rcpp::stop("ERROR in taylorf: 'alpha' must be between 0 and 1 strictly.");
//     }
//
//     int n = 10;
//     double inv_x = 1.0 / x;
//     double inv_x2 = inv_x * inv_x;
//     double xalpha = exp(alpha * log(x));
//
//     arma::vec evens = alpha - 2.0 * arma::regspace(1, n);
//     arma::vec odds = alpha - 2.0 * arma::regspace(1, n) + 1.0;
//     arma::vec summands = arma::cumprod(- inv_x2 * evens % odds);
//
//     return (xalpha * inv_x) * (1.0 + arma::accu(summands));
//
// }
//
// double taylorg( double x, double alpha ) {
//
//     if (x < 20.0)
//         Rcpp::stop("ERROR in taylorg: 'x' must be above 20 for correct approximation.");
//
//     if (alpha <= 0.0 || alpha >= 1.0) {
//         Rcpp::stop("ERROR in taylorg: 'alpha' must be between 0 and 1 strictly.");
//     }
//
//     int n = 10;
//     double inv_x = 1.0 / x;
//     double inv_x2 = inv_x * inv_x;
//     double xalpha = exp(alpha * log(x));
//
//     arma::vec evens = alpha - 2.0 * arma::regspace(1, n-1);
//     arma::vec odds = alpha - 2.0 * arma::regspace(1, n-1) - 1.0;
//     arma::vec summands = arma::cumprod(- inv_x2 * evens % odds);
//
//     return (1.0 - alpha) * (xalpha * inv_x2) * (1.0 + arma::accu(summands));
//
// }

arma::cx_double inc_gamma_imag( double x, double alpha ) {

    if (x < 0.0)
        Rcpp::stop("ERROR in inc_gamma_imag: 'x' cannot be negative.");

    if (x == 0.0)
        return exp(-.5*i*arma::datum::pi*alpha) * boost::math::tgamma(alpha);

    if (alpha <= 0.0 || alpha >= 1.0) {
        Rcpp::stop("ERROR in inc_gamma_imag: 'alpha' must be between 0 and 1 strictly.");
    }

    // Taylor expansion of Ci and Si around x = 0, up to 34 terms calculated for x = 20
    // This approximation is of order 10^{-8}
    if (x < 20.0) {
        int n = std::ceil(7.0 + 1.36 * x); // Change constant term 'a' to get an approximation to order 10^{-a-1}
        double x2 = x * x;
        double xalpha = exp(alpha * log(x));

        arma::vec evens = 2 * arma::regspace(1, n);
        arma::vec oddsm1 = 2 * arma::regspace(1, n) - 1;
        arma::vec oddsp1 = 2 * arma::regspace(1, n) + 1;
        arma::vec Ci_summands = arma::cumprod(-x2 / (evens % oddsm1)) / (evens + alpha);
        arma::vec Si_summands = arma::cumprod(-x2 / (evens % oddsp1)) / (oddsp1 + alpha);

        double Ci = xalpha * (1.0 / alpha + arma::accu(Ci_summands));
        double Si = x * xalpha * (1.0 / (1.0 + alpha) + arma::accu(Si_summands));

        return exp(-.5*i*arma::datum::pi*alpha) * boost::math::tgamma(alpha) - Ci + i*Si;
    }

    // Taylor expansion of f and g around x -> infty, only 10 terms calculated
    // This approximation is at most of order 10^{-8}
    // Can't be used for x < 20 because then it is of lower order
    int n = 10;
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    double xalpha = exp(alpha * log(x));

    arma::vec evens = alpha - 2.0 * arma::regspace(1, n);
    arma::vec oddsf = alpha - 2.0 * arma::regspace(1, n) + 1.0;
    arma::vec oddsg = alpha - 2.0 * arma::regspace(1, n-1) - 1.0;
    arma::vec summands_f = arma::cumprod(- inv_x2 * evens % oddsf);
    arma::vec summands_g = arma::cumprod(- inv_x2 * evens.subvec(0, n-2) % oddsg);

    double taylorf = (xalpha * inv_x) * (1.0 + arma::accu(summands_f));
    double taylorg = (1.0 - alpha) * (xalpha * inv_x2) * (1.0 + arma::accu(summands_g));

    double Ci = taylorg * cos(x) - taylorf * sin(x);
    double Si = taylorg * sin(x) + taylorf * cos(x);

    return Ci - i*Si;

}

// Powers of 10
double quick_pow10(int n)
{
    static double pow10[14] = {
        1.0, 10.0, 100.0, 1000.0, 10000.0,
        100000.0, 1000000.0, 10000000.0, 100000000.0, 1000000000.0,
        10000000000.0, 100000000000.0, 1000000000000.0, 10000000000000.0
    };

    return pow10[n];
}

double quick_negpow10(int n)
{
    static double pow10[19] = {
        1.0, 0.1, 0.01, 0.001, 0.0001,
        0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001,
        0.0000000001, 0.00000000001, 0.000000000001, 0.0000000000001, 0.00000000000001,
        0.000000000000001, 0.0000000000000001, 0.00000000000000001, 0.000000000000000001
    };

    return pow10[n];
}

// Padé approximants
// cf. https://en.wikipedia.org/wiki/Trigonometric_integral
// Rowe, B. T. P., Jarvis, M., Mandelbaum, R., Bernstein, G. M., Bosch, J., Simet, M., … Gill, M. S. S. (2015). GalSim: The modular galaxy image simulation toolkit. Astronomy and Computing, 10, 121–150. https://doi.org/10.1016/j.ascom.2015.02.002
double padef( double x ) {
    double inv_x = 1.0 / x,
           inv_x2 = inv_x * inv_x;

    double num = 1 +
        inv_x2 * (7.44437068161936700618 * quick_pow10(2) +
        inv_x2 * (1.96396372895146869801 * quick_pow10(5) +
        inv_x2 * (2.37750310125431834034 * quick_pow10(7) +
        inv_x2 * (1.43073403821274636888 * quick_pow10(9) +
        inv_x2 * (4.33736238870432522765 * quick_pow10(10) +
        inv_x2 * (6.40533830574022022911 * quick_pow10(11) +
        inv_x2 * (4.20968180571076940208 * quick_pow10(12) +
        inv_x2 * (1.00795182980368574617 * quick_pow10(13) +
        inv_x2 * (4.94816688199951963482 * quick_pow10(12) +
        inv_x2 * (- 4.94701168645415959931 * quick_pow10(11)
    ))))))))));

    double denom = 1 +
        inv_x2 * (7.46437068161927678031 * quick_pow10(2) +
        inv_x2 * (1.97865247031583951450 * quick_pow10(5) +
        inv_x2 * (2.41535670165126845144 * quick_pow10(7) +
        inv_x2 * (1.47478952192985464958 * quick_pow10(9) +
        inv_x2 * (4.58595115847765779830 * quick_pow10(10) +
        inv_x2 * (7.08501308149515401563 * quick_pow10(11) +
        inv_x2 * (5.06084464593475076774 * quick_pow10(12) +
        inv_x2 * (1.43468549171581016479 * quick_pow10(13) +
        inv_x2 * (1.11535493509914254097 * quick_pow10(13)
    )))))))));

    return inv_x * num / denom;
}

double padeg( double x ) {
    double inv_x = 1.0 / x,
           inv_x2 = inv_x * inv_x;

    double num = 1 +
        inv_x2 * (8.1359520115168615 * quick_pow10(2) +
        inv_x2 * (2.35239181626478200 * quick_pow10(5) +
        inv_x2 * (3.12557570795778731 * quick_pow10(7) +
        inv_x2 * (2.06297595146763354 * quick_pow10(9) +
        inv_x2 * (6.83052205423625007 * quick_pow10(10) +
        inv_x2 * (1.09049528450362786 * quick_pow10(12) +
        inv_x2 * (7.57664583257834349 * quick_pow10(12) +
        inv_x2 * (1.81004487464664575 * quick_pow10(13) +
        inv_x2 * (6.43291613143049485 * quick_pow10(12) +
        inv_x2 * (- 1.36517137670871689 * quick_pow10(12)
    ))))))))));


    double denom = 1 +
        inv_x2 * (8.19595201151451564 * quick_pow10(2) +
        inv_x2 * (2.40036752835578777 * quick_pow10(5) +
        inv_x2 * (3.26026661647090822 * quick_pow10(7) +
        inv_x2 * (2.23355543278099360 * quick_pow10(9) +
        inv_x2 * (7.87465017341829930 * quick_pow10(10) +
        inv_x2 * (1.39866710696414565 * quick_pow10(12) +
        inv_x2 * (1.17164723371736605 * quick_pow10(13) +
        inv_x2 * (4.01839087307656620 * quick_pow10(13) +
        inv_x2 * (3.99653257887490811 * quick_pow10(13)
    )))))))));

    return inv_x2 * num / denom;
}

double Ci( double x ) {

    if (x < 0)
        Rcpp::stop("ERROR in Ci: 'x' cannot be negative.");

    if (x <= 4) {
        double x2 = x * x;

        double num = -0.25 +
            x2 * (7.51851524438898291 * quick_negpow10(3) +
            x2 * (- 1.27528342240267686 * quick_negpow10(4) +
            x2 * (1.05297363846239184 * quick_negpow10(6) +
            x2 * (- 4.68889508144848019 * quick_negpow10(9) +
            x2 * (1.06480802891189243 * quick_negpow10(11) +
            x2 * (- 9.93728488857585407 * quick_negpow10(15)
        ))))));

        double denom = 1 +
            x2 * (1.1592605689110735 * quick_negpow10(2) +
            x2 * (6.72126800814254432 * quick_negpow10(5) +
            x2 * (2.55533277086129636 * quick_negpow10(7) +
            x2 * (6.97071295760958946 * quick_negpow10(10) +
            x2 * (1.38536352772778619 * quick_negpow10(12) +
            x2 * (1.89106054713059759 * quick_negpow10(15) +
            x2 * (1.39759616731376855 * quick_negpow10(18)
        )))))));

        return arma::datum::euler + log(x) + x2 * num / denom;
    }

    return padef(x) * sin(x) - padeg(x) * cos(x);

}

double Si( double x ) {

    if (x < 0)
        Rcpp::stop("ERROR in Si: 'x' cannot be negative.");

    if (x <= 4) {
        double x2 = x * x;

        double num = 1 +
            x2 * (- 4.54393409816329991 * quick_negpow10(2) +
            x2 * (1.15457225751016682 * quick_negpow10(3) +
            x2 * (- 1.41018536821330254 * quick_negpow10(5) +
            x2 * (9.43280809438713025 * quick_negpow10(8) +
            x2 * (- 3.53201978997168357 * quick_negpow10(10) +
            x2 * (7.08240282274875911 * quick_negpow10(13) +
            x2 * (- 6.05338212010422477 * quick_negpow10(16)
        )))))));

        double denom = 1 +
            x2 * (1.01162145739225565 * quick_negpow10(2) +
            x2 * (4.99175116169755106 * quick_negpow10(5) +
            x2 * (1.55654986308745614 * quick_negpow10(7) +
            x2 * (3.28067571055789734 * quick_negpow10(10) +
            x2 * (4.5049097575386581 * quick_negpow10(13) +
            x2 * (3.21107051193712168 * quick_negpow10(16)
        ))))));

        return x * num / denom;
    }

    return 0.5 * arma::datum::pi - padef(x) * cos(x) - padeg(x) * sin(x);

}

arma::cx_double E1_imaginary( double x ) {

    if (x < 0)
        Rcpp::stop("ERROR in E1_imaginary: 'x' cannot be negative.");

    return i * (- 0.5 * arma::datum::pi + Si(x)) - Ci(x);

}

// // Contour integration for the incomplete gamma function with imaginary argument
// // Integrate a continuous function on interval [a, b]
// double integral_midpoint(double(*f)(double x), double a, double b, int n) {
//     double step = (b - a) / n;  // width of each small rectangle
//     double area = 0.0;  // signed area
//     for (int i = 0; i < n; i ++) {
//         area += f(a + (i + 0.5) * step); // sum up each small rectangle
//     }
//     return (area * step);
// }
//
// double integral_midpoint(double(*f)(double x, double nu), double a, double b, int n, double nu) {
//     double step = (b - a) / n;  // width of each small rectangle
//     double area = 0.0;  // signed area
//     for (int i = 0; i < n; i ++) {
//         area += f(a + (i + 0.5) * step, nu); // sum up each small rectangle
//     }
//     return (area * step);
// }
//
// double integral_simpson(double(*f)(double x), double a, double b, int n) {
//     if (n % 2 == 1)
//         n++;
//     double step = (b - a) / n;  // width of each small rectangle
//     double area = f(a) + f(b);  // first and last indices
//     for (int i = 1; i < n; i++) {
//         if (i % 2 == 1)
//             area += 4 * f(a + i * step);
//         else
//             area += 2 * f(a + i * step);
//     }
//     return (area * step / 3.0);
// }
//
// double integral_simpson(double(*f)(double x, double nu), double a, double b, int n, double nu) {
//     if (n % 2 == 1)
//         n++;
//     double step = (b - a) / n;  // width of each small rectangle
//     double area = f(a, nu) + f(b, nu);  // first and last indices
//     for (int i = 1; i < n; i++) {
//         if (i % 2 == 1)
//             area += 4 * f(a + i * step, nu);
//         else
//             area += 2 * f(a + i * step, nu);
//     }
//     return (area * step / 3.0);
// }
//
// double integral_simpson(double(*f)(double x, double nu, double r), double a, double b, int n, double nu, double r) {
//     if (n % 2 == 1)
//         n++;
//     double step = (b - a) / n;  // width of each small rectangle
//     double area = f(a, nu, r) + f(b, nu, r);  // first and last indices
//     for (int i = 1; i < n; i++) {
//         if (i % 2 == 1)
//             area += 4 * f(a + i * step, nu, r);
//         else
//             area += 2 * f(a + i * step, nu, r);
//     }
//     return (area * step / 3.0);
// }
//
// // Real and imaginary part for the contour integral on the quadrant with radius r > 0
// double quadrant_real(double x, double nu, double r) {
//     return exp(-r*cos(x)) * cos(x*nu - r*sin(x));
// }
//
// double quadrant_imag(double x, double nu, double r) {
//     return exp(-r*cos(x)) * sin(x*nu - r*sin(x));
// }
//
// arma::cx_double contour_quadrant(double nu, double r) {
//     double real_part = integral_simpson(quadrant_real, 0.0, .5 * arma::datum::pi, 100, nu, r);
//     double imag_part = integral_simpson(quadrant_imag, 0.0, .5 * arma::datum::pi, 100, nu, r);
//     return exp(nu*log(r)) * arma::cx_double(real_part, imag_part);
// }
//
// arma::cx_double inc_gamma_imag(double nu, double r) {
//
//     if (r < 0)
//         Rcpp::stop("ERROR in inc_gamma_imag: 'r' cannot be negative.");
//     if (nu <= 0)
//         Rcpp::stop("ERROR in inc_gamma_imag: 'nu' must be between 0 and 1.");
//     if (nu >= 1)
//         Rcpp::stop("ERROR in inc_gamma_imag: 'nu' must be between 0 and 1.");
//
//     arma::cx_double term1 = exp(-.5*i*arma::datum::pi*nu) * boost::math::tgamma(nu, r);
//     arma::cx_double term2 = exp(-.5*i*arma::datum::pi*(nu-1)) * contour_quadrant(nu, r);
//     return  term1 - term2;
//
// }
