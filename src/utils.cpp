#include "utils.hpp"
#include <string>
#include <iostream>

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
};

// Powers of 10
double quick_pow10(int n)
{
    static double pow10[14] = {
        1.0, 10.0, 100.0, 1000.0, 10000.0,
        100000.0, 1000000.0, 10000000.0, 100000000.0, 1000000000.0,
        10000000000.0, 100000000000.0, 1000000000000.0, 10000000000000.0
    };

    return pow10[n];
};

double quick_negpow10(int n)
{
    static double pow10[19] = {
        1.0, 0.1, 0.01, 0.001, 0.0001,
        0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001,
        0.0000000001, 0.00000000001, 0.000000000001, 0.0000000000001, 0.00000000000001,
        0.000000000000001, 0.0000000000000001, 0.00000000000000001, 0.000000000000000001
    };

    return pow10[n];
};

// Padé approximants
// cf. https://en.wikipedia.org/wiki/Trigonometric_integral
// Rowe, B. T. P., Jarvis, M., Mandelbaum, R., Bernstein, G. M., Bosch, J., Simet, M., … Gill, M. S. S. (2015). GalSim: The modular galaxy image simulation toolkit. Astronomy and Computing, 10, 121–150. https://doi.org/10.1016/j.ascom.2015.02.002
double padef( double x ) {
    double inv_x = 1.0 / x,
           inv_x2 = inv_x * inv_x,
           inv_x4 = inv_x2 * inv_x2,
           inv_x6 = inv_x4 * inv_x2,
           inv_x8 = inv_x6 * inv_x2,
           inv_x10 = inv_x8 * inv_x2,
           inv_x12 = inv_x10 * inv_x2,
           inv_x14 = inv_x12 * inv_x2,
           inv_x16 = inv_x14 * inv_x2,
           inv_x18 = inv_x16 * inv_x2,
           inv_x20 = inv_x18 * inv_x2;

    double num = 1 + 7.44437068161936700618 * quick_pow10(2) * inv_x2
                   + 1.96396372895146869801 * quick_pow10(5) * inv_x4
                   + 2.37750310125431834034 * quick_pow10(7) * inv_x6
                   + 1.43073403821274636888 * quick_pow10(9) * inv_x8
                   + 4.33736238870432522765 * quick_pow10(10) * inv_x10
                   + 6.40533830574022022911 * quick_pow10(11) * inv_x12
                   + 4.20968180571076940208 * quick_pow10(12) * inv_x14
                   + 1.00795182980368574617 * quick_pow10(13) * inv_x16
                   + 4.94816688199951963482 * quick_pow10(12) * inv_x18
                   - 4.94701168645415959931 * quick_pow10(11) * inv_x20;

    double denom = 1 + 7.46437068161927678031 * quick_pow10(2) * inv_x2
                     + 1.97865247031583951450 * quick_pow10(5) * inv_x4
                     + 2.41535670165126845144 * quick_pow10(7) * inv_x6
                     + 1.47478952192985464958 * quick_pow10(9) * inv_x8
                     + 4.58595115847765779830 * quick_pow10(10) * inv_x10
                     + 7.08501308149515401563 * quick_pow10(11) * inv_x12
                     + 5.06084464593475076774 * quick_pow10(12) * inv_x14
                     + 1.43468549171581016479 * quick_pow10(13) * inv_x16
                     + 1.11535493509914254097 * quick_pow10(13) * inv_x18;

    return inv_x * num / denom;
};

double padeg( double x ) {
    double inv_x = 1.0 / x,
           inv_x2 = inv_x * inv_x,
           inv_x4 = inv_x2 * inv_x2,
           inv_x6 = inv_x4 * inv_x2,
           inv_x8 = inv_x6 * inv_x2,
           inv_x10 = inv_x8 * inv_x2,
           inv_x12 = inv_x10 * inv_x2,
           inv_x14 = inv_x12 * inv_x2,
           inv_x16 = inv_x14 * inv_x2,
           inv_x18 = inv_x16 * inv_x2,
           inv_x20 = inv_x18 * inv_x2;

    double num = 1 + 8.1359520115168615 * quick_pow10(2) * inv_x2
                   + 2.35239181626478200 * quick_pow10(5) * inv_x4
                   + 3.12557570795778731 * quick_pow10(7) * inv_x6
                   + 2.06297595146763354 * quick_pow10(9) * inv_x8
                   + 6.83052205423625007 * quick_pow10(10) * inv_x10
                   + 1.09049528450362786 * quick_pow10(12) * inv_x12
                   + 7.57664583257834349 * quick_pow10(12) * inv_x14
                   + 1.81004487464664575 * quick_pow10(13) * inv_x16
                   + 6.43291613143049485 * quick_pow10(12) * inv_x18
                   - 1.36517137670871689 * quick_pow10(12) * inv_x20;


    double denom = 1 + 8.19595201151451564 * quick_pow10(2) * inv_x2
                     + 2.40036752835578777 * quick_pow10(5) * inv_x4
                     + 3.26026661647090822 * quick_pow10(7) * inv_x6
                     + 2.23355543278099360 * quick_pow10(9) * inv_x8
                     + 7.87465017341829930 * quick_pow10(10) * inv_x10
                     + 1.39866710696414565 * quick_pow10(12) * inv_x12
                     + 1.17164723371736605 * quick_pow10(13) * inv_x14
                     + 4.01839087307656620 * quick_pow10(13) * inv_x16
                     + 3.99653257887490811 * quick_pow10(13) * inv_x18;

    return inv_x2 * num / denom;
};

double Ci( double x ) {
    if (x < 0)
        throw "ERROR in Ci: 'x' cannot be negative.";
    if (x <= 4) {
        double x2 = x * x,
               x4 = x2 * x2,
               x6 = x4 * x2,
               x8 = x6 * x2,
               x10 = x8 * x2,
               x12 = x10 * x2,
               x14 = x12 * x2;

        double num = -0.25 + 7.51851524438898291 * quick_negpow10(3) * x2
                           - 1.27528342240267686 * quick_negpow10(4) * x4
                           + 1.05297363846239184 * quick_negpow10(6) * x6
                           - 4.68889508144848019 * quick_negpow10(9) * x8
                           + 1.06480802891189243 * quick_negpow10(11) * x10
                           - 9.93728488857585407 * quick_negpow10(15) * x12;

        double denom = 1 + 1.1592605689110735 * quick_negpow10(2) * x2
                         + 6.72126800814254432 * quick_negpow10(5) * x4
                         + 2.55533277086129636 * quick_negpow10(7) * x6
                         + 6.97071295760958946 * quick_negpow10(10) * x8
                         + 1.38536352772778619 * quick_negpow10(12) * x10
                         + 1.89106054713059759 * quick_negpow10(15) * x12
                         + 1.39759616731376855 * quick_negpow10(18) * x14;

        return arma::datum::euler + log(x) + x2 * num / denom;
    }

    return padef(x) * sin(x) - padeg(x) * cos(x);
};

double Si( double x ) {
    if (x < 0)
        throw "ERROR in Si: 'x' cannot be negative.";
    if (x <= 4) {
        double x2 = x * x,
            x4 = x2 * x2,
            x6 = x4 * x2,
            x8 = x6 * x2,
            x10 = x8 * x2,
            x12 = x10 * x2,
            x14 = x12 * x2;

        double num = 1 - 4.54393409816329991 * quick_negpow10(2) * x2
            + 1.15457225751016682 * quick_negpow10(3) * x4
            - 1.41018536821330254 * quick_negpow10(5) * x6
            + 9.43280809438713025 * quick_negpow10(8) * x8
            - 3.53201978997168357 * quick_negpow10(10) * x10
            + 7.08240282274875911 * quick_negpow10(13) * x12
            - 6.05338212010422477 * quick_negpow10(16) * x14;

        double denom = 1 + 1.01162145739225565 * quick_negpow10(2) * x2
            + 4.99175116169755106 * quick_negpow10(5) * x4
            + 1.55654986308745614 * quick_negpow10(7) * x6
            + 3.28067571055789734 * quick_negpow10(10) * x8
            + 4.5049097575386581 * quick_negpow10(13) * x10
            + 3.21107051193712168 * quick_negpow10(16) * x12;

        return x * num / denom;
    }

    return 0.5 * arma::datum::pi - padef(x) * cos(x) - padeg(x) * sin(x);
};

arma::cx_vec E1_imaginary( arma::vec x ) {
    arma::cx_vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::cx_vec::iterator it_y = y.begin();

    try {
        // Loop on x
        for (; it_x != it_x_end; ++it_x, ++it_y) {
            *it_y = i * (- 0.5 * arma::datum::pi + Si(*it_x)) - Ci(*it_x);
        }
    } catch (const char* msg) {
        std::cerr << msg << std::endl;
    }

    return y;
};
