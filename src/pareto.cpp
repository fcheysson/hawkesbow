#include "model.h"
#include "utils.h"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) * theta * param(2)^theta * t^(-theta-1)
// and theta = 3
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate /in [0,1]
// param(2): scale parameter of the Pareto distribution function

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec Pareto3::h( arma::vec x ) {
    arma::vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::vec::iterator it_y = y.begin();

    // Temp variables
    double term = 3.0 * param(1) * param(2) * param(2) * param(2);
    double x2;

    // Loop on x
    for (; it_x != it_x_end; ++it_x, ++it_y) {
        if (*it_x < param(2))
            *it_y = 0.0;
        else {
            x2 = *it_x * *it_x;
            *it_y = term / (x2 * x2);
        }
    }
    return y;
}

arma::cx_vec Pareto3::H( arma::vec xi ) {
    arma::cx_vec y(xi.n_elem);

    // Iterators
    arma::vec::iterator it_xi = xi.begin();
    arma::vec::iterator it_xi_end = xi.end();
    arma::cx_vec::iterator it_y = y.begin();

    // Temp variables
    double xia, xia2, xia3;
    arma::cx_double term1, term2;

    // Loop on xi
    for (; it_xi != it_xi_end; ++it_xi, ++it_y) {
        if (*it_xi >= 0) {
            xia = *it_xi * param(2);
            xia2 = xia * xia;
            xia3 = xia2 * xia;
            term1 = param(1) * exp(- i * xia) * (1.0 - 0.5 * i * xia - 0.5 * xia2);
            term2 = 0.5 * param(1) * i * xia3 * E1_imaginary(xia);
            *it_y = term1 + term2;
        }
        if (*it_xi < 0) {   // take conjugate
            xia = - *it_xi * param(2);
            xia2 = xia * xia;
            xia3 = xia2 * xia;
            term1 = param(1) * exp(- i * xia) * (1.0 - 0.5 * i * xia - 0.5 * xia2);
            term2 = 0.5 * param(1) * i * xia3 * E1_imaginary(xia);
            *it_y = std::conj(term1 + term2);
        }
    }

    return y;
}


////////////////////////////////////////////////
// Pareto2

arma::vec Pareto2::h( arma::vec x ) {
    arma::vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::vec::iterator it_y = y.begin();

    // Temp variables
    double term = 2.0 * param(1) * param(2) * param(2);

    // Loop on x
    for (; it_x != it_x_end; ++it_x, ++it_y) {
        if (*it_x < param(2))
            *it_y = 0.0;
        else {
            *it_y = term / (*it_x * *it_x * *it_x);
        }
    }
    return y;
}

arma::cx_vec Pareto2::H( arma::vec xi ) {
    arma::cx_vec y(xi.n_elem);

    // Iterators
    arma::vec::iterator it_xi = xi.begin();
    arma::vec::iterator it_xi_end = xi.end();
    arma::cx_vec::iterator it_y = y.begin();

    // Temp variables
    double xia, xia2;
    arma::cx_double term1, term2;

    // Loop on xi
    for (; it_xi != it_xi_end; ++it_xi, ++it_y) {
        if (*it_xi >= 0) {
            xia = *it_xi * param(2);
            xia2 = xia * xia;
            term1 = param(1) * exp(- i * xia) * (1.0 - i * xia);
            term2 = param(1) * xia2 * E1_imaginary(xia);
            *it_y = term1 - term2;
        }
        if (*it_xi < 0) {   // take conjugate
            xia = - *it_xi * param(2);
            xia2 = xia * xia;
            term1 = param(1) * exp(- i * xia) * (1.0 - i * xia);
            term2 = param(1) * xia2 * E1_imaginary(xia);
            *it_y = std::conj(term1 - term2);
        }
    }

    return y;
}

////////////////////////////////////////////////
// Pareto1

arma::vec Pareto1::h( arma::vec x ) {
    arma::vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::vec::iterator it_y = y.begin();

    // Temp variables
    double term = param(1) * param(2);

    // Loop on x
    for (; it_x != it_x_end; ++it_x, ++it_y) {
        if (*it_x < param(2))
            *it_y = 0.0;
        else {
            *it_y = term / (*it_x * *it_x);
        }
    }
    return y;
}

arma::cx_vec Pareto1::H( arma::vec xi ) {
    arma::cx_vec y(xi.n_elem);

    // Iterators
    arma::vec::iterator it_xi = xi.begin();
    arma::vec::iterator it_xi_end = xi.end();
    arma::cx_vec::iterator it_y = y.begin();

    // Temp variables
    double xia;
    arma::cx_double term1, term2;

    // Loop on xi
    for (; it_xi != it_xi_end; ++it_xi, ++it_y) {
        if (*it_xi >= 0) {
            xia = *it_xi * param(2);
            term1 = param(1) * exp(- i * xia);
            term2 = param(1) * i * xia * E1_imaginary(xia);
            *it_y = term1 - term2;
        }
        if (*it_xi < 0) {   // take conjugate
            xia = - *it_xi * param(2);
            term1 = param(1) * exp(- i * xia);
            term2 = param(1) * i * xia * E1_imaginary(xia);
            *it_y = std::conj(term1 - term2);
        }
    }

    return y;
}
