#include "model.hpp"
#include "utils.hpp"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) * theta * param(2)^theta * t^(-theta-1)
// and theta = 3
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate /in [0,1]
// param(2): scale parameter of the Pareto distribution function

////////////////////////////////////////////////
// Methods for long term mean and its derivative

double Pareto3::mean() {
    return param(0) / ( 1.0 - param(1) );
}

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec Pareto3::h( arma::vec x ) {
    arma::vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::vec::iterator it_y = y.begin();

    // Temp variables
    double term = 3 * param(1) * param(2) * param(2) * param(2);
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
    arma::vec xia = xi * param(2),
              xia2 = xia % xia,
              xia3 = xia2 % xia;
    arma::cx_vec term1 = param(1) * exp(- i * xia) % (1.0 - 0.5 * i * xia - 0.5 * xia2);
    arma::cx_vec term2 = 0.5 * param(1) * i * xia3 % E1_imaginary(xia);
    return term1 + term2;
}
