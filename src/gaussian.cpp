#include "model.h"
#include "utils.h"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) / sqrt(2 * pi * param(3)) * exp(-(t-param(2))^2 / (2 * param(3)))
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate /in [0,1]
// param(2): mean parameter of the Gaussian distribution function
// param(3): variance parameter of the Gaussian distribution function

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec Gaussian::h( arma::vec x ) {
    return param(1) * arma::exp(- .5 * (x - param(2)) % (x - param(2)) / param(3)) / sqrt(2 * arma::datum::pi * param(3));
}

arma::cx_vec Gaussian::H( arma::vec xi ) {
    return param(1) * arma::exp(- i * param(2) * xi) % arma::exp(- .5 * param(3) * xi % xi);
}
