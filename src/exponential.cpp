#include "model.hpp"
#include "utils.hpp"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) * param(2) * exp( -param(2) * t )
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate /in [0,1]
// param(2): time scale parameter of the exponential fertility distribution

////////////////////////////////////////////////
// Methods for long term mean and its derivative

double Exponential::mean() {
    return param(0) / ( 1.0 - param(1) );
}

arma::vec Exponential::dmean() {
    double denom = 1.0 / ( 1.0 - param(1) );
    arma::vec grad = { denom, param(0) * denom * denom, 0 };
    return grad;
}

arma::mat Exponential::ddmean() {
    double denom = 1.0 / ( 1.0 - param(1) );
    double denom2 = denom * denom;
    arma::mat hess = {  {   0.0,                  denom2, 0.0},
                        {denom2, 2*param(0)*denom2*denom, 0.0},
                        {   0.0,                     0.0, 0.0} };
    return hess;
}

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec Exponential::h( arma::vec x ) {
    return param(1) * param(2) * arma::exp( - param(2) * x );
}

arma::cx_vec Exponential::H( arma::vec xi ) {
    arma::vec factor = param(1) * param(2) / ( param(2)*param(2) + xi%xi );
    arma::cx_vec zeta = arma::cx_vec( factor * param(2), - factor % xi );
    return zeta;
}

arma::cx_mat Exponential::dH( arma::vec xi ) {
    arma::cx_mat grad = arma::zeros<arma::cx_mat>( xi.n_elem, param.n_elem );
    arma::cx_vec denom = 1.0 / (param(2) + i * xi);
    grad.col(1) = param(2) * denom;
    grad.col(2) = i * param(1) * (xi % denom % denom);
    return grad;
}

arma::cx_cube Exponential::ddH( arma::vec xi ) {
    arma::cx_cube hess = arma::zeros<arma::cx_cube>( param.n_elem, param.n_elem, xi.n_elem );
    arma::cx_vec denom = 1.0 / (param(2) + i * xi);
    arma::cx_vec grad12 = i * (xi % denom % denom);
    hess.tube(1,2) = grad12;
    hess.tube(2,1) = grad12;
    hess.tube(2,2) = -2.0 * param(1) * grad12 % denom;
    return hess;
}
