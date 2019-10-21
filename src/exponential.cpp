#include "model.hpp"
#include "utils.hpp"

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec Exponential::h( arma::vec x ) {
    return param(1) * param(2) * arma::exp( - param(2) * x );
};

arma::cx_vec Exponential::H( arma::vec xi ) {
    arma::vec factor = param(1) * param(2) / ( param(2)*param(2) + xi%xi );
    arma::cx_vec zeta = arma::cx_vec( factor * param(2), - factor % xi );
    return zeta;
};

arma::cx_mat Exponential::dH( arma::vec xi ) {
    arma::cx_mat grad = arma::zeros<arma::cx_mat>( xi.n_elem, param.n_elem );
    arma::cx_vec denom = 1.0 / (param(2) + i * xi);
    grad.col(1) = param(2) * denom;
    grad.col(2) = i * param(1) * (xi % denom % denom);
    return grad;
};

arma::cx_cube Exponential::ddH( arma::vec xi ) {
    arma::cx_cube hess = arma::zeros<arma::cx_cube>( param.n_elem, param.n_elem, xi.n_elem );
    arma::cx_vec denom = 1.0 / (param(2) + i * xi);
    arma::cx_vec grad12 = i * (xi % denom % denom);
    hess.tube(1,2) = grad12;
    hess.tube(2,1) = grad12;
    hess.tube(2,2) = -2.0 * param(1) * grad12 % denom;
    return hess;
};
