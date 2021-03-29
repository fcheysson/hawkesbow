#include "model.h"
#include "utils.h"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) * param(2) * exp( -param(2) * t )
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate /in [0,1]
// param(2): time scale parameter of the exponential fertility distribution

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

////////////////////////////////////////////
// Methods for maximum likelihood estimation

double Exponential::loglik( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double beta = param(2);

    // Sum log \lambda*
    double A = 0.0;
    double part1 = log(eta);
    for (arma::uword i = 1; i < n; i++) {
        A = exp(- beta * (events(i) - events(i-1))) * (1.0 + A);
        part1 += log(eta + mu * beta * A);
    }

    // Int \lambda* dt
    double part2 = eta * end + mu * ((double)n - exp(-beta * (end - events(n-1))) * (1.0 + A));

    return part1 - part2;
}

arma::vec Exponential::dloglik( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double beta = param(2);

    // Fill for i = 0
    arma::vec grad = { 1.0/eta, 0.0, 0.0 };

    // Iterate on arrival times
    double A = 0.0;
    double C = 0.0;
    double B;
    double denom, expint;
    for (arma::uword i = 1; i < n; i++) {
        expint = exp(- beta * (events(i) - events(i-1)));
        A = expint * (1.0 + A);
        C = expint * (events(i-1) + C);
        B = events(i) * A - C;
        denom = 1.0 / (eta + mu * beta * A);
        grad(0) += 1.0 * denom;
        grad(1) += beta * A * denom;
        grad(2) += (mu * A - mu * beta * B) * denom;
    }

    expint = exp(- beta * (end - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = end * A - C;
    grad(0) -= end;
    grad(1) -= ( (double)n - A );
    grad(2) -= mu * B;

    return grad;
}

arma::mat Exponential::ddloglik( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double mu2 = mu * mu;
    const double beta = param(2);
    const double beta2 = beta * beta;

    // Fill for i = 0
    arma::mat hess = { {-1.0/eta, 0.0, 0.0},
    {     0.0, 0.0, 0.0},
    {     0.0, 0.0, 0.0} };

    // Iterate on arrival times
    double A = 0.0;
    double C = 0.0;
    double B;
    double E = 0.0;
    double D;
    double denom, denom2, expint;
    for (arma::uword i = 1; i < n; i++) {
        expint = exp(- beta * (events(i) - events(i-1)));
        A = expint * (1.0 + A);
        C = expint * (events(i-1) + C);
        B = events(i) * A - C;
        E = expint * (events(i-1) * events(i-1) + E);
        D = events(i) * events(i) * A + E - 2.0 * events(i) * C;
        denom = 1.0 / (eta + mu * beta * A);
        denom2 = denom * denom;
        hess(0,0) -= 1.0 * denom2;
        hess(0,1) -= beta * A * denom2;
        hess(0,2) -= mu * (A - beta * B) * denom2;
        hess(1,1) -= beta2 * A * A * denom2;
        hess(1,2) += eta * (A - beta * B) * denom2;
        hess(2,2) += (eta * mu * (beta * D - 2 * B) + mu2 * (beta2 * (A * E - C * C) - A * A)) * denom2;
    }

    expint = exp(- beta * (end - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = end * A - C;
    E = expint * (events(n-1) * events(n-1) + E);
    D = end * end * A + E - 2.0 * end * C;
    hess(1,2) -= B;
    hess(2,2) += mu * D;

    // Symm
    hess(1,0) = hess(0,1);
    hess(2,0) = hess(0,2);
    hess(2,1) = hess(1,2);

    return hess;
}

Rcpp::List Exponential::loglikngrad( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double beta = param(2);

    // Fill for i = 0
    double lik = log(eta);
    arma::vec grad = { 1.0/eta, 0.0, 0.0 };

    // Iterate on arrival times
    double A = 0.0;
    double C = 0.0;
    double B;
    double denom, expint;
    for (arma::uword i = 1; i < n; i++) {
        expint = exp(- beta * (events(i) - events(i-1)));
        A = expint * (1.0 + A);
        C = expint * (events(i-1) + C);
        B = events(i) * A - C;
        denom = 1.0 / (eta + mu * beta * A);
        lik += log(eta + mu * beta * A);
        grad(0) += 1.0 * denom;
        grad(1) += beta * A * denom;
        grad(2) += (mu * A - mu * beta * B) * denom;
    }

    // Likelihood of non occurrence
    expint = exp(- beta * (end - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = end * A - C;
    lik -= eta * end + mu * ((double)n - A);
    grad(0) -= end;
    grad(1) -= ( (double)n - A );
    grad(2) -= mu * B;

    return Rcpp::List::create(Rcpp::Named("objective") = lik, Rcpp::Named("gradient") = grad);
}


////////////////////////////////////////////////
// SymmetricExponential

arma::vec SymmetricExponential::h( arma::vec x ) {
    return .5 * param(1) * param(2) * arma::exp( - param(2) * abs(x) );
}

arma::cx_vec SymmetricExponential::H( arma::vec xi ) {
    double beta2 = param(2) * param(2);
    arma::cx_vec zeta = arma::cx_vec( param(1) * beta2 / (beta2 + xi % xi), arma::zeros<arma::vec>(xi.n_elem) );
    return zeta;
}
