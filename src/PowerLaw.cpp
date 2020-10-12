#include "model.hpp"
#include "utils.hpp"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) * param(2) * param(3) ^ param(2) * (param(3) + t) ^ (- param(2) - 1)
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate \in [0,1]
// param(2): shape parameter of the power law kernel
// param(3): scale parameter of the power law kernel

////////////////////////////////////////////////
// Methods for long term mean and its derivative

double PowerLaw::mean() {
    return param(0) / ( 1.0 - param(1) );
}

arma::vec PowerLaw::dmean() {
    double denom = 1.0 / ( 1.0 - param(1) );
    arma::vec grad = { denom, param(0) * denom * denom, 0 };
    return grad;
}

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec PowerLaw::h( arma::vec x ) {
    return param(1) * param(2) * exp(param(2) * log(param(3))) *
        arma::exp(-(param(2) + 1) * arma::log(param(3) + x));
}

arma::cx_vec PowerLaw::H( arma::vec xi ) {              //todo
    arma::vec factor = param(1) * param(2) / ( param(2)*param(2) + xi%xi );
    arma::cx_vec zeta = arma::cx_vec( factor * param(2), - factor % xi );
    return zeta;
}

////////////////////////////////////////////
// Methods for maximum likelihood estimation

double PowerLaw::loglik() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
    const double eta = param(0);
    const double mu = param(1);
    const double theta = param(2);
    const double a = param(3);
    const double atheta = exp(theta * log(a));

    // Sum log \lambda*
    double part1 = log(eta);
    double hsum;
    for (arma::uword i = 1; i < n; i++) {
        hsum = arma::accu(arma::exp(-(theta+1) * arma::log(a + events.subvec(0, i-1))));
        part1 += log(eta + mu * theta * atheta * hsum);
    }

    // Int \lambda* dt
    hsum = arma::accu(arma::exp(theta * (log(a) - arma::log(a + T - events))));
    double part2 = eta * T + mu * (double)n - mu * hsum;

    return part1 - part2;
};

arma::vec PowerLaw::dloglik() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
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

    expint = exp(- beta * (T - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = T * A - C;
    grad(0) -= T;
    grad(1) -= ( (double)n - A );
    grad(2) -= mu * B;

    return grad;
};

Rcpp::List PowerLaw::loglikngrad() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
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
    expint = exp(- beta * (T - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = T * A - C;
    lik -= eta * T + mu * ((double)n - A);
    grad(0) -= T;
    grad(1) -= ( (double)n - A );
    grad(2) -= mu * B;

    return Rcpp::List::create(Rcpp::Named("objective") = lik, Rcpp::Named("gradient") = grad);
};
