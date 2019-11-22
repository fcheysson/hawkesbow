#include "model.hpp"
#include "utils.hpp"

// Likelihood methods
double Exponential::loglik() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
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
    double part2 = eta * T + mu * ((double)n - exp(-beta * (T - events(n-1))) * (1.0 + A));

    return part1 - part2;
};

arma::vec Exponential::dloglik() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
    const double eta = param(0);
    const double alpha = param(1)*param(2);
    const double beta = param(2);
    const double inv_beta = 1.0/beta;

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
        denom = 1.0 / (eta + alpha * A);
        grad(0) += 1.0 * denom;
        grad(1) += A * denom;
        grad(2) -= alpha * B * denom;
    }

    expint = exp(- beta * (T - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = T * A - C;
    grad(0) -= T;
    grad(1) -= inv_beta * ( (double)n - A );
    grad(2) += alpha * inv_beta * ( inv_beta * ( (double)n - A ) - B );

    return grad;
};

arma::mat Exponential::ddloglik() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
    const double eta = param(0);
    const double alpha = param(1)*param(2);
    const double beta = param(2);
    const double inv_beta = 1.0/beta;

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
        denom = 1.0 / (eta + alpha * A);
        denom2 = denom * denom;
        hess(0,0) -= 1.0 * denom2;
        hess(0,1) -= A * denom2;
        hess(0,2) += alpha * B * denom2;
        hess(1,1) -= A * A * denom2;
        hess(1,2) -= eta * B * denom2;
        hess(2,2) += alpha * denom * (D - alpha * denom * B * B);
    }

    expint = exp(- beta * (T - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = T * A - C;
    E = expint * (events(n-1) * events(n-1) + E);
    D = T * T * A + E - 2.0 * T * C;
    hess(1,2) += inv_beta * (inv_beta * ((double)n - A) - B);
    hess(2,2) += alpha * inv_beta * (D + 2.0 * inv_beta * (B - inv_beta * ((double)n  - A)));

    // Symm
    hess(1,0) = hess(0,1);
    hess(2,0) = hess(0,2);
    hess(2,1) = hess(1,2);

    return hess;
};

Rcpp::List Exponential::loglikngrad() {

    // Constants
    const arma::vec events = data->getEvents();
    const arma::uword n = events.n_elem;
    const double T = data->getTimeEnd() - data->getTimeBegin();
    const double eta = param(0);
    const double alpha = param(1)*param(2);
    const double beta = param(2);
    const double inv_beta = 1.0/beta;

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
        denom = 1.0 / (eta + alpha * A);
        lik += log(eta + alpha * A);
        grad(0) += 1.0 * denom;
        grad(1) += A * denom;
        grad(2) -= alpha * B * denom;
    }

    // Likelihood of non occurrence
    expint = exp(- beta * (T - events(n-1)));
    A = expint * (1.0 + A);
    C = expint * (events(n-1) + C);
    B = T * A - C;
    lik -= eta * T + (alpha / beta) * ((double)n - A);
    grad(0) -= T;
    grad(1) -= inv_beta * ( (double)n - A );
    grad(2) += alpha * inv_beta * ( inv_beta * ( (double)n - A ) - B );

    return Rcpp::List::create(Rcpp::Named("objective") = lik, Rcpp::Named("gradient") = grad);
};
