#include "model.h"
#include "utils.h"

// Model is given by its intensity :
// Intensity = param(0) + int h(t-u) dN(u)
// where h(t) = param(1) * param(2) * param(3) ^ param(2) * (param(3) + t) ^ (- param(2) - 1)
// Then:
// param(0): baseline intensity of the Hawkes process
// param(1): reproduction rate \in [0,1]
// param(2): shape parameter of the power law kernel
// param(3): scale parameter of the power law kernel

//////////////////////////////////////////////////////////////
// Methods for time- and frequency-domain excitation functions

arma::vec PowerLaw::h( arma::vec x ) {
    return param(1) * param(2) * exp(param(2) * log(param(3))) *
        arma::exp(-(param(2) + 1.0) * arma::log(param(3) + x));
}

arma::cx_vec PowerLaw::H( arma::vec xi ) {

    const double mu = param(1);
    const double theta = param(2);
    const double a = param(3);
    arma::vec xia = a * arma::abs(xi);

    arma::cx_vec y = mu * (1 + Etheta_imaginary(theta, xia));

    // Iterators
    arma::vec::iterator it_xi = xi.begin();
    arma::vec::iterator it_xi_end = xi.end();
    arma::cx_vec::iterator it_y = y.begin();

    // For xi < 0, take the conjugate
    for (; it_xi != it_xi_end; ++it_xi, ++it_y) {
        if (*it_xi == 0.0)
            *it_y = mu;
        else if (*it_xi < 0.0)
            *it_y = std::conj(*it_y);
    }

    return y;

}

arma::cx_mat PowerLaw::dH( arma::vec xi ) {

    // not functional

    const double mu = param(1);
    const double theta = param(2);
    const double a = param(3);
    arma::vec xia = a * arma::abs(xi);

    // We have H = mu * (1 + almostH)
    arma::cx_vec almostH = Etheta_imaginary(theta, xia);

    arma::cx_mat grad = arma::zeros<arma::cx_mat>( xi.n_elem, param.n_elem );
    grad.col(1) = 1.0 + almostH;
    grad.col(2) = mu * (i * xia % (1.0 + almostH)); // This is so wrong - finite-difference approximation instead?
    grad.col(3) = mu * (i * xi + (i * xi + theta / a) % almostH);

    //// TO FINISH: conjugate values for xi < 0
    // Iterators
    arma::vec::iterator it_xi = xi.begin();
    arma::vec::iterator it_xi_end = xi.end();
    arma::uword k = 0;

    arma::cx_rowvec deriv_at_0 = { 0.0, 1.0, 0.0, 0.0 };

    // For xi < 0, take the conjugate
    for (; it_xi != it_xi_end; ++it_xi, ++k) {
        if (*it_xi == 0.0)
            grad.row(k) = deriv_at_0;
        else if (*it_xi < 0.0)
            grad.row(k) = arma::conj(grad.row(k));
    }

    return grad;

}

////////////////////////////////////////////
// Methods for maximum likelihood estimation

double PowerLaw::loglik( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double theta = param(2);
    const double a = param(3);
    const double log_a = log(a);
    const double atheta = exp(theta * log_a);

    // Sum log \lambda*
    double part1 = log(eta);
    double hsum;
    for (arma::uword i = 1; i < n; i++) {
        hsum = arma::accu(arma::exp(-(theta+1.0) * arma::log(a + events(i) - events.subvec(0, i-1))));
        part1 += log(eta + mu * theta * atheta * hsum);
    }

    // Int \lambda* dt
    hsum = arma::accu(arma::exp(- theta * arma::log(a + end - events)));
    double part2 = eta * end + mu * (double)n - mu * atheta * hsum;

    return part1 - part2;
}

arma::vec PowerLaw::dloglik( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double theta = param(2);
    const double a = param(3);
    const double log_a = log(a);
    const double atheta = exp(theta * log_a);

    // Fill for i = 0
    arma::vec grad = { 1.0/eta, 0.0, 0.0, 0.0 };

    // Iterate on arrival times
    double hsum, inv_lambda;
    for (arma::uword i = 1; i < n; i++) {
        arma::vec diffs = a + events(i) - events.subvec(0, i-1);
        arma::vec log_diffs = arma::log(diffs);
        arma::vec diffsmthetam1 = arma::exp(-theta * log_diffs) / diffs;
        hsum = arma::accu(diffsmthetam1);
        inv_lambda = 1.0 / (eta + mu * theta * atheta * hsum);

        grad(0) += inv_lambda;
        grad(1) += inv_lambda * theta * atheta * hsum;
        grad(2) += inv_lambda * mu * atheta *
            arma::accu(
                diffsmthetam1 % (
                        1.0 + theta * (log_a - log_diffs)
                )
            );
        grad(3) += inv_lambda * mu * theta * (atheta / a) *
            arma::accu(
                (theta * (events(i) - events.subvec(0, i-1)) - a) %
                    diffsmthetam1 / diffs
            );
    }

    arma::vec log_diffs = arma::log(a + end - events);
    grad(0) -= end;
    grad(1) -= n - atheta * arma::accu(arma::exp(- theta * log_diffs));
    grad(2) += mu * atheta * arma::accu(
        (log_a - log_diffs) % arma::exp(-theta * log_diffs)
    );
    grad(3) += mu * theta * exp((theta - 1.0) * log_a) * arma::accu(
        (end - events) % arma::exp(-(theta + 1.0) * log_diffs)
    );

    return grad;
}

Rcpp::List PowerLaw::loglikngrad( const arma::vec& events, double end ) {

    // Constants
    const arma::uword n = events.n_elem;
    const double eta = param(0);
    const double mu = param(1);
    const double theta = param(2);
    const double a = param(3);
    const double log_a = log(a);
    const double atheta = exp(theta * log_a);

    // Fill for i = 0
    double lik = log(eta);
    arma::vec grad = { 1.0/eta, 0.0, 0.0, 0.0 };

    // Iterate on arrival times
    double hsum, inv_lambda;
    for (arma::uword i = 1; i < n; i++) {
        arma::vec diffs = a + events(i) - events.subvec(0, i-1);
        arma::vec log_diffs = arma::log(diffs);
        arma::vec diffsmthetam1 = arma::exp(-(theta+1.0) * log_diffs);
        hsum = arma::accu(diffsmthetam1);
        inv_lambda = 1.0 / (eta + mu * theta * atheta * hsum);

        lik += log(eta + mu * theta * atheta * hsum);

        grad(0) += inv_lambda;
        grad(1) += inv_lambda * theta * atheta * hsum;
        grad(2) += inv_lambda * mu * atheta *
            arma::accu(
                diffsmthetam1 % (
                        1.0 + theta * (log_a - log_diffs)
                )
            );
        grad(3) += inv_lambda * mu * theta * (atheta / a) *
            arma::accu(
                (theta * (events(i) - events.subvec(0, i-1)) - a) %
                    diffsmthetam1 / diffs
            );
    }

    // Compensator part
    arma::vec log_diffs = arma::log(a + end - events);
    hsum = arma::accu(arma::exp(- theta * log_diffs));

    lik -= eta * end + mu * (double)n - mu * atheta * hsum;

    grad(0) -= end;
    grad(1) -= n - atheta * arma::accu(arma::exp(- theta * log_diffs));
    grad(2) += mu * atheta * arma::accu(
        (log_a - log_diffs) % arma::exp(-theta * log_diffs)
    );
    grad(3) += mu * theta * exp((theta - 1.0) * log_a) * arma::accu(
        (end - events) % arma::exp(-(theta + 1.0) * log_diffs)
    );

    return Rcpp::List::create(Rcpp::Named("objective") = lik, Rcpp::Named("gradient") = grad);
}
