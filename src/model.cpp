#include "model.h"
#include "utils.h"

////////////////////////////////////////////////
// Methods for long term mean and its derivative

double Model::mean() {
    return param(0) / ( 1.0 - param(1) );
}

arma::vec Model::dmean() {
    double denom = 1.0 / ( 1.0 - param(1) );
    arma::vec grad = { denom, param(0) * denom * denom, 0 };
    return grad;
}

arma::mat Model::ddmean() {
    double denom = 1.0 / ( 1.0 - param(1) );
    double denom2 = denom * denom;
    arma::mat hess = {  {   0.0,                  denom2, 0.0},
                        {denom2, 2*param(0)*denom2*denom, 0.0},
                        {   0.0,                     0.0, 0.0} };
    return hess;
}

//////////////////////////////////////////////////////////////////
// Methods for continuous- and discretized-time spectral densities

// G(w) = |1-H(w)|^{-2}

arma::vec Model::G( arma::vec xi ) {
    arma::cx_vec temp = arma::cx_double(1.0, 0.0) - H( xi );
    return 1.0 / arma::conv_to<arma::vec>::from( temp % arma::conj(temp) );
}

arma::mat Model::dG( arma::vec xi ) {
    arma::mat grad( xi.n_elem, param.n_elem );
    arma::cx_mat dHxi = dH(xi);
    arma::cx_vec term1 = 1.0 - arma::conj(H(xi));
    arma::vec Gxi = G(xi);
    arma::vec Gxi2 = Gxi % Gxi;
    for (arma::uword k = 0; k < param.n_elem; k++) {
        grad.col(k) = 2.0 * Gxi2 % arma::real( term1 % dHxi.col(k) );
    }
    return grad;
}

arma::cube Model::ddG( arma::vec xi ) {
    arma::cube hess( param.n_elem, param.n_elem, xi.n_elem );
    arma::vec Gxi = G(xi);
    arma::vec Gxi2 = Gxi % Gxi;
    arma::cx_vec term0 = 1.0 - arma::conj(H(xi));
    arma::cx_mat dHxi = dH(xi);
    arma::cx_cube ddHxi = ddH(xi);
    arma::vec term1(xi.n_elem);
    arma::vec term2(xi.n_elem);
    arma::cx_vec tube(xi.n_elem);
    for (arma::uword i = 0; i < param.n_elem; i++) {
        for (arma::uword j = 0; j < param.n_elem; j++) {
            tube = ddHxi(arma::span(i),arma::span(j), arma::span::all);
            term1 = arma::real( term0 % tube - dHxi.col(i) % arma::conj(dHxi.col(j)) );
            term2 = arma::real( term0 % dHxi.col(i) ) % arma::real( term0 % dHxi.col(j) );
            hess.tube(i, j) = 2 * Gxi2 % (term1 + 4 * Gxi % term2);
        }
    }
    return hess;
}

// f(w) = m * binsize * sinc²(w/2) * G(w/binsize)

arma::vec Model::f( arma::vec xi ) {
    arma::vec term1 = _sinc( .5 * xi );
    return mean() * binsize * term1 % term1 % G(xi / binsize);
}

arma::mat Model::df( arma::vec xi ) {
    arma::mat grad( xi.n_elem, param.n_elem );

    arma::vec term0 = _sinc( .5 * xi );
    arma::vec term1 = term0 % term0;

    double m = mean();
    arma::vec dm = dmean();
    arma::vec Gxi = G( xi / binsize );
    arma::mat dGxi = dG( xi / binsize );

    for (arma::uword k = 0; k < param.n_elem; k++) {
        grad.col(k) = binsize * term1 % ( dm(k) * Gxi + m * dGxi.col(k) );
    }
    return grad;
}

arma::cube Model::ddf( arma::vec xi ) {
    arma::cube hess( param.n_elem, param.n_elem, xi.n_elem );

    arma::vec term0 = _sinc( .5 * xi );
    arma::vec term1 = term0 % term0;

    double m = mean();
    arma::vec dm = dmean();
    arma::mat ddm = ddmean();
    arma::vec Gxi = G( xi / binsize );
    arma::mat dGxi = dG( xi / binsize );
    arma::cube ddGxi = ddG( xi / binsize );

    arma::vec tube(xi.n_elem);
    for (arma::uword i = 0; i < param.n_elem; i++) {
        for (arma::uword j = 0; j < param.n_elem; j++) {
            tube = ddGxi(arma::span(i),arma::span(j), arma::span::all);
            hess.tube(i, j) = binsize * term1 % ( ddm(i, j) * Gxi + dm(i) * dGxi.col(j) + dm(j) * dGxi.col(i) + m * tube );
        }
    }
    return hess;
}

// f1(w) = sum_{k=-trunc}^{+trunc} f(w + 2*k*pi)

arma::vec Model::f1( arma::vec xi, int trunc ) {
    arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
    arma::vec y(xi.n_elem);
    for (arma::uword k = 0; k < xi.n_elem; k++) {
        y(k) = arma::sum( f(xi(k) + omega) );
    }
    return y;
}

arma::mat Model::df1( arma::vec xi, int trunc ) {
    arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
    arma::mat y(xi.n_elem, param.n_elem);
    for (arma::uword k = 0; k < xi.n_elem; k++) {
        y.row(k) = arma::sum( df(xi(k) + omega), 0 );
    }
    return y;
}

arma::cube Model::ddf1( arma::vec xi, int trunc ) {
    arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
    arma::cube y(param.n_elem, param.n_elem, xi.n_elem);
    for (arma::uword k = 0; k < xi.n_elem; k++) {
        y.slice(k) = arma::sum( ddf(xi(k) + omega), 2 );
    }
    return y;
}

// Whittle likelihood estimation methods
double Model::whittle( arma::vec& periodogram, int trunc ) {
    arma::uword n = periodogram.n_elem + 1;     // + 1 because we removed element 0 in 'whittle.R'
    arma::uword n2 = n / 2;
    arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(1, n2) / (double)n;
    arma::vec spectrum = f1( omega, trunc );
    return arma::sum( arma::log(spectrum) + periodogram.subvec(0, n2-1) / spectrum );
}
