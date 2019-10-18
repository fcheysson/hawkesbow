#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' @export Model
class Model {

private:
    arma::vec param;

public:
    virtual ~Model() {};

    // Virtual methods for time- and frequency-domain excitation functions
    virtual arma::vec h( arma::vec x ) { return arma::zeros<arma::vec>(x.n_elem); };
    virtual arma::cx_vec H( arma::vec xi ) { return arma::zeros<arma::cx_vec>(xi.n_elem); };
    virtual arma::cx_mat dH( arma::vec xi ) { return arma::zeros<arma::cx_mat>(param.n_elem, xi.n_elem); };
    virtual arma::cx_cube ddH( arma::vec xi ) { return arma::zeros<arma::cx_cube>(param.n_elem, param.n_elem, xi.n_elem); };

    // Methods for continuous- and discretized-time spectral densities
    arma::vec G( arma::vec xi );    // G(w) = |1-H(w)|^{-2}
    arma::mat dG( arma::vec xi );
    arma::cube ddG( arma::vec xi );

    arma::vec f( arma::vec xi );    // f(w) = m * binsize * sinc²(w/2) * G(w/binsize)
    arma::mat gradf( arma::vec xi );
    arma::cube hessf( arma::vec xi );

    arma::vec f_aliasing( arma::vec xi, int trunc );  // f_aliasing(w) = sum_{k=-trunc}^{+trunc} f(w + 2*k*pi)
    arma::mat gradf_aliasing( arma::vec xi, int trunc );
    arma::cube hessf_aliasing( arma::vec xi, int trunc );

};

// //' @export Exponential
// class Exponential: public Model {
//
// private:
//
//
// public:
//
//
// };
