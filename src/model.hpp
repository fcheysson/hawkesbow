#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' @export Model
class Model {

protected:
    arma::vec param;

public:
    virtual ~Model() {};

    // Methods for long term mean and its derivatives
    double mean();
    arma::vec dmean();
    arma::mat ddmean();

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
    arma::mat df( arma::vec xi );
    arma::cube ddf( arma::vec xi );

    arma::vec f1( arma::vec xi, int trunc );  // f1(w) = sum_{k=-trunc}^{+trunc} f(w + 2*k*pi)
    arma::mat df1( arma::vec xi, int trunc );
    arma::cube ddf1( arma::vec xi, int trunc );

    // Property for param
    void setParam( arma::vec param_ ) { param = param_; };
    arma::vec getParam() { return param; };

};

//' @export Exponential
class Exponential: public Model {

private:


public:
    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );
    arma::cx_mat dH( arma::vec xi );
    arma::cx_cube ddH( arma::vec xi );

};