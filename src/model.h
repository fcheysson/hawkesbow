#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export Model
class Model {

protected:
    arma::vec param;
    double binsize;

public:
    Model() {};
    Model( arma::vec param ) : param(param), binsize(1.0) {};
    Model( double binsize ) : param(arma::zeros<arma::vec>(1)), binsize(binsize) {};
    Model( arma::vec param, double binsize ) : param(param), binsize(binsize) {};
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

    // Likelihood estimation methods
    virtual double loglik( const arma::vec& events, double end ) { return 0.0; };
    virtual arma::vec dloglik( const arma::vec& events, double end ) { return arma::zeros<arma::vec>(param.n_elem); };
    virtual arma::mat ddloglik( const arma::vec& events, double end ) { return arma::zeros<arma::mat>(param.n_elem, param.n_elem); };
    virtual Rcpp::List loglikngrad( const arma::vec& events, double end ) { return Rcpp::List::create(); };

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

    // Whittle likelihood estimation methods
    double whittle( arma::vec& I, int trunc );

    // Property for param
    void setParam( arma::vec param_ ) { param = param_; };
    arma::vec getParam() { return param; };

    // Property for binsize
    void setBinsize( double binsize_ ) { binsize = binsize_; };
    double getBinsize() { return binsize; };

};

//' @export Exponential
class Exponential: public Model {

private:


public:
    Exponential() : Model(arma::vec({1.0, 0.5, 1.0})) {};
    Exponential( arma::vec param ) : Model(param) {};
    Exponential( double binsize ) : Model(arma::vec({1.0, 0.5, 1.0}), binsize) {};
    Exponential( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );
    arma::cx_mat dH( arma::vec xi );
    arma::cx_cube ddH( arma::vec xi );

    // Likelihood methods
    double loglik( const arma::vec& events, double end );
    arma::vec dloglik( const arma::vec& events, double end );
    arma::mat ddloglik( const arma::vec& events, double end );
    Rcpp::List loglikngrad( const arma::vec& events, double end );

};

//' @export SymmetricExponential
class SymmetricExponential: public Model {

private:


public:
    SymmetricExponential() : Model(arma::vec({1.0, 0.5, 1.0})) {};
    SymmetricExponential( arma::vec param ) : Model(param) {};
    SymmetricExponential( double binsize ) : Model(arma::vec({1.0, 0.5, 1.0}), binsize) {};
    SymmetricExponential( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );
    // arma::cx_mat dH( arma::vec xi );
    // arma::cx_cube ddH( arma::vec xi );

    // Likelihood methods
    // double loglik( const arma::vec& events, double end );
    // arma::vec dloglik( const arma::vec& events, double end );
    // arma::mat ddloglik( const arma::vec& events, double end );
    // Rcpp::List loglikngrad( const arma::vec& events, double end );

};

//' @export PowerLaw
class PowerLaw: public Model {

private:


public:
    PowerLaw() : Model(arma::vec({1.0, 0.5, 3.0, 1.0})) {};
    PowerLaw( arma::vec param ) : Model(param) {};
    PowerLaw( double binsize ) : Model(arma::vec({1.0, 0.5, 3.0, 1.0}), binsize) {};
    PowerLaw( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );
    arma::cx_mat dH( arma::vec xi );

    // Likelihood methods
    double loglik( const arma::vec& events, double end );
    arma::vec dloglik( const arma::vec& events, double end );
    Rcpp::List loglikngrad( const arma::vec& events, double end );

};

//' @export Pareto3
class Pareto3: public Model {

private:


public:
    Pareto3() : Model(arma::vec({1.0, 0.5, 1.0})) {};
    Pareto3( arma::vec param ) : Model(param) {};
    Pareto3( double binsize ) : Model(arma::vec({1.0, 0.5, 1.0}), binsize) {};
    Pareto3( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );

};

//' @export Pareto2
class Pareto2: public Model {

private:


public:
    Pareto2() : Model(arma::vec({1.0, 0.5, 1.0})) {};
    Pareto2( arma::vec param ) : Model(param) {};
    Pareto2( double binsize ) : Model(arma::vec({1.0, 0.5, 1.0}), binsize) {};
    Pareto2( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );

};

//' @export Pareto1
class Pareto1: public Model {

private:


public:
    Pareto1() : Model(arma::vec({1.0, 0.5, 1.0})) {};
    Pareto1( arma::vec param ) : Model(param) {};
    Pareto1( double binsize ) : Model(arma::vec({1.0, 0.5, 1.0}), binsize) {};
    Pareto1( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );

};

//' @export Gaussian
class Gaussian: public Model {

private:


public:
    Gaussian() : Model(arma::vec({1.0, 0.5, 0.0, 1.0})) {};
    Gaussian( arma::vec param ) : Model(param) {};
    Gaussian( double binsize ) : Model(arma::vec({1.0, 0.5, 0.0, 1.0}), binsize) {};
    Gaussian( arma::vec param, double binsize ) : Model(param, binsize) {};

    // Methods for time- and frequency-domain excitation functions
    arma::vec h( arma::vec x );
    arma::cx_vec H( arma::vec xi );

};
