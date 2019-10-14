#include "data.hpp"
#include "utils.hpp"

// Discretize continuous data to discrete data
struct DiscreteData ContinuousData::toDiscrete( int length ) {
    arma::vec bounds = arma::linspace( timeBegin, timeEnd, length + 1 );
    arma::vec counts = arma::zeros( length );
    arma::uword index = 1;
    for (arma::uword )
};
