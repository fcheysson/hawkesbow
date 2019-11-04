#include "data.hpp"

// Discretize continuous data to discrete data
DiscreteData ContinuousData::binl( unsigned int length ) {

    // Bounds and counts of bins
    arma::vec               bounds = arma::linspace( timeBegin, timeEnd, length + 1 );
    arma::Col<unsigned int> counts = arma::zeros<arma::Col<unsigned int>>( length );

    // Associated iterators
    arma::vec::iterator                 it_bounds = bounds.begin();
    arma::Col<unsigned int>::iterator   it_counts = counts.begin();

    // Get binsize and prepare for loop
    double binsize = *(it_bounds + 1) - *it_bounds;
    ++it_bounds;    // First element is zero, start at the end of first bin instead

    // Iterators for events member
    arma::vec::iterator it_events = events.begin();
    arma::vec::iterator it_events_end = events.end();

    // Loop on events
    // Complexity = O(events.n_elem + length)
    while ( it_events != it_events_end ) {
        if (*it_events < *it_bounds) {
            *it_counts += 1;
            ++it_events;
        }
        else {
            ++it_bounds;
            ++it_counts;
        }
    }

    DiscreteData discreteData( counts, binsize, timeBegin );

    return discreteData;
}

// Discretize continuous data to discrete data
DiscreteData ContinuousData::bins( double binsize ) {

    // Bounds and counts of bins
    arma::vec               bounds = arma::regspace( timeBegin, binsize, timeEnd );
    arma::Col<unsigned int> counts = arma::zeros<arma::Col<unsigned int>>( bounds.n_elem - 1 );

    // Associated iterators
    arma::vec::iterator                 it_bounds = bounds.begin();
    arma::Col<unsigned int>::iterator   it_counts = counts.begin();

    // Prepare for loop
    ++it_bounds;    // First element is zero, start at the end of first bin instead

    // Iterators for events member
    arma::vec::iterator it_events = events.begin();
    arma::vec::iterator it_events_end = events.end();

    // Loop on events
    while ( it_events != it_events_end ) {
        if (*it_events < *it_bounds) {
            *it_counts += 1;
            ++it_events;
        }
        else {
            ++it_bounds;
            ++it_counts;
        }
    }

    DiscreteData discreteData( counts, binsize, timeBegin );

    return discreteData;
}
