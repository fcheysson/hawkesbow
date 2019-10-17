#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' @export DiscreteData
class DiscreteData {

private:
    arma::Col<unsigned int> counts;                 // Vector of counts
    double                  binsize;
    double                  timeBegin, timeEnd;     // Time bounds

public:
    // // Default constructor           // Removed since there are no setter methods
    // DiscreteData() {};

    // Constructor from vector of occurrences and binsize
    DiscreteData( arma::Col<unsigned int> counts, double binsize ) :
        counts( counts ),
        binsize( binsize ),
        timeBegin( 0 ),
        timeEnd( binsize*counts.n_elem )
        {};

    // Constructor from vector of occurrences, binsize and timeBegin
    DiscreteData( arma::Col<unsigned int> counts, double binsize, double timeBegin ) :
        counts( counts ),
        binsize( binsize ),
        timeBegin( timeBegin ),
        timeEnd( timeBegin + binsize*counts.n_elem )
        {};

    // Getters
    arma::Col<unsigned int> getCounts()     { return counts; };
    double                  getBinsize()    { return binsize; };
    double                  getTimeBegin()  { return timeBegin; };
    double                  getTimeEnd()    { return timeEnd; };
    arma::vec               getTimeRange()  { return { timeBegin, timeEnd }; };

};

//' @export ContinuousData
class ContinuousData {

private:
    arma::vec   events;               // Vector of occurrences
    double      timeBegin, timeEnd;      // Time bounds

public:
    // // Default constructor
    // ContinuousData() {};

    // Constructor from vector of occurrences and time upper bound
    ContinuousData( arma::vec events, double timeEnd ) :
        events( events ),
        timeBegin( 0 ),
        timeEnd( timeEnd )
        {};

    // Constructor from vector of occurrences and time bounds
    ContinuousData( arma::vec events, double timeBegin, double timeEnd ) :
        events( events ),
        timeBegin( timeBegin ),
        timeEnd( timeEnd )
        {};

    // Getters
    arma::vec   getEvents()     { return events; };
    double      getTimeBegin()  { return timeBegin; };
    double      getTimeEnd()    { return timeEnd; };
    arma::vec   getTimeRange()  { return { timeBegin, timeEnd }; };

    // Discretize from length or binsize
    DiscreteData toDiscrete( unsigned int length );
    DiscreteData toDiscrete( double binsize );

};

