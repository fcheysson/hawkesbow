#pragma once
#include <RcppArmadillo.h>


//' @export DiscreteData
class DiscreteData {
    arma::vec counts;               // Vector of counts
    double binsize;
    int n_elem;
    double timeBegin, timeEnd;      // Time bounds

    // Default constructor
    DiscreteData() {};

    // Constructor from vector of occurrences and binsize
    DiscreteData( arma::vec counts, int binsize ) :
        counts( counts ), binsize( binsize ), n_elem( counts.n_elem ), timeBegin( 0 ), timeEnd( binsize*n_elem ) {};
};

//' @export ContinuousData
class ContinuousData {

private:
    arma::vec events;               // Vector of occurrences
    double timeBegin, timeEnd;      // Time bounds

public:
    // Default constructor
    ContinuousData() {};

    // Constructor from vector of occurrences and time upper bound
    ContinuousData( arma::vec events, double timeEnd ) :
        events( events ), timeBegin( 0 ), timeEnd( timeEnd ) {};

    // Constructor from vector of occurrences and time bounds
    ContinuousData( arma::vec events, double timeBegin, double timeEnd ) :
        events( events ), timeBegin( timeBegin ), timeEnd( timeEnd ) {};

    // Getters and setters
    void setEvents( arma::vec events_ ) { events = events_; };
    arma::vec getEvents() { return events; };

    void setTimeBegin( double timeBegin_ ) { timeBegin = timeBegin_; };
    double getTimeBegin() { return timeBegin; };

    void setTimeEnd( double timeEnd_ ) { timeEnd = timeEnd_; };
    double getTimeEnd() { return timeEnd; };

    void setTimeRange( arma::vec timeRange ) { timeBegin = timeRange(0); timeEnd = timeRange(1); };
    arma::vec getTimeRange() { return { timeBegin, timeEnd }; };

    // Discretize from length
    struct DiscreteData toDiscrete( int length );

};

