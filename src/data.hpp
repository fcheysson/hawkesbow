#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

class Data {

protected:
    double timeBegin, timeEnd;     // Time bounds

public:
    Data() {};
    Data( double timeBegin, double timeEnd ) : timeBegin(timeBegin), timeEnd(timeEnd) {};

    virtual ~Data() {};

    double    getTimeBegin()  { return timeBegin; };
    double    getTimeEnd()    { return timeEnd; };
    arma::vec getTimeRange()  { return { timeBegin, timeEnd }; };

    // Virtual methods to call a derived class member from a pointer to base class object
    virtual arma::Col<unsigned int> getCounts()     { return arma::zeros<arma::Col<unsigned int>>(1); };
    virtual double                  getBinsize()    { return 0.0; };
    virtual arma::vec               getEvents()     { return arma::zeros<arma::vec>(1); };

};

//' @export DiscreteData
class DiscreteData: public Data {

private:
    arma::Col<unsigned int> counts;                 // Vector of counts
    double                  binsize;

public:
    // // Default constructor           // Removed since there are no setter methods
    // DiscreteData() {};

    // Constructor from vector of occurrences and binsize
    DiscreteData( arma::Col<unsigned int> counts, double binsize ) :
        Data( 0, binsize*counts.n_elem),
        counts( counts ),
        binsize( binsize )
        {};

    // Constructor from vector of occurrences, binsize and timeBegin
    DiscreteData( arma::Col<unsigned int> counts, double binsize, double timeBegin ) :
        Data( timeBegin, timeBegin + binsize*counts.n_elem ),
        counts( counts ),
        binsize( binsize )
        {};

    // Getters
    arma::Col<unsigned int> getCounts()     { return counts; };
    double                  getBinsize()    { return binsize; };

};

//' @export ContinuousData
class ContinuousData: public Data {

private:
    arma::vec events;               // Vector of occurrences

public:
    // // Default constructor
    // ContinuousData() {};

    // Constructor from vector of occurrences and time upper bound
    ContinuousData( arma::vec events, double timeEnd ) :
        Data( 0, timeEnd ),
        events( events )
        {};

    // Constructor from vector of occurrences and time bounds
    ContinuousData( arma::vec events, double timeBegin, double timeEnd ) :
        Data( timeBegin, timeEnd ),
        events( events )
        {};

    // Getters
    arma::vec   getEvents()     { return events; };

    // Discretize from length or binsize
    DiscreteData binl( unsigned int length );
    DiscreteData bins( double binsize );

};

