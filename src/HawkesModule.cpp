#include "model.hpp"
#include "data.hpp"

// For overloaded methods: help disambiguate which one to use
// Here: select method if argument is an integer
bool get_int_valid(SEXP* args, int nargs){
    if( nargs != 1 ) return false ;
    if( TYPEOF(args[0]) != INTSXP ) return false ;
    return ( LENGTH(args[0]) == 1 ) ;
}

RCPP_EXPOSED_CLASS(DiscreteData);

RCPP_MODULE(HawkesModule) {
    using namespace Rcpp;

    class_<Model>("Model")
        .default_constructor() // This exposes the default constructor
    ;

    class_<DiscreteData>("DiscreteData")
        // .default_constructor()
        .constructor<arma::Col<unsigned int>,double>()
        .constructor<arma::Col<unsigned int>,double,double>()
        .property("counts", &DiscreteData::getCounts)
        .property("binsize", &DiscreteData::getBinsize)
        .property("timeBegin", &DiscreteData::getTimeBegin)
        .property("timeEnd", &DiscreteData::getTimeEnd)
        .property("timeRange", &DiscreteData::getTimeRange)
    ;

    // // Help the compiler disambiguate things
    // DiscreteData (ContinuousData::*toDiscrete_length)(unsigned int) = &ContinuousData::toDiscrete;
    // DiscreteData (ContinuousData::*toDiscrete_binsize)(double) = &ContinuousData::toDiscrete;

    class_<ContinuousData>("ContinuousData")
        // .default_constructor()
        .constructor<arma::vec,double>()
        .constructor<arma::vec,double,double>()
        .property("events", &ContinuousData::getEvents)
        .property("timeBegin", &ContinuousData::getTimeBegin)
        .property("timeEnd", &ContinuousData::getTimeEnd)
        .property("timeRange", &ContinuousData::getTimeRange)
        .method("toDiscrete", ( DiscreteData (ContinuousData::*)(unsigned int) )(&ContinuousData::toDiscrete), "ContinuousData", &get_int_valid )
        .method("toDiscrete", ( DiscreteData (ContinuousData::*)(double) )      (&ContinuousData::toDiscrete))
    ;
}
