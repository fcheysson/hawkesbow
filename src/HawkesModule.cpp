#include "model.hpp"
#include "data.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

RCPP_MODULE(HawkesModule) {
    using namespace Rcpp;

    class_<Model>("Model")
        .default_constructor() // This exposes the default constructor
    ;

    class_<ContinuousData>("ContinuousData")
        .default_constructor()
        .constructor<arma::vec,double>()
        .constructor<arma::vec,double,double>()
        .property("events", &ContinuousData::getEvents, &ContinuousData::setEvents)
        .property("timeBegin", &ContinuousData::getTimeBegin, &ContinuousData::setTimeBegin)
        .property("timeEnd", &ContinuousData::getTimeEnd, &ContinuousData::setTimeEnd)
        .property("timeRange", &ContinuousData::getTimeRange, &ContinuousData::setTimeRange)
    ;
}
