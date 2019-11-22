#include "model.hpp"
#include "data.hpp"

// // For overloaded methods: help disambiguate which one to use
// // Here: select method if argument is an integer
// bool get_int_valid(SEXP* args, int nargs){
//     if( nargs != 1 ) return false ;
//     if( TYPEOF(args[0]) != INTSXP ) return false ;
//     return ( LENGTH(args[0]) == 1 ) ;
// }

RCPP_EXPOSED_CLASS(Data);
RCPP_EXPOSED_CLASS(DiscreteData);

RCPP_MODULE(HawkesModule) {
    using namespace Rcpp;

    class_<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("G", &Model::G)
        .method("dG", &Model::dG)
        .method("ddG", &Model::ddG)
        .method("f", &Model::f)
        .method("df", &Model::df)
        .method("ddf", &Model::ddf)
        .method("f1", &Model::f1)
        .method("df1", &Model::df1)
        .method("ddf1", &Model::ddf1)
        .method("whittle", &Model::whittle)
        .method("attach", &Model::attach)
        .method("detach", &Model::detach)
        .property("param", &Model::getParam, &Model::setParam)
        .property("data", &Model::getData)
    ;

    class_<Exponential>("Exponential")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("mean", &Exponential::mean)
        .method("dmean", &Exponential::dmean)
        .method("ddmean", &Exponential::ddmean)
        .method("h", &Exponential::h)
        .method("H", &Exponential::H)
        .method("dH", &Exponential::dH)
        .method("ddH", &Exponential::ddH)
        .method("loglik", &Exponential::loglik)
        .method("dloglik", &Exponential::dloglik)
        .method("ddloglik", &Exponential::ddloglik)
        .method("loglikngrad", &Exponential::loglikngrad)
    ;

    class_<Data>("Data")
        // .default_constructor()
        .constructor<double,double>()
        .property("timeBegin", &DiscreteData::getTimeBegin)
        .property("timeEnd", &DiscreteData::getTimeEnd)
        .property("timeRange", &DiscreteData::getTimeRange)
    ;

    class_<DiscreteData>("DiscreteData")
        .derives<Data>("Data")
        // .default_constructor()
        .constructor<arma::Col<unsigned int>,double>()
        .constructor<arma::Col<unsigned int>,double,double>()
        .property("counts", &DiscreteData::getCounts)
        .property("binsize", &DiscreteData::getBinsize)
    ;

    // // Help the compiler disambiguate things
    // DiscreteData (ContinuousData::*toDiscrete_length)(unsigned int) = &ContinuousData::toDiscrete;
    // DiscreteData (ContinuousData::*toDiscrete_binsize)(double) = &ContinuousData::toDiscrete;

    class_<ContinuousData>("ContinuousData")
        .derives<Data>("Data")
        // .default_constructor()
        .constructor<arma::vec,double>()
        .constructor<arma::vec,double,double>()
        .property("events", &ContinuousData::getEvents)
        .method("binl", &ContinuousData::binl)
        .method("bins", &ContinuousData::bins)
        // .method("toDiscrete", ( DiscreteData (ContinuousData::*)(unsigned int) )(&ContinuousData::toDiscrete), "ContinuousData", &get_int_valid )
        // .method("toDiscrete", ( DiscreteData (ContinuousData::*)(double) )      (&ContinuousData::toDiscrete))
    ;
}
