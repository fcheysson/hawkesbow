#include "model.h"

RCPP_MODULE(HawkesModule) {
    using namespace Rcpp;

    class_<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("mean", &Model::mean)
        .method("dmean", &Model::dmean)
        .method("ddmean", &Model::ddmean)
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
        .property("param", &Model::getParam, &Model::setParam)
        .property("binsize", &Model::getBinsize, &Model::setBinsize)
    ;

    class_<Exponential>("Exponential")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &Exponential::h)
        .method("H", &Exponential::H)
        .method("dH", &Exponential::dH)
        .method("ddH", &Exponential::ddH)
        .method("loglik", &Exponential::loglik)
        .method("dloglik", &Exponential::dloglik)
        .method("ddloglik", &Exponential::ddloglik)
        .method("loglikngrad", &Exponential::loglikngrad)
    ;

    class_<SymmetricExponential>("SymmetricExponential")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &SymmetricExponential::h)
        .method("H", &SymmetricExponential::H)
    ;

    class_<PowerLaw>("PowerLaw")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &PowerLaw::h)
        .method("H", &PowerLaw::H)
        // .method("dH", &PowerLaw::dH)
        .method("loglik", &PowerLaw::loglik)
        .method("dloglik", &PowerLaw::dloglik)
        .method("loglikngrad", &PowerLaw::loglikngrad)
    ;

    class_<Pareto3>("Pareto3")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &Pareto3::h)
        .method("H", &Pareto3::H)
    ;

    class_<Pareto2>("Pareto2")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &Pareto2::h)
        .method("H", &Pareto2::H)
    ;

    class_<Pareto1>("Pareto1")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &Pareto1::h)
        .method("H", &Pareto1::H)
    ;

    class_<Gaussian>("Gaussian")
        .derives<Model>("Model")
        .default_constructor() // This exposes the default constructor
        .method("h", &Gaussian::h)
        .method("H", &Gaussian::H)
    ;
}
