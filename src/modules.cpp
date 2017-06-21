#include "global.h"
#include "Component.h"
#include "Mixture.h"

RCPP_MODULE(RcppComponent){
  using namespace Rcpp;

  // exposing the class in R as "Component"
  class_<Component>("Component")

  // constructor
  .constructor<int>("Create new component")
  .constructor<int, arma::mat>("Create new component and initialise")

  .field("N", &Component::N)
  .field("D", &Component::D)
  .field("nu", &Component::nu)
  .field("kappa", &Component::kappa)
  .field("m", &Component::m)
  .field("S", &Component::S)
  .field("L", &Component::L)

  .method("is_empty", &Component::is_empty)
  .method("add_sample", &Component::add_sample)
  .method("rm_sample", &Component::rm_sample)
  .method("marginal_loglik", &Component::marginal_loglik)
  .method("posterior_predictive", &Component::posterior_predictive)
  ;
}

RCPP_MODULE(RcppMixture){
  using namespace Rcpp;

  // exposing the class in R as "Component"
  class_<Mixture>("Mixture")

  // constructor
  .constructor<arma::mat, arma::ivec>("Create new Mixture")

  .field("K", &Mixture::K)
  .field("z", &Mixture::z)

  .method("add_sample", &Mixture::add_sample)
  .method("rm_sample", &Mixture::rm_sample)
  .method("collapsed_gibbs", &Mixture::collapsed_gibbs)
  .method("split_merge", &Mixture::split_merge)
  ;
}
