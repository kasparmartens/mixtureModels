#include "global.h"
#include "Component.h"


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
