#include "global.h"
#include "Component.h"
#include "Mixture.h"

RCPP_EXPOSED_CLASS(Component);

RCPP_MODULE(module_Component){
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
  .field("mu", &Component::mu)
  .field("Sigma", &Component::Sigma)

  .method("is_empty", &Component::is_empty)
  .method("add_sample", &Component::add_sample)
  .method("rm_sample", &Component::rm_sample)
  .method("marginal_loglik", &Component::marginal_loglik)
  .method("posterior_predictive", &Component::posterior_predictive)
  .method("update_IW_pars", &Component::update_IW_pars)
  .method("get_S", &Component::get_S)
  .method("reinitialise", &Component::reinitialise)
  ;
}

RCPP_EXPOSED_CLASS(Mixture);

RCPP_MODULE(module_Mixture){
  using namespace Rcpp;

  // exposing the class in R as "Component"
  class_<Mixture>("Mixture")

  // constructor
  .constructor<arma::mat, arma::ivec, bool, bool>("Create new Mixture")

  .field("K", &Mixture::K)
  .field("N", &Mixture::N)
  .field("z", &Mixture::z)
  .field("X", &Mixture::X)
  .field("alpha", &Mixture::alpha)

  .method("add_sample", &Mixture::add_sample)
  .method("rm_sample", &Mixture::rm_sample)
  .method("collapsed_gibbs", &Mixture::collapsed_gibbs)
  .method("collapsed_gibbs_obs_i", &Mixture::collapsed_gibbs_obs_i)
  .method("split_merge", &Mixture::split_merge)
  .method("get_marginal_loglik", &Mixture::get_marginal_loglik)
  .method("update_X", &Mixture::update_X)
  .method("generate_sample", &Mixture::generate_sample)
  .method("get_z", &Mixture::get_z)
  .method("get_component", &Mixture::get_component)
  .method("update_alpha", &Mixture::update_alpha)
  ;
}

