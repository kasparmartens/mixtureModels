#ifndef MIXTURE_H
#define MIXTURE_H

#include "global.h"
#include "Component.h"

#include <vector>
typedef std::vector<std::unique_ptr<Component>> container;

class Mixture {
public:
  int N;
  int D;
  int K;
  bool is_DPM;
  bool is_MFM;
  arma::ivec z;
  arma::mat X;
  // only for MFM
  arma::vec log_V_n;

  double kappa0;
  double nu0;
  arma::vec m0;
  arma::mat L0;
  arma::mat S0;
  double alpha;

  container components;

  Component empty_component;

  Mixture(arma::mat X, arma::ivec z, bool is_DPM, bool is_MFM);

  void add_sample(int i, int k);
  void rm_sample(int i);
  void update_X(arma::mat X);

  void add_component();
  void add_component(arma::uvec X);
  void rm_component(int k);

  int collapsed_gibbs_obs_i(int i);
  void collapsed_gibbs();

  void split_merge();
  void propose_split(int i, int j);
  void propose_merge(int i, int j);
  
  void update_alpha(int n_steps);

  double get_loglik(arma::rowvec x);
  NumericVector get_marginal_loglik();
  Rcpp::List generate_sample(int n);
  Component get_component(int k);

  NumericVector get_z();

  void check_empty_clusters();
};


#endif
