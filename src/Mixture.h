#ifndef MIXTURE_H
#define MIXTURE_H

#include "global.h"
#include "Component.h"

#include <memory>
#include <vector>
typedef std::vector<std::unique_ptr<Component>> container;

class Mixture {
public:
  int N;
  int D;
  int K;
  arma::ivec z;
  arma::mat X;

  double kappa0;
  double nu0;
  arma::vec m0;
  arma::mat L0;
  arma::mat S0;
  double alpha;

  container components;

  Component empty_component;

  Mixture(arma::mat X, arma::ivec z);

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
};


#endif
