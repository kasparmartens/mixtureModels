#ifndef COMPONENT_H
#define COMPONENT_H

#include "global.h"

class Component {
public:
  int N;
  int D;
  double kappa0;
  double kappa;
  double nu0;
  double nu;
  arma::vec m0;
  arma::vec m;
  arma::mat L0;
  arma::mat L;
  arma::mat S0;
  arma::mat S;
  double alpha;

  Component(int D);
  Component(int D, arma::mat X);

  bool is_empty();
  void add_sample(arma::vec x);
  void rm_sample(arma::vec x);

  double marginal_loglik();
  double posterior_predictive(arma::vec x);
};

#endif