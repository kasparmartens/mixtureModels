#include "Mixture.h"

#include <RcppArmadilloExtensions/sample.h>

Mixture::Mixture(arma::mat XX, arma::ivec zz) : empty_component(XX.n_cols){
  // data matrix X
  N = XX.n_rows;
  D = XX.n_cols;
  X = XX;
  // cluster allocations [0, ..., K-1]
  z = zz - 1;
  K = z.max()+1;
  alpha = 0.1;

  for(int k = 0; k < K; k++){
    arma::mat X_k = X.rows(find(z == k));

    components.push_back(make_unique<Component>(D, X_k));
    printf("k = %d, N = %d\n", k, components[k]->get_N());
  }
}

void Mixture::update_X(arma::mat XX){
  // update the X coordinates, leave cluster assignments unchanged
  for(int k = 0; k < K; k++){
    arma::mat X_k = XX.rows(find(z == k));
    components[k]->reinitialise(X_k);
  }
}

void Mixture::add_sample(int i, int k){
  // update cluster label
  z[i] = k;
  // is it necessary to create a new component?
  if(k > K-1){
    K = K + 1;
    components.push_back(make_unique<Component>(D));
  }
  // add X[i, ] to component k
  arma::vec x = arma::conv_to<arma::vec>::from(X.row(i));
  components[k]->add_sample(x);
}

void Mixture::rm_sample(int i){
  int k = z[i];
  // remove X[i, ] from component k
  arma::vec x = arma::conv_to<arma::vec>::from(X.row(i));
  components[k]->rm_sample(x);
  // if the component is empty now, remove it
  if(components[k]->is_empty()){
    components.erase(components.begin()+k);
    z.elem(find(z > k)) -= 1;
    K = K - 1;
  }
}

int Mixture::collapsed_gibbs_obs_i(int i){
  arma::vec x = arma::conv_to<arma::vec>::from(X.row(i));
  arma::vec logprobs(K + 1);
  int count;
  // int prev_z_i = z[i];
  z[i] = -1;
  // existing clusters [0, ..., K-1]
  for(int k = 0; k < K; k++){
    // int count = sum(z == k);
    arma::uvec temp = (z == k);
    logprobs[k] = log(sum(temp)) + components[k]->posterior_predictive(x);
  }
  // new cluster K
  logprobs[K] = log(alpha) + empty_component.posterior_predictive(x);
  arma::vec probs = softmax(logprobs);
  // arma::vec probs_unchanged(probs.begin(), probs.size(), true);

  IntegerVector sequence = seq_len(K+1)-1;
  int k0 = Rcpp::as<int>(Rcpp::RcppArmadillo::sample_main(sequence, 1, false, probs));
  return k0;
}

void Mixture::collapsed_gibbs(){
  for(int i=0; i<N; i++){
    rm_sample(i);
    int k0 = collapsed_gibbs_obs_i(i);
    add_sample(i, k0);
  }
}
