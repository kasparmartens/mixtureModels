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
  }
}

void Mixture::update_X(arma::mat XX){
  X = XX;
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
    add_component();
  }
  // add X[i, ] to component k
  components[k]->add_sample(arma::conv_to<arma::vec>::from(X.row(i)));
}

void Mixture::rm_sample(int i){
  int k = z[i];
  // remove X[i, ] from component k
  components[k]->rm_sample(arma::conv_to<arma::vec>::from(X.row(i)));
  // if the component is empty now, remove it
  if(components[k]->is_empty()){
    rm_component(k);
  }
}

void Mixture::add_component(){
  K = K + 1;
  components.push_back(make_unique<Component>(D));
}

void Mixture::add_component(arma::uvec ind){
  z.elem(ind).fill(K);
  components.push_back(make_unique<Component>(D, X.rows(ind)));
  K = K + 1;
}

void Mixture::rm_component(int k){
  components.erase(components.begin()+k);
  z.elem(find(z > k)) -= 1;
  K = K - 1;
}

NumericVector Mixture::get_marginal_loglik(){
  NumericVector res(K);
  for(int k=0; k<K; k++){
    res[k] = components[k]->marginal_loglik();
  }
  return res;
}

int Mixture::collapsed_gibbs_obs_i(int i){
  arma::vec x = arma::conv_to<arma::vec>::from(X.row(i));
  arma::vec logprobs(K + 1);
  // int prev_z_i = z[i];
  z[i] = -1;
  // existing clusters [0, ..., K-1]
  for(int k = 0; k < K; k++){
    arma::uvec temp = find(z == k);
    logprobs[k] = log(temp.size()) + components[k]->posterior_predictive(x);
    // Component comp(D, X.rows(temp));
    // double diffS = accu(abs(comp.get_S() - components[k]->get_S()));
    // printf("diffS %g \n", diffS);
  }
  // new cluster K
  logprobs[K] = log(alpha) + empty_component.posterior_predictive(x);
  arma::vec probs = softmax(logprobs);
  // arma::vec probs_unchanged(probs.begin(), probs.size(), true);

  IntegerVector sequence = seq_len(K+1)-1;
  int k0 = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(sequence, 1, false, probs));
  return k0;
}

void Mixture::check_empty_clusters(){
  for(int k=0; k<K; k++){
    if(components[k]->is_empty()){
      printf("k = %d is empty \n", k);
    }
  }
}

void Mixture::collapsed_gibbs(){
  for(int i=0; i<N; i++){
    rm_sample(i);
    check_empty_clusters();
    int k0 = collapsed_gibbs_obs_i(i);
    check_empty_clusters();
    add_sample(i, k0);
  }
}

void Mixture::split_merge(){
  IntegerVector sequence = seq_len(N)-1;
  IntegerVector indexes = Rcpp::RcppArmadillo::sample(sequence, 2, false);
  int i = indexes[0];
  int j = indexes[1];
  if(z[i] == z[j]){
    propose_split(i, j);
  } else{
    propose_merge(i, j);
  }
}

void Mixture::propose_split(int i, int j){
  Component S_i(D);
  Component S_j(D);
  S_i.add_sample(X.row(i).t());
  S_j.add_sample(X.row(j).t());
  arma::uvec S_ind = find(z == z[i]);
  int n_elements = S_ind.size();
  arma::uvec permutation = Rcpp::RcppArmadillo::sample(S_ind, n_elements, false);
  arma::ivec temp_z;
  temp_z.zeros(N);
  temp_z[i] = 1;
  temp_z[j] = 2;
  double MH_logratio = 0.0;
  for(int k=0; k<n_elements; k++){
    int index = permutation[k];
    if(index == i || index == j){
      // do nothing
    } else{
      arma::vec x = arma::conv_to<arma::vec>::from(X.row(index));
      double p_i = S_i.get_N() * exp(S_i.posterior_predictive(x));
      double p_j = S_j.get_N() * exp(S_j.posterior_predictive(x));
      double prob_i = p_i / (p_i + p_j);
      if(R::runif(0, 1) < prob_i){
        S_i.add_sample(x);
        temp_z[index] = 1;
        MH_logratio += log(prob_i);
      } else{
        S_j.add_sample(x);
        temp_z[index] = 2;
        MH_logratio += log(1-prob_i);
      }
    }
  }
  double logprob_proposed = S_i.marginal_loglik() + S_j.marginal_loglik();
  double logprob_current = components[z[i]]->marginal_loglik();
  MH_logratio = logprob_proposed - logprob_current - MH_logratio;
  MH_logratio += log(alpha) + Rf_lgammafn(S_i.get_N()) + Rf_lgammafn(S_j.get_N()) - Rf_lgammafn(S_i.get_N() + S_j.get_N());
  if(R::runif(0, 1) < exp(MH_logratio)){
    int prev_z_i = z[i];
    add_component(find(temp_z == 1));
    add_component(find(temp_z == 2));
    rm_component(prev_z_i);
  }
}

void Mixture::propose_merge(int i, int j){
  arma::uvec S_ind = find((z == z[i]) + (z == z[j]) == 1);
  Component S_merged(D, X.rows(S_ind));

  Component S_i(D);
  Component S_j(D);
  S_i.add_sample(X.row(i).t());
  S_j.add_sample(X.row(j).t());
  // arma::uvec S_ind = find(z == z[i]);
  int n_elements = S_ind.size();
  arma::uvec permutation = Rcpp::RcppArmadillo::sample(S_ind, n_elements, false);
  // arma::ivec temp_z(N).zeros();
  // temp_z[i] = 1;
  // temp_z[j] = 2;
  double MH_logratio = 0.0;
  for(int k=0; k<n_elements; k++){
    int index = permutation[k];
    if(index == i || index == j){
      // do nothing
    } else{
      arma::vec x = arma::conv_to<arma::vec>::from(X.row(index));
      double p_i = S_i.get_N() * exp(S_i.posterior_predictive(x));
      double p_j = S_j.get_N() * exp(S_j.posterior_predictive(x));
      double prob_i = p_i / (p_i + p_j);
      if(z[index] == z[i]){
        S_i.add_sample(x);
        // temp_z[index] = 1;
        MH_logratio += log(prob_i);
      } else if(z[index] == z[j]){
        S_j.add_sample(x);
        // temp_z[index] = 2;
        MH_logratio += log(1-prob_i);
      } else{
        printf("something went wrong\n");
      }
    }
  }
  double logprob_proposed = S_merged.marginal_loglik();
  double logprob_current = S_i.marginal_loglik() + S_j.marginal_loglik();
  MH_logratio = logprob_proposed - logprob_current + MH_logratio;
  MH_logratio += -log(alpha) - Rf_lgammafn(S_i.get_N()) - Rf_lgammafn(S_j.get_N()) + Rf_lgammafn(S_merged.get_N());
  if(R::runif(0, 1) < exp(MH_logratio)){
    int prev1 = std::min(z[i], z[j]);
    int prev2 = std::max(z[i], z[j]);
    add_component(S_ind);
    rm_component(prev2);
    rm_component(prev1);
  }
}

Rcpp::List Mixture::generate_sample(int n){
  arma::mat out(n, D);
  arma::vec cluster_probs(K);
  for(int k=0; k<K; k++){
    cluster_probs[k] = components[k]->get_N();
  }
  // IntegerVector seq = seq_len(K)-1;
  arma::uvec seq = arma::linspace<arma::uvec>(0, K-1, K);
  arma::uvec cluster_alloc = Rcpp::RcppArmadillo::sample(seq, n, true, cluster_probs);
  for(int k=0; k<K; k++){
    arma::uvec which_ind = find(cluster_alloc == k);
    int n_k = which_ind.size();
    if(n_k > 0){
      components[k]->update_IW_pars();
      out.rows(which_ind) = rmvnorm_arma(n_k, components[k]->get_mu(), components[k]->get_Sigma());
    }
  }
  // compute p(x | mu, Sigma) for each point
  arma::vec loglik;
  loglik.zeros(n);
  for(int k=0; k<K; k++){
    arma::uvec which_ind = find(cluster_alloc == k);
    int n_k = which_ind.size();
    if(n_k > 0){
      double pi_k = (double) n_k / N;
      loglik += log(pi_k) + dmvnrm_arma(out, components[k]->get_mu().t(), components[k]->get_Sigma(), true);
    }
  }
  return Rcpp::List::create(Named("X") = out,
                            Named("z") = 1 + IntegerVector(cluster_alloc.begin(), cluster_alloc.end()),
                            Named("loglik") = NumericVector(loglik.begin(), loglik.end()));
}

NumericVector Mixture::get_z(){
  return NumericVector(z.begin(), z.end()) + 1;
}

Component Mixture::get_component(int k){
  return *components[k];
}
