#include "global.h"

double logdet_chol_fast(arma::mat L){
  return 2*sum(log(L.diag()));
}

double loglik_marginal_NIW_fast(int N, int D, double kappa, double nu, arma::mat S_chol, double kappa0, double nu0, arma::mat S0_chol){
  const double pi = 3.141592653589793238462643383280;
  double res = -0.5*N*D*log(pi) - 0.5*D*(log(kappa) - log(kappa0)) - 0.5*(nu*logdet_chol_fast(S_chol) - nu0*logdet_chol_fast(S0_chol));
  for(int d=1; d<=D; d++){
    res += Rf_lgammafn(0.5*(nu - d + 1)) - Rf_lgammafn(0.5*(nu0 - d + 1));
  }
  return res;
}

arma::vec softmax(arma::vec logx){
  logx -= logx.max();
  arma::vec out = exp(logx) / sum(exp(logx));
  return out;
}

double dDP(double alpha, arma::ivec z, int K){
  double logprob = K * log(alpha);
  for(int k=0; k<K; k++){
    arma::uvec which_ind = find(z == k);
    int N_k = which_ind.size();
    logprob += Rf_lgammafn(N_k);
  }
  int N = z.size();
  logprob -= (Rf_lgammafn(alpha + N) - Rf_lgammafn(alpha));
  return logprob;
}

double RWMH_log_scale(double alpha, arma::ivec z, int K, double sd, double a0, double b0){
  double logprob_current = R::dgamma(alpha, a0, 1.0/b0, true) + log(alpha) + dDP(alpha, z, K);
  double alpha_proposal = exp(log(alpha) + R::rnorm(0, sd));
  double logprob_proposal = R::dgamma(alpha_proposal, a0, 1.0/b0, true) + log(alpha_proposal) + dDP(alpha_proposal, z, K);
  if(R::runif(0, 1) < exp(logprob_proposal - logprob_current)){
    alpha = alpha_proposal;
  }
  return alpha;
}

double logsumexp(double a, double b){
  double m = max(a, b);
  if(arma::is_finite(m)){
    return(log(exp(a-m) + exp(b-m)) + m);
  } else{
    return(m);
  }
}

// [[Rcpp::export]]
arma::vec calculate_log_V_n(double alpha, int N, int how_many){
  how_many = min(how_many, N);
  arma::vec log_V_n(how_many);
  double a, b, c;
  for(int t=1; t<=how_many; t++){
    a = 0;
    c = -1 * arma::datum::inf;
    int k = t;
    while(abs(a-c) > 1e-12){
      a = c;
      b = Rf_lgammafn(t+1) - Rf_lgammafn(k-t+1) - Rf_lgammafn(k*alpha+N) + Rf_lgammafn(k*alpha) - log(exp(1)-1) - Rf_lgammafn(k+1);
      c = logsumexp(a, b);
      k += 1;
    }
    log_V_n[t-1] = c;
  }
  return log_V_n;
}
