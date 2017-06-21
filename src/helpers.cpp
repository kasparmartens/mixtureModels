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
