#include "global.h"

double logdet_chol_fast(arma::mat L, int D){
  double sum = 0.0;
  for(int i=0; i<D; i++){
    sum += log(L(i, i));
  }
  return 2*sum;
}

double loglik_marginal_NIW_fast(int N, int D, double kappa, double nu, arma::mat S_chol, double kappa0, double nu0, arma::mat S0_chol){
  const double pi = 3.141592653589793238462643383280;
  double res = -0.5*N*D*log(pi) - 0.5*D*(log(kappa) - log(kappa0)) - 0.5*(nu*logdet_chol_fast(S_chol, D) - nu0*logdet_chol_fast(S0_chol, D));
  for(int d=0; d<D; d++){
    res += Rf_lgammafn(0.5*(nu - d + 1)) - Rf_lgammafn(0.5*(nu0 - d + 1));
  }
  return res;
}


arma::mat chol_update_arma(arma::mat LL, arma::vec xx, int D) {
  arma::mat L(LL.memptr(), D, D, true);
  arma::vec x(xx.memptr(), D, true);
  int n = x.size();
  for(int k=0; k<n; k++){
    double r = sqrt(L(k, k)*L(k, k) + x[k]*x[k]);
    double c = r / L(k, k);
    double s = x[k] / L(k, k);
    L(k, k) = r;
    for(int j=k+1; j<n; j++){
      L(k, j) = (L(k, j) + s*x[j])/c;
      x[j] = c*x[j] - s*L(k, j);
    }
  }
  return L;
}

arma::mat chol_downdate(arma::mat LL, arma::vec xx, int D) {
  arma::mat L(LL.memptr(), D, D, true);
  arma::vec x(xx.memptr(), D, true);
  int n = x.size();
  for(int k=0; k<n; k++){
    double r = sqrt(L(k, k)*L(k, k) - x[k]*x[k]);
    double c = r / L(k, k);
    double s = x[k] / L(k, k);
    L(k, k) = r;
    for(int j=k+1; j<n; j++){
      L(k, j) = (L(k, j) - s*x[j])/c;
      x[j] = c*x[j] - s*L(k, j);
    }
  }
  return L;
}
