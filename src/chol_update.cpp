#include "global.h"

// Cholesky rank one update and downdate to be used in cpp

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
    if(isnan(r)){
      printf("Error: chol downdate problem!\n");
    }
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


// Cholesky rank one update and downdate to be used in R

//' @export
// [[Rcpp::export]]
NumericMatrix chol_update(NumericMatrix LL, NumericVector xx) {
  NumericMatrix L = Rcpp::clone(LL);
  NumericVector x = Rcpp::clone(xx);
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

//' @export
// [[Rcpp::export]]
NumericMatrix chol_downdate(NumericMatrix LL, NumericVector xx) {
  NumericMatrix L = Rcpp::clone(LL);
  NumericVector x = Rcpp::clone(xx);
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
