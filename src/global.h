#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

double logdet_chol_fast(arma::mat L, int D);
double loglik_marginal_NIW_fast(int N, int D, double kappa, double nu, arma::mat S_chol, double kappa0, double nu0, arma::mat S0_chol);

// cholesky rank one update and downdate
arma::mat chol_update_arma(arma::mat LL, arma::vec xx, int D);
arma::mat chol_downdate(arma::mat LL, arma::vec xx, int D);

arma::vec softmax(arma::vec logx);

double logdet_chol_fast(arma::mat L);

// generate from distributions
arma::mat riwishart(unsigned int df, const arma::mat& S);
arma::mat rmvnorm_arma(int n, arma::vec mu, arma::mat sigma);
arma::vec dmvnrm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd);
