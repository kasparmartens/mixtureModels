#include "global.h"

// obtained from: https://github.com/coatless/r-to-armadillo/blob/master/src/distributions.cpp

arma::mat rwishart(unsigned int df, const arma::mat& S){
  // Dimension of returned wishart
  unsigned int m = S.n_rows;

  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);

  // Fill the diagonal
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }

  // Fill the lower matrix with random guesses
  for(unsigned int j = 0; j < m; j++){
    for(unsigned int i = j+1; i < m; i++){
      Z(i,j) = R::rnorm(0,1);
    }
  }

  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);

  // Return random wishart
  return C.t()*C;
}

arma::mat riwishart(unsigned int df, const arma::mat& S){
  return rwishart(df,S.i()).i();
}

// obtained from http://gallery.rcpp.org/articles/simulate-multivariate-normal/

arma::mat rmvnorm_arma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
