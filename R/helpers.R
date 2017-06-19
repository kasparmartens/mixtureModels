loglik_marginal_NIW <- function(N, D, kappa, nu, S_chol, kappa0, nu0, S0_chol){
  -0.5*N*D*log(pi) - 0.5*D*(log(kappa) - log(kappa0)) - 0.5*(nu*logdet_chol(S_chol) - nu0*logdet_chol(S0_chol)) + sum(lgamma(0.5*(nu - 1:D + 1)) - lgamma(0.5*(nu0 - 1:D + 1)))
}
