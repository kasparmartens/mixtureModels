logdet_chol = function(L) 2*sum(log(diag(L)))

loglik_marginal_NIW <- function(N, D, kappa, nu, S_chol, kappa0, nu0, S0_chol){
  -0.5*N*D*log(pi) - 0.5*D*(log(kappa) - log(kappa0)) - 0.5*(nu*logdet_chol(S_chol) - nu0*logdet_chol(S0_chol)) + sum(lgamma(0.5*(nu - 1:D + 1)) - lgamma(0.5*(nu0 - 1:D + 1)))
}

add_sample_and_get_cholS_helper <- function(x, m, kappa, nu, L){
  m2 <- m + (x - m) / (kappa + 1L)
  chol_update(L, sqrt((kappa + 1L)/kappa) * (x - m2))
}

softmax <- function(...){
  logx <- unlist(list(...))
  normx <- logx - max(logx)
  exp(normx) / sum(exp(normx))
}
