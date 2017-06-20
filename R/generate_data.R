generate_from_2D_GMM = function(N, mu_list, Sigma_list = NULL){
  D = 2
  K = length(mu_list)
  X = matrix(0, N, D)
  z = sample(1:K, N, replace=TRUE)
  if(is.null(Sigma_list)){
    Sigma_list <- lapply(mu_list, function(x)diag(length(x)))
  }
  for(k in 1:length(mu_list)){
    subset <- z == k
    N_k <- sum(subset)
    if(N_k > 0){
      X[subset, ] <- mvtnorm::rmvnorm(N_k, mu_list[[k]], Sigma_list[[k]])
    }
  }
  list(X = X, z = z)
}
