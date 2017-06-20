Component <- setRefClass("Component",
                         fields = list(
                           N = "integer",
                           D = "integer",
                           kappa = "numeric",
                           nu = "numeric",
                           m = "numeric",
                           L = "matrix",
                           S = "matrix",
                           kappa0 = "numeric",
                           nu0 = "numeric",
                           m0 = "numeric",
                           L0 = "matrix",
                           Sigma = "matrix",
                           mu = "numeric"),

                         methods = list(
                           initialize = function(init_pars, X = NULL){
                             D <<- as.integer(init_pars$D)
                             kappa0 <<- init_pars$kappa0
                             nu0 <<- init_pars$nu0
                             m0 <<- init_pars$m0
                             S0 <- init_pars$S0
                             L0 <<- chol(S0 + kappa0 * outer(m0, m0))
                             N <<- 0L
                             kappa <<- kappa0
                             nu <<- nu0
                             m <<- m0
                             L <<- L0
                             S <<- diag(D)
                             if(!is.null(X)){
                               N_k <- nrow(X)
                               kappa_k <- kappa + N_k
                               m_k <- (kappa*m + colSums(X)) / kappa_k
                               S_k <- S + kappa * outer(m, m) + t(X) %*% X
                               L <<- chol(S_k - kappa_k * outer(m_k, m_k))
                               if(any(is.nan(L))) stop("init problem")
                               N <<- N + N_k
                               nu <<- nu + N_k
                               kappa <<- kappa_k
                               m <<- m_k
                               S <<- S_k
                             }
                           },
                           is_empty = function(){
                             ifelse(N == 0, TRUE, FALSE)
                           },
                           add_sample = function(x){
                             kappa <<- kappa + 1L
                             m <<- ((kappa - 1)*m + x) / kappa
                             nu <<- nu + 1L
                             L <<- chol_update(L, sqrt(kappa/(kappa - 1)) * (x - m))
                             if(any(is.nan(L))) stop("Add sample")
                             N <<- N + 1L
                             S <<- S + outer(x, x)
                             testthat::expect_equal(chol(S - kappa * outer(m, m)), L)
                           },
                           rm_sample = function(x){
                             L <<- chol_downdate(L, sqrt(kappa/(kappa - 1)) * (x - m))
                             if(any(is.nan(L))){

                               stop("Rm sample")
                             }
                             kappa <<- kappa - 1L
                             m <<- ((kappa + 1)*m - x) / kappa
                             nu <<- nu - 1L
                             N <<- N - 1L
                             S <<- S - outer(x, x)
                             testthat::expect_equal(chol(S - kappa * outer(m, m)), L)
                             # if(mean(abs(chol(S - kappa * outer(m, m)) - L)) > 1e-5) stop("We have a problem")
                           },
                           get_cholS = function(){
                             # subtracts the cluster means from the covariance matrix
                             # and returns Cholesky decomposition
                             # LL <- chol(S - kappa * outer(m, m))
                             return(L)
                           },
                           get_S = function(){
                             S - kappa * outer(m, m)
                           },
                           add_sample_and_get_cholS = function(x){
                             # call the generic function
                             return(add_sample_and_get_cholS_helper(x, m, kappa, nu, L))
                           },
                           marginal_loglik = function(){
                             loglik_marginal_NIW(N, D, kappa, nu, L, kappa0, nu0, L0)
                           },
                           posterior_predictive = function(x){
                             # prediction when data point x were in the current cluster
                             L_updated <- add_sample_and_get_cholS(x)
                             # calculate loglikelihood
                             loglik <- loglik_marginal_NIW(1, D, kappa + 1L, nu + 1L, L_updated, kappa, nu, L)
                             # loglik <- -0.5*D*log(pi) - 0.5*D*log((kappa+1)/kappa) - 0.5*(nu+1)*logdet_chol(L_updated) + lgamma(0.5*(nu + 1)) + 0.5*nu*logdet_chol(L_current) - lgamma(0.5*(nu + 1 - D))
                             loglik
                           },
                           update_IW_pars = function(){
                             Sigma <<- solve(rWishart(1, nu, solve(get_S()))[, , 1])
                             mu <<- as.numeric(mvtnorm::rmvnorm(1, mean = m, sigma = 1/kappa * Sigma))
                           }
                         ))
