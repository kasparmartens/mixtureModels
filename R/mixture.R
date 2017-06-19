Rcpp::loadModule("RcppComponent", TRUE)

Mixture <- setRefClass("Mixture",
                       fields = list(
                         components = "list",
                         K = "integer",
                         D = "integer",
                         X = "matrix",
                         init_pars = "list",
                         N_total = "integer",
                         z = "integer",
                         kappa0 = "numeric",
                         nu0 = "numeric",
                         m0 = "numeric",
                         L0 = "matrix",
                         S0 = "matrix",
                         alpha = "numeric"),

                       methods = list(
                         initialize = function(K, D, X, z, kappa0 = 1.0, nu0 = D+2, m0 = rep(0, D), S0 = diag(rep(1, D)), alpha = 0.1){
                           K <<- as.integer(K)
                           D <<- as.integer(D)
                           N_total <<- nrow(X)
                           z <<- as.integer(z)
                           X <<- X
                           alpha <<- alpha

                           kappa0 <<- kappa0
                           nu0 <<- nu0
                           m0 <<- m0
                           S0 <<- S0
                           L0 <<- chol(S0 + kappa0 * outer(m0, m0))

                           init_pars <<- list(D = D, kappa0 = kappa0, nu0 = nu0, m0 = m0, S0 = S0, L0 = L0)
                           for(k in 1:K){
                             components[[k]] <<- Component$new(init_pars, X = X[z == k, , drop=FALSE])
                           }
                         },
                         update_X = function(X){
                           # given a new X value, but conditioning on the current clustering z,
                           # we need to update m and L for each component
                           # (the rest will remain constant)
                           X <<- X
                           for(k in 1:K){
                             X_k <- X[z == k, , drop=FALSE]
                             N_k <- sum(z == k)
                             kappa_k <- (kappa0 + N_k)
                             m_k <- (kappa0 * m0 + colSums(X_k)) / kappa_k
                             S_k <- S0 + kappa0 * outer(m0, m0) + t(X_k) %*% X_k
                             components[[k]]$m <<- m_k
                             components[[k]]$L <<- chol(S_k - kappa_k * outer(m_k, m_k))
                             components[[k]]$S <<- S_k
                           }
                         },
                         add_sample = function(i, k){
                           x <- X[i, ]
                           z[i] <<- k
                           if(k > K){
                             new_component()
                           }
                           components[[k]]$add_sample(x)
                         },
                         rm_sample = function(i){
                           x <- X[i, ]
                           k <- z[i]
                           components[[k]]$rm_sample(x)
                           if(components[[k]]$is_empty()){
                             rm_component(k)
                           }
                         },
                         new_component = function(which_ind = NULL){
                           K <<- K + 1L
                           if(is.null(which_ind)){
                             components[[K]] <<- Component$new(init_pars)
                           }else{
                             components[[K]] <<- Component$new(init_pars, X = X[which_ind, , drop=FALSE])
                             z[which_ind] <<- K
                           }
                         },
                         rm_component = function(k){
                           components[[k]] <<- NULL
                           K <<- K - 1L
                           # update the cluster numbers
                           z <<- ifelse(z > k, z-1L, z)
                         },
                         collapsed_gibbs_for_obs_i = function(i){
                           x <- X[i, ]
                           logprobs <- rep(0, K+1)
                           # for existing clusters
                           for(k in 1:K){
                             logprior <- log(sum(z[-i] == k))
                             loglik <- components[[k]]$posterior_predictive(x)
                             logprobs[k] <- logprior + loglik
                             # cat("k", k, "logprior", logprior, "loglik", loglik, "\n")
                           }
                           # cat("logprobs" ,logprobs, "\n")
                           # for a new cluster
                           Ltemp <- add_sample_and_get_cholS_helper(x, m0, kappa0, nu0, L0)
                           # temp <- helper_likelihood(1, D, kappa0+1, nu0+1, Ltemp) - helper_likelihood(0, D, kappa0, nu0, L0)
                           temp <- -0.5*D*log(pi) - 0.5*D*log((kappa0+1)/kappa0) - 0.5*(nu0+1)*logdet_chol(Ltemp) + lgamma(0.5*(nu0 + 1)) +
                             0.5*nu0*logdet_chol(L0) - lgamma(0.5*(nu0 + 1 - D))
                           logprobs[K+1] <- log(alpha) + temp
                           k0 = base::sample(1:(K+1), 1, prob = softmax(logprobs))
                           k0
                         },
                         collapsed_gibbs = function(){
                           for(i in 1:N_total){
                             rm_sample(i)
                             k0 <- collapsed_gibbs_for_obs_i(i)
                             add_sample(i, k0)
                           }
                         },
                         merge_split = function(){
                           i <- sample(1:N_total, 1)
                           j <- sample(setdiff(1:N_total, i), 1)
                           if(z[i] == z[j]){
                             propose_split(i, j, components[[z[i]]])
                           } else{
                             propose_merge(i, j)
                           }
                         },
                         propose_split = function(i, j, S_current){
                           S_i <- Component$new(init_pars)
                           S_j <- Component$new(init_pars)
                           S_i$add_sample(X[i, ])
                           S_j$add_sample(X[j, ])
                           S_ind <- which(z %in% c(z[i], z[j]))
                           # temp_z <- rep(0L, length(S_ind))
                           # temp_z[S_ind == i] <- 1L
                           # temp_z[S_ind == j] <- 2L
                           temp_z <- rep(0L, N_total)
                           temp_z[c(i, j)] <- c(1L, 2L)
                           MH_logratio <- 0
                           if(length(S_ind) > 2){
                             ind_perm <- sample((S_ind))
                             # temp cluster allocations within S_ind
                             for(kk in setdiff(ind_perm, c(i, j))){
                               # kk <- S_ind[k]
                               x <- X[kk, ]
                               # choose whether add observation k to S_i or S_j
                               p_i <- S_i$N * exp(S_i$posterior_predictive(x))
                               p_j <- S_j$N * exp(S_j$posterior_predictive(x))
                               prob_i <- p_i / (p_i + p_j)
                               if(runif(1) < prob_i){
                                 S_i$add_sample(x)
                                 temp_z[kk] <- 1L
                                 MH_logratio <- MH_logratio + log(prob_i)
                               } else{
                                 S_j$add_sample(x)
                                 temp_z[kk] <- 2L
                                 MH_logratio <- MH_logratio + log(1-prob_i)
                               }
                             }
                           }
                           logprob_proposed <- S_j$marginal_loglik() + S_i$marginal_loglik()
                           logprob_current <- S_current$marginal_loglik()
                           MH_logratio <- - MH_logratio + logprob_proposed - logprob_current + log(alpha) + lgamma(S_i$N) + lgamma(S_j$N) - lgamma(S_i$N + S_j$N)
                           # accept or reject the constructed proposal
                           if(runif(1) < exp(MH_logratio)){
                             # cat("accepted split cluster", z[i], "with prob", exp(MH_logratio), "\n")
                             rm_component(z[i])
                             # cat("splitting elements", S_ind, "into", S_ind[temp_z == 1L], "and", S_ind[temp_z == 2L], "\n")
                             # new_component(which_ind = which(temp_z == 1L))
                             components[[K+1]] <<- S_i
                             components[[K+2]] <<- S_j
                             z[temp_z == 1L] <<- K+1L
                             z[temp_z == 2L] <<- K+2L
                             K <<- K + 2L
                             # new_component(which_ind = which(temp_z == 2L))
                             # cat("new splits:", which(temp_z == 1L), "and", which(temp_z == 2L), "\n")
                           }
                         },
                         propose_merge = function(i, j){
                           S_ind <- which(z %in% c(z[i], z[j]))
                           S_merged <- Component$new(init_pars, X = X[S_ind, , drop=FALSE])
                           S_i <- Component$new(init_pars)
                           S_j <- Component$new(init_pars)
                           S_i$add_sample(X[i, ])
                           S_j$add_sample(X[j, ])
                           MH_logratio <- 0
                           if(length(S_ind) > 2){
                             # imaginary clusters S_i and S_j
                             ind_perm <- sample(setdiff(S_ind, c(i, j)))
                             for(k in 1:length(ind_perm)){
                               kk <- ind_perm[k]
                               x <- X[kk, ]
                               # choose whether add observation k to S_i or S_j
                               p_i <- S_i$N * exp(S_i$posterior_predictive(x))
                               p_j <- S_j$N * exp(S_j$posterior_predictive(x))
                               prob_i <- p_i / (p_i + p_j)
                               if(z[kk] == z[i]){
                                 S_i$add_sample(x)
                                 MH_logratio <- MH_logratio + log(prob_i)
                               } else{
                                 S_j$add_sample(x)
                                 MH_logratio <- MH_logratio + log(1-prob_i)
                               }
                             }
                           }
                           logprob_current <- S_j$marginal_loglik() + S_i$marginal_loglik()
                           logprob_proposed <- S_merged$marginal_loglik()
                           MH_logratio <- MH_logratio + logprob_proposed - logprob_current - log(alpha) - lgamma(S_i$N) - lgamma(S_j$N) + lgamma(S_i$N + S_j$N)
                           if(runif(1) < exp(MH_logratio)){
                             # cat("accepted merge with prob", exp(MH_logratio), "\n")
                             rm_component(z[i])
                             rm_component(z[j])
                             new_component(which_ind = S_ind)
                           }
                           rm(S_i, S_j, S_merged)
                         },
                         generate_sample = function(n){
                           cluster_probs <- sapply(components, function(x)x$N)
                           cluster_allocations <- base::sample(1:K, n, replace=T, prob = cluster_probs)
                           out <- matrix(0, n, D)
                           # for all existing clusters
                           mu_list <- list()
                           Sigma_list <- list()
                           for(k in 1:K){
                             components[[k]]$update_IW_pars()
                             Sigma_list[[k]] <- components[[k]]$Sigma
                             mu_list[[k]] <- components[[k]]$mu
                             subset <- (cluster_allocations == k)
                             if(sum(subset) > 0){
                               out[subset, ] <- mvtnorm::rmvnorm(sum(subset), mean = mu_list[[k]], sigma = Sigma_list[[k]])
                             }
                           }
                           # # for the new cluster
                           # subset <- (cluster_allocations == k+1)
                           # if(sum(subset) > 0){
                           #   Sigma <- solve(rWishart(1, nu0, solve(S0))[, , 1])
                           #   mu <- as.numeric(mvtnorm::rmvnorm(1, mean = m0, sigma = 1/kappa0 * Sigma))
                           #   cat("Started a new cluster, with mean", mu, "and covariance ", Sigma, "\n")
                           #   out[subset, ] <- mvtnorm::rmvnorm(sum(subset), mean = mu, sigma = Sigma)
                           # }
                           list(X = out, z = cluster_allocations, mu = mu_list, Sigma = Sigma_list)
                         },
                         check_L = function(){
                           for(k in 1:K){
                             X_k <- X[z == k, , drop=FALSE]
                             N_k <- nrow(X_k)
                             kappa_k <- kappa0 + N_k
                             m_k <- (kappa0*m0 + colSums(X_k)) / kappa_k
                             S_k <- S0 + t(X_k) %*% X_k
                             L_k <- chol(S_k - kappa_k * outer(m_k, m_k))
                             testthat::expect_equal(L_k, components[[k]]$L)
                           }
                         },
                         plot = function(){
                           # r1 <- extendrange(X[, 1])
                           # r2 <- extendrange(X[, 2])
                           # xval <- seq(r1[1], r1[2], 0.02)
                           # yval <- seq(r2[1], r2[2], 0.02)
                           # dftemp <- melt(outer(xval , yval))
                           # x <- xval[df.temp$Var1]
                           # y <- yval[df.temp$Var2]
                           # mat <- cbind(x, y)
                           # df <- data.frame(c())
                           mu_list <- list()
                           Sigma_list <- list()
                           for(k in 1:K){
                             components[[k]]$update_IW_pars()
                             Sigma_list[[k]] <- components[[k]]$Sigma
                             mu_list[[k]] <- components[[k]]$mu
                           }
                           plot_mixture(mu_list, Sigma_list, X, z)
                         }
                       ))

add_sample_and_get_cholS_helper = function(x, m, kappa, nu, L){
  m2 <- m + (x - m) / (kappa + 1L)
  chol_update(L, sqrt((kappa + 1L)/kappa) * (x - m2))
}
