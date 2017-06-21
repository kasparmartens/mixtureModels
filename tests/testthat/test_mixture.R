context("Testing DP mixture: Gibbs sampler and split-merge moves")

source("relabelling_helper.R")

library(dplyr)

alpha <- 0.1
D <- 2
N <- 5
# X <- matrix(rnorm(N*D), N, D)
X <- generate_from_2D_GMM(N, list(c(-1, -1), c(0, 2), c(3, 0)))$X

partition <- setparts(N)
log_posterior <- rep(NA, ncol(partition))
for(j in 1:ncol(partition)){
  z <- partition[, j]
  K <- max(z)
  N_k_all <- table(z)
  components <- list()
  loglik <- rep(NA, K)
  for(k in 1:K){
    N_k <- sum(z == k)
    components[[k]] <- Component$new(D, X = X[z == k, , drop=FALSE])
    loglik[k] <- components[[k]]$marginal_loglik()
  }
  log_prior <- K * log(alpha) + sum(lgamma(N_k_all)) - sum(log(alpha+1:N-1))
  log_posterior[j] <- log_prior + sum(loglik)
}
probs <- softmax(log_posterior)
pattern <- apply(as.matrix(t(partition)), 1, paste_pattern)
df_correct <- data.frame(pattern = pattern, prob = probs) %>%
  mutate(pattern = as.character(pattern)) %>%
  arrange(desc(prob))

test_that("collapsed gibbs sampler is correct", {
  z <- 1:N
  mixture <- Mixture$new(X, z)
  n_iter <- 10000
  clustering <- matrix(0L, n_iter, N)
  for(i in 1:n_iter){
    for(ii in 1:3){
      mixture$collapsed_gibbs()
    }
    clustering[i, ] <- as.numeric(mixture$z) + 1
  }
  df_empirical <- helper_summarise_partitions(clustering, ref_patterns = df_correct$pattern)
  df_both <- merge(df_correct, df_empirical)
  expect_equal(df_both$prob, df_both$freq, tolerance = 1e-2)
})

test_that("merge split sampler is correct", {
  z <- 1:N
  mixture <- Mixture$new(X, z)
  n_iter <- 10000
  clustering <- matrix(0L, n_iter, N)
  for(i in 1:n_iter){
    for(ii in 1:3){
      mixture$split_merge()
    }
    clustering[i, ] <- as.numeric(mixture$z) + 1
  }
  df_empirical <- helper_summarise_partitions(clustering, ref_patterns = df_correct$pattern)
  df_both <- merge(df_correct, df_empirical)
  expect_equal(df_both$prob, df_both$freq, tolerance = 1e-2)
})
