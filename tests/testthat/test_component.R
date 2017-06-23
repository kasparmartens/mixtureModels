context("Testing single component")

test_that("sequential updates work correctly", {
  N <- 5
  D <- 2
  X <- matrix(rnorm(N*D), N, D)
  comp1 <- Component$new(D = D)
  # sequential update
  for(i in 1:N){
    comp1$add_sample(X[i, ])
  }
  # full batch update
  comp2 <- Component$new(D = D, X = X)

  expect_equal(comp1$m, comp2$m)
  expect_equal(comp1$L, comp2$L)
  expect_equal(comp1$get_S(), comp2$get_S())
})

test_that("cholesky downdate works as well", {
  # generate some data
  N <- 5
  D <- 2
  X <- matrix(rnorm(N*D), N, D)
  # initialise components
  comp1 <- Component$new(D = D, X = X)
  comp2 <- Component$new(D = D, X = X)
  # generate one more data point
  x <- runif(D)
  # add and remove this data point
  comp1$add_sample(x)
  comp1$rm_sample(x)

  expect_equal(comp1$m, comp2$m)
  expect_equal(comp1$L, comp2$L)
  expect_equal(comp1$L, chol(comp1$get_S()))
})
