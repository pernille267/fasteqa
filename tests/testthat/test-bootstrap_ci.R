library(fasteqa)
library(testthat)

# Required functions
double_check_resample_zeta <- function(data, x){
  bdata <- resample_samples(data)
  return(estimate_zeta_ols(bdata)$zeta)
}

double_check_resample_lambda <- function(data, x){
  bdata <- resample_samples(data)
  return(global_precision_estimates(bdata)$lambda)
}

# Test Confidence Intervals of Estimated zeta
test_that(desc = "Test bootstrap confidence intervals of estimated zeta", code = {
  orig_zeta <- estimate_zeta_ols(test_data)$zeta
  resampled_zetas <- replicate(n = 1e4, expr = resample_zeta_estimates(test_data)$zeta)
  loo_zetas <- sapply(1:20, FUN = function(x) loo_zeta_estimates(test_data, loo_id = x)$zeta, simplify = TRUE)
  
  lower_bcis <- c(1.514, 1.624, 0.612, 1.341)
  upper_bcsi <- c(4.355, 4.328, 3.316, 3.749)
  
  expected_1 <- c(1.514, 4.355)
  expected_2 <- c(1.624, 4.328)
  expected_3 <- c(0.612, 3.316)
  expected_4 <- c(1.341, 3.749)
  
  actual_1 <- round(bootstrap_ci(resampled_zetas, loo_zetas, orig_zeta, 1), 3)
  actual_2 <- round(bootstrap_ci(resampled_zetas, loo_zetas, orig_zeta, 2), 3)
  actual_3 <- round(bootstrap_ci(resampled_zetas, loo_zetas, orig_zeta, 3), 3)
  actual_4 <- round(bootstrap_ci(resampled_zetas, loo_zetas, orig_zeta, 4), 3)
  
  expect_equal(object = actual_1, expected = expected_1, tolerance = 0.05)
  expect_equal(object = actual_2, expected = expected_2, tolerance = 0.05)
  expect_equal(object = actual_3, expected = expected_3, tolerance = 0.05)
  expect_equal(object = actual_4, expected = expected_4, tolerance = 0.05)
  
})

# Test Confidence Intervals for Estimated lambda
test_that(desc = "Test bootstrap confidence intervals of estimated lambda", code = {
  orig_lambda <- global_precision_estimates(test_data)$lambda
  resampled_lambdas <- replicate(n = 1e4, expr = resample_global_precision_estimates(test_data)$lambda)
  loo_lambdas <- sapply(1:20, FUN = function(x) loo_global_precision_estimates(test_data, loo_id = x)$lambda, simplify = TRUE)
  
  expected_1 <- c(1.520, 3.387)
  expected_2 <- c(1.422, 3.307)
  expected_3 <- c(1.636, 3.521)
  expected_4 <- c(1.660, 3.561)
  
  actual_1 <- round(bootstrap_ci(resampled_lambdas, loo_lambdas, orig_lambda, 1), 3)
  actual_2 <- round(bootstrap_ci(resampled_lambdas, loo_lambdas, orig_lambda, 2), 3)
  actual_3 <- round(bootstrap_ci(resampled_lambdas, loo_lambdas, orig_lambda, 3), 3)
  actual_4 <- round(bootstrap_ci(resampled_lambdas, loo_lambdas, orig_lambda, 4), 3)
  
  expect_equal(object = actual_1, expected = expected_1, tolerance = 0.05)
  expect_equal(object = actual_2, expected = expected_2, tolerance = 0.05)
  expect_equal(object = actual_3, expected = expected_3, tolerance = 0.05)
  expect_equal(object = actual_4, expected = expected_4, tolerance = 0.05)
  
})
