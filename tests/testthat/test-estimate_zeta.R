library(fasteqa)
library(testthat)

# Required functions for testing
double_check_estimate_zeta <- function(data){
  SampleID <- data$SampleID
  MP_A <- data$MP_A
  MP_B <- data$MP_B
  var_MP_A <- tapply(X = MP_A, SampleID, function(x) if(length(x) >= 2){var(x, na.rm = TRUE)}else{NA_real_}, simplify = TRUE)
  var_MP_B <- tapply(X = MP_B, SampleID, function(x) if(length(x) >= 2){var(x, na.rm = TRUE)}else{NA_real_}, simplify = TRUE)
  pooled_var_MP_A <- mean(var_MP_A, na.rm = TRUE)
  pooled_var_MP_B <- mean(var_MP_B, na.rm = TRUE)
  if(is.na(pooled_var_MP_A) | is.na(pooled_var_MP_B)){
    return(list(zeta = NA_real_))
  }
  
  get_zeta <- function(x, y, pvarx, pvary){
    na_pairs <- is.na(x) | is.na(y)
    x <- x[!na_pairs]
    y <- y[!na_pairs]
    n <- length(x)
    mx <- mean(x)
    my <- mean(y)
    sxx <- var(x) * (n - 1)
    sxy <- var(x, y) * (n - 1)
    b1 <- sxy / sxx
    b0 <- my - b1 * mx
    mse <- 1 / (n - 2) * sum((y - (b0 + b1 * x))**2)
    evarpred <- mse * (n + 2) / n
    return(evarpred / (pvary + (b1**2) * pvarx))
  }
  
  if(pooled_var_MP_A < 0.5 * pooled_var_MP_B){
    return(list(zeta = get_zeta(MP_A, MP_B, pooled_var_MP_A, pooled_var_MP_B)))
  }
  else if (pooled_var_MP_A >= 0.5 * pooled_var_MP_B){
    return(list(zeta = get_zeta(MP_B, MP_A, pooled_var_MP_B, pooled_var_MP_A)))
  }
  return(list(zeta = NA_real_))
  
}


# Checking for Possible Programming Errors
test_that(desc = "Checking for Programming errors", code = {
  test_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10)
  test_data <- simulate_eqa_data(parameters = test_parameters)
  
  actual <- estimate_zeta_ols(test_data)
  expected <- double_check_estimate_zeta(test_data)
  
  expect_equal(actual$zeta, expected$zeta)
})

# Check handling of NA values
test_that(desc = "Checking Handling of NA values", code = {
  test_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10)
  test_data <- simulate_eqa_data(parameters = test_parameters)
  test_data$MP_A[c(1:2, 9)] <- NA_real_
  
  actual <- estimate_zeta_ols(test_data)
  expected <- double_check_estimate_zeta(test_data)
  
  expect_equal(actual$zeta, expected$zeta)
})









