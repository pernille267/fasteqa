library(testthat)
library(fasteqa)

# Required Functions for Testing
double_check_predict_eqa <- function(data, new_data, imprecision_estimates, R = 3, R_ratio = 1, method = "fg",
                                     level = 0.99, allow_reverse_regression = FALSE, rounding = 3){
  
  lambda <- imprecision_estimates$lambda
  Var_B <- imprecision_estimates$Var_B
  Var_A <- imprecision_estimates$Var_A
  
  n <- length(data$MP_B)
  
  mx <- mean(data$MP_B)
  my <- mean(data$MP_A)
  msxx <- var(data$MP_B)
  msyy <- var(data$MP_A)
  msxy <- var(data$MP_B, data$MP_A)
  
  sub_expression_1 <- msyy - lambda * msxx
  sub_expression_2 <- sqrt((msyy - lambda * msxx) ** 2 + 4 * lambda * msxy ** 2)
  sub_expression_3 <- 2 * msxy
  
  b1_deming <- (sub_expression_1 + sub_expression_2) / sub_expression_3
  b0_deming <- my - b1_deming * mx
  b1_ols <- msxy / msxx
  b0_ols <- my - b1_ols * mx
  
  var_b1_deming <- (b1_deming ** 2) / (n * msxy ** 2) * (msxx * msyy - msxy ** 2)
  Var_B_estimate <- (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda)
  Var_A_estimate <- lambda * Var_B_estimate
  
  latent <- data$MP_B * (lambda / (lambda + b1_deming ** 2)) + (b1_deming / (lambda + b1_deming ** 2)) * (data$MP_A - b0_deming)
  mu_hat <- mean(latent)
  
  predicted_MP_A_fg_clsi <- b0_deming + b1_deming * new_data$MP_B
  predicted_MP_A_ols <- b0_ols + b1_ols * new_data$MP_B
  predicted_MP_B_fg <- mu_hat * (1 - msxy / b1_deming / msxx) + predicted_MP_A_fg_clsi * msxy / b1_deming / msxx
  
  residuals_ols <- data$MP_A - (b0_ols + b1_ols * data$MP_B)
  mmse <- var(residuals_ols) * (n - 1) / (n - 2)
  
  t_quantile_fg_ols <- qt(1 - (1 - level) / 2, n - 2)
  t_quantile_clsi <- qt(1 - (1 - level) / 2, n * (R - 1))
  
  var_pred_error_fg <- (1 + 1 / (n - 2)) * (var_b1_deming * (predicted_MP_B_fg - mu_hat) ** 2 + var_b1_deming * Var_B_estimate * R_ratio + (1 + 1 / n) * (Var_A_estimate + (b1_deming ** 2) * Var_B_estimate) * R_ratio)
  var_pred_error_clsi <- var_b1_deming * (new_data$MP_B - mx) ** 2 + (1 + 1 / n) * (Var_A + (b1_deming ** 2) * Var_B) / R
  var_pred_error_ols <- mmse * (1 + 1 / n + (new_data$MP_B - mx) ** 2 / (msxx * (n - 1)))
  
  lwr_fg <- predicted_MP_A_fg_clsi - t_quantile_fg_ols * sqrt(var_pred_error_fg)
  upr_fg <- predicted_MP_A_fg_clsi + t_quantile_fg_ols * sqrt(var_pred_error_fg)
  lwr_clsi <- predicted_MP_A_fg_clsi - t_quantile_clsi * sqrt(var_pred_error_clsi)
  upr_clsi <- predicted_MP_A_fg_clsi + t_quantile_clsi * sqrt(var_pred_error_clsi)
  lwr_ols <- predicted_MP_A_ols - t_quantile_fg_ols * sqrt(var_pred_error_ols)
  upr_ols <- predicted_MP_A_ols + t_quantile_fg_ols * sqrt(var_pred_error_ols)
  
  if(method == "fg"){
    output <- list("SampleID" = new_data$SampleID,
                   "MP_B" = new_data$MP_B,
                   "MP_A" = new_data$MP_A,
                   "prediction" = predicted_MP_A_fg_clsi,
                   "lwr" = lwr_fg,
                   "upr" = upr_fg,
                   "inside" = ifelse(new_data$MP_A >= lwr_fg & new_data$MP_A <= upr_fg, 1, 0))
  }
  else if(method == "clsi"){
    output <- list("SampleID" = new_data$SampleID,
                   "MP_B" = new_data$MP_B,
                   "MP_A" = new_data$MP_A,
                   "prediction" = predicted_MP_A_fg_clsi,
                   "lwr" = lwr_clsi,
                   "upr" = upr_clsi,
                   "inside" = ifelse(new_data$MP_A >= lwr_clsi & new_data$MP_A <= upr_clsi, 1, 0))
  }
  else if(method == "ols"){
    output <- list("SampleID" = new_data$SampleID,
                   "MP_B" = new_data$MP_B,
                   "MP_A" = new_data$MP_A,
                   "prediction" = predicted_MP_A_ols,
                   "lwr" = lwr_ols,
                   "upr" = upr_ols,
                   "inside" = ifelse(new_data$MP_A >= lwr_ols & new_data$MP_A <= upr_ols, 1, 0))
  }
  else{
    output <- list("Output" = NA_character_)
  }
  
  return(output)
  
}

simulate_cs_data_eq_data <- function(parameters, b0 = 0, b1 = 1, obs_tau = 5){
  parameters <- c(parameters, list(b0 = b0, b1 = b1))
  test_cs_data <- sim_eqa_data(parameters = parameters,
                               type = 0,
                               AR = TRUE,
                               include_parameters = TRUE)
  test_cs_data$parameters$n <- 1
  if(is.null(obs_tau)){
    test_eq_data <- sim_eqa_data(parameters = c(list(b0 = b0, b1 = b1), test_cs_data$parameters),
                                 type = 0,
                                 include_parameters = FALSE)
  }
  else{
    test_eq_data <- sim_eqa_data(parameters = c(list(obs_tau = obs_tau, b0 = b0, b1 = b1),
                                                test_cs_data$parameters),
                                 type = 0,
                                 include_parameters = FALSE)  
  }
  test_eq_data$SampleID <- paste("EQAM", test_eq_data$SampleID)
  test_cs_data <- test_cs_data$simulated_data
  impr_cs_data <- global_precision_estimates(test_cs_data)
  test_cs_data_mor <- fun_of_replicates(test_cs_data)
  return(list(test_cs_data_mor, impr_cs_data, test_eq_data))
}

get_inside_pi <- function(parameters, b0 = 0, b1 = 1, obs_tau = 5){
  
  cs_n <- parameters$n
  
  sim_data <- simulate_cs_data_eq_data(parameters, b0, b1, obs_tau)
  
  actual_1 <- predict_eqa(data = sim_data[[1]],
                          new_data = sim_data[[3]],
                          imprecision_estimates = sim_data[[2]],
                          R = parameters$R,
                          R_ratio = 1,
                          method = "fg",
                          level = 0.95,
                          rounding = 12L)
  
  actual_2 <- predict_eqa(data = sim_data[[1]],
                          new_data = sim_data[[3]],
                          imprecision_estimates = sim_data[[2]],
                          R = parameters$R,
                          R_ratio = 1,
                          method = "clsi",
                          level = 0.95,
                          rounding = 12L)
  
  actual_3 <- predict_eqa(data = sim_data[[1]],
                          new_data = sim_data[[3]],
                          imprecision_estimates = sim_data[[2]],
                          R = parameters$R,
                          R_ratio = 1,
                          method = "ols",
                          level = 0.95,
                          rounding = 12L)
  
  return(c("fg" = actual_1$inside,
           "clsi" = actual_2$inside,
           "ols" = actual_3$inside))
}


# TEST 1 (Checking for unexpected programming errors ...)
test_that(desc = "Check for programming errors", code = {
  
  parameters <- list(n = 25, R = 3, cvx = 0.02, cvy = 0.01, cil = 2, ciu = 10)
  b0 <- 0
  b1 <- 1
  obs_tau <- NULL
  
  # Test if custom R functions give same results as Rcpp function
  test_data <- simulate_cs_data_eq_data(parameters, b0, b1, obs_tau)
  
  # From R function
  expected_1 <- double_check_predict_eqa(test_data[[1]], test_data[[3]], test_data[[2]],
                                         method = "fg", allow_reverse_regression = FALSE, level = 0.95)
  expected_2 <- double_check_predict_eqa(test_data[[1]], test_data[[3]], test_data[[2]],
                                         method = "clsi", allow_reverse_regression = FALSE, level = 0.95)
  expected_3 <- double_check_predict_eqa(test_data[[1]], test_data[[3]], test_data[[2]],
                                         method = "ols", allow_reverse_regression = FALSE, level = 0.95)
  
  # From Rcpp function
  actual_1 <- predict_eqa(test_data[[1]], test_data[[3]], test_data[[2]],
                          method = "fg", allow_reverse_regression = FALSE, level = 0.95, rounding = 12)
  actual_2 <- predict_eqa(test_data[[1]], test_data[[3]], test_data[[2]],
                          method = "clsi", allow_reverse_regression = FALSE, level = 0.95, rounding = 12)
  actual_3 <- predict_eqa(test_data[[1]], test_data[[3]], test_data[[2]],
                          method = "ols", allow_reverse_regression = FALSE, level = 0.95, rounding = 12)
  
  test_obj_1 <- sapply(X = 2:length(expected_1),
                       FUN = function(stat) abs(actual_1[[stat]] - expected_1[[stat]]) < 1e-4)
  test_obj_2 <- sapply(X = 2:length(expected_2),
                       FUN = function(stat) abs(actual_2[[stat]] - expected_2[[stat]]) < 1e-4)
  test_obj_3 <- sapply(X = 2:length(expected_3),
                       FUN = function(stat) abs(actual_3[[stat]] - expected_3[[stat]]) < 1e-4)
  
  expect_equal(test_obj_1, rep(TRUE, length(test_obj_1)))
  expect_equal(test_obj_2, rep(TRUE, length(test_obj_2)))
  expect_equal(test_obj_3, rep(TRUE, length(test_obj_3)))
  
})

# TEST 2 (Check if ACF coverage is suitable)
test_that(desc = "Check ACF coverage", code = {
  
  parameters <- list(n = 25, R = 3, cvx = 0.02, cvy = 0.01, cil = 2, ciu = 10)
  b0 <- 0
  b1 <- 1
  obs_tau <- NULL
  
  pi_insides <- t(replicate(n = 2e3,
                            expr = get_inside_pi(parameters, b0, b1, obs_tau),
                            simplify = TRUE))
  
  pi_acf <- apply(pi_insides, 2, function(x) round(mean(x) * 100, 3L))
  
  expect_true(max(abs(pi_acf - 95)) < 2, label = "Empirical PI Conf. level in [93, 97]")
  
})

# TEST 3 (Check if PCF coverage is suitable)
test_that(desc = "Check ACF coverage", code = {
  
  parameters <- list(n = 25, R = 3, cvx = 0.02, cvy = 0.01, cil = 2, ciu = 10)
  b0 <- 0
  b1 <- 1
  obs_tau <- 3
  
  pi_insides <- t(replicate(n = 2e3,
                            expr = get_inside_pi(parameters, b0, b1, obs_tau),
                            simplify = TRUE))
  
  pi_pcf <- apply(pi_insides, 2, function(x) round(mean(x) * 100, 3L))
  
  expect_true(max(abs(pi_pcf - 95)) < 2, label = "Empirical PI Conf. level in [93, 97]")
  
})




