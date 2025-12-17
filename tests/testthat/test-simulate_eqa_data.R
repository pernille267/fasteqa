library(testthat)
library(fasteqa)

# Required functions for testing
calculate_cv <- function(data, MP_A = FALSE){
  if(MP_A){
    grand_mean <- mean(data$MP_A, na.rm = TRUE)
    samplewise_var <- tapply(X = data$MP_A,
                             INDEX = data$SampleID,
                             FUN = var,
                             na.rm = TRUE,
                             simplify = TRUE)
    pooled_var <- mean(samplewise_var, na.rm = TRUE)
    return(sqrt(pooled_var) / grand_mean)
  }
  grand_mean <- mean(data$MP_B, na.rm = TRUE)
  samplewise_var <- tapply(X = data$MP_B,
                           INDEX = data$SampleID,
                           FUN = var,
                           na.rm = TRUE,
                           simplify = TRUE)
  pooled_var <- mean(samplewise_var, na.rm = TRUE)
  return(sqrt(pooled_var) / grand_mean) 
}

get_random_generated_parameters <- function(){
  return(sim_eqa_data(parameters = list(b0 = 0),
                      type = 0,
                      AR = TRUE,
                      include_parameters = TRUE)$parameters)
}

test_that(desc = "Test if base parameters do what they are expected to do", {
  
  # Simulate data using sim_eqa_data()
  sim_data <- sim_eqa_data(parameters = list(n = 500,
                                             R = 100,
                                             cvx = 0.02,
                                             cvy = 0.01,
                                             cil = 2,
                                             ciu = 10,
                                             b0 = -0.1,
                                             b1 = 1.1),
                           type = 0L,
                           AR = TRUE,
                           include_parameters = FALSE)
  
  # Check if simulated data have the desired structure
  
  # 1. Check if output is a list
  expect_true(object = is.list(sim_data))
  
  # 2. Check if n = 500
  expect_length(object = unique(sim_data$SampleID),
                n = 500)
  
  # 3. Check if SampleID contains 1, 2, ..., 499, 500
  expect_setequal(object = unique(sim_data$SampleID),
                  expected = 1L:500L)
  
  # 4. Check if R = 100
  expect_length(object = unique(sim_data$ReplicateID),
                n = 100)
  
  # 5. Check if ReplicateID contains 1, 2, ..., 99, 100
  expect_setequal(object = unique(sim_data$ReplicateID),
                  expected = 1L:100L)
  
  
  # 6. Check if cvx = 0.02 [Acceptable interval: (0.018, 0.022)]
  expect_equal(object = calculate_cv(sim_data, FALSE),
               expected = 0.02,
               tolerance = 0.1)
  
  # 7. Check if cvy = 0.01 [Acceptable interval: (0.008, 0.012)]
  expect_equal(object = calculate_cv(sim_data, TRUE),
               expected = 0.01,
               tolerance = 0.2)
  
  # Calculate MOR data
  sim_data_mor <- fun_of_replicates(sim_data)
  
  # 8. Check if cil = 2 [Acceptable interval: (1.8, 2.2)]
  expect_equal(object = min(sim_data_mor$MP_B, na.rm = TRUE),
               expected = 2,
               tolerance = 0.1)
  
  # 9. Check if ciu = 10 [Acceptable interval: (9.8, 10.2)]
  expect_equal(object = max(sim_data_mor$MP_B, na.rm = TRUE),
               expected = 10,
               tolerance = 0.02)
  
  # Estimate beta_0 and beta_1 [Acceptable interval: (-0.12, -0.08)]
  b0b1 <- unname(coef(lm(sim_data_mor$MP_A ~ sim_data_mor$MP_B)))
  
  # 10. Check if beta_0 = -0.1
  expect_equal(object = b0b1[1],
               expected = -0.1,
               tolerance = 0.2)
  
  # 11. Check if beta_1 = 1.1 [Acceptable interval: (1.045, 1.155)]
  expect_equal(object = b0b1[2],
               expected = 1.1,
               tolerance = 0.05)
  
  
  
})

test_that(desc = "Test random generation of parameters", code = {
  
  rpg <- replicate(n = 1e4,
                   expr = as.data.frame(get_random_generated_parameters()),
                   simplify = FALSE)
  rpg_bound <- rpg[[1]]
  for(i in 2:1e4){
    rpg_bound <- rbind(rpg_bound, rpg[[i]])
  }
  
  # 1. Check generation of n
  expect_equal(object = max(rpg_bound$n, na.rm = TRUE),
               expected = 30)
  expect_equal(object = min(rpg_bound$n, na.rm = TRUE),
               expected = 20)
  
  # 2. Check generation of R
  expect_equal(object = max(rpg_bound$R, na.rm = TRUE),
               expected = 4)
  expect_equal(object = min(rpg_bound$R, na.rm = TRUE),
               expected = 2)
  expect_equal(object = sum(rpg_bound$R == 2),
               expected = 1e3,
               tolerance = 0.1)
  expect_equal(object = sum(rpg_bound$R == 3),
               expected = 8.5e3,
               tolerance = 0.01)
  expect_equal(object = sum(rpg_bound$R == 4),
               expected = 5e2,
               tolerance = 0.2)
  
  # 3. Check generation of cvx
  expect_equal(object = quantile(rpg_bound$cvx, probs = c(0.25, 0.50, 0.75)),
               expected = qbeta(c(0.25, 0.50, 0.75), 2, 5) / 10,
               tolerance = 0.1,
               ignore_attr = TRUE)
  
  # 4. Check generation of cvy
  expect_equal(object = quantile(rpg_bound$cvy, probs = c(0.25, 0.50, 0.75)),
               expected = qbeta(c(0.25, 0.50, 0.75), 2, 5) / 10,
               tolerance = 0.1,
               ignore_attr = TRUE)
  
  # 5. Check generation of multiplier = ciu / cil - 1
  multiplier <- rpg_bound$ciu / rpg_bound$cil  - 1
  expect_equal(object = quantile(multiplier, probs = c(0.25, 0.50, 0.75)),
               expected = qbeta(c(0.25, 0.50, 0.75), 0.78, 11) * 44,
               tolerance = 0.1,
               ignore_attr = TRUE)
  
  
  
})

test_that(desc = "Test dist parameter", code = {
  
  # Simulate data using dist = 'norm'
  sim_data_norm_dist <- sim_eqa_data(parameters = list(n = 1e5,
                                                       R = 5,
                                                       cvx = 0.00001,
                                                       cvy = 0.00001,
                                                       cil = 2,
                                                       ciu = 10,
                                                       b0 = 0,
                                                       b1 = 1,
                                                       dist = "norm"),
                           type = 0L,
                           AR = FALSE,
                           include_parameters = FALSE)
  
  # Check if cil = tau_0.01 and ciu = tau_0.99 (5% margin of error accepted)
  expect_equal(object = quantile(sim_data_norm_dist$MP_B, probs = c(0.01, 0.99)),
               expected = qnorm(c(0.01, 0.99), 6, 0.5 * 8 / qnorm(0.99)),
               tolerance = 0.05,
               ignore_attr = TRUE)
  
  # Check if other quantiles match (5% margin of error accepted)
  expect_equal(object = quantile(sim_data_norm_dist$MP_B, probs = c(0.25, 0.50, 0.75)),
               expected = qnorm(c(0.25, 0.50, 0.75), 6, 0.2149292 * 8),
               tolerance = 0.05,
               ignore_attr = TRUE)
  
  # Simulate data using dist = 'lst'
  sim_data_lst_dist <- sim_eqa_data(parameters = list(n = 1e5,
                                                      R = 5,
                                                      cvx = 0.00001,
                                                      cvy = 0.00001,
                                                      cil = 2,
                                                      ciu = 10,
                                                      b0 = 0,
                                                      b1 = 1,
                                                      dist = "lst"),
                                     type = 0L,
                                     AR = FALSE,
                                     include_parameters = FALSE)
  
  # Check if cil = tau_0.01 and ciu = tau_0.99 (5% margin of error accepted)
  expect_equal(object = quantile(sim_data_lst_dist$MP_B, probs = c(0.01, 0.99)),
               expected = 6 + 0.5 * 8 * (1 / qt(0.99, 5)) * qt(c(0.01, 0.99), 5),
               tolerance = 0.05,
               ignore_attr = TRUE)
  
  # Check if other quantiles match (5% margin of error accepted)
  expect_equal(object = quantile(sim_data_lst_dist$MP_B, probs = c(0.25, 0.50, 0.75)),
               expected = 6 + 0.5 * 8 * (1 / qt(0.99, 5)) * qt(c(0.25, 0.50, 0.75), 5),
               tolerance = 0.05,
               ignore_attr = TRUE)
  
  # Simulate data using dist = 'lnorm'
  sim_data_lnorm_dist <- sim_eqa_data(parameters = list(n = 1e5,
                                                        R = 5,
                                                        cvx = 0.00001,
                                                        cvy = 0.00001,
                                                        cil = 2,
                                                        ciu = 10,
                                                        b0 = 0,
                                                        b1 = 1,
                                                        dist = "lnorm"),
                                    type = 0L,
                                    AR = FALSE,
                                    include_parameters = FALSE)
  
  # Check if cil = tau_0.01 and ciu = tau_0.99 (5% margin of error accepted)
  mulog <- 0.5 * (log(10 * 2))
  sdlog <- 0.5 * log(10 / 2) / qnorm(0.99)
  expect_equal(object = quantile(sim_data_lnorm_dist$MP_B, probs = c(0.01, 0.99)),
               expected = qlnorm(c(0.01, 0.99), mulog, sdlog),
               tolerance = 0.05,
               ignore_attr = TRUE)
  
  # Check if other quantiles match (5% margin of error accepted)
  expect_equal(object = quantile(sim_data_lnorm_dist$MP_B, probs = c(0.25, 0.50, 0.75)),
               expected = qlnorm(c(0.25, 0.50, 0.75), mulog, sdlog),
               tolerance = 0.05,
               ignore_attr = TRUE)
  
})


# --- Tests for the improved sim_eqa_data---

test_that("Default parameters run without error", {
  expect_no_error(sim_eqa_data2())
})

test_that("Output structure is correct (AR=TRUE, default)", {
  n_def <- 25
  R_def <- 3
  result <- sim_eqa_data2(AR = TRUE)
  
  expect_true(is.list(result))
  expect_named(result, c("SampleID", "ReplicateID", "MP_A", "MP_B"))
  expect_length(result$SampleID, n_def * R_def)
  expect_length(result$ReplicateID, n_def * R_def)
  expect_length(result$MP_A, n_def * R_def)
  expect_length(result$MP_B, n_def * R_def)
  
  expect_type(result$SampleID, "integer")
  expect_type(result$ReplicateID, "integer")
  expect_type(result$MP_A, "double")
  expect_type(result$MP_B, "double")
})

test_that("Output structure is correct (AR=FALSE)", {
  n_test <- 10
  result <- sim_eqa_data2(parameters = list(n = n_test), AR = FALSE)
  
  expect_true(is.list(result))
  expect_named(result, c("SampleID", "MP_A", "MP_B")) # No ReplicateID
  expect_length(result$SampleID, n_test)
  expect_length(result$MP_A, n_test)
  expect_length(result$MP_B, n_test)
  
  expect_type(result$SampleID, "integer")
  expect_type(result$MP_A, "double")
  expect_type(result$MP_B, "double")
})

test_that("include_parameters = TRUE works", {
  result_noparams <- sim_eqa_data2(AR = TRUE, include_parameters = FALSE)
  result_params <- sim_eqa_data2(AR = TRUE, include_parameters = TRUE)
  
  expect_false("parameters" %in% names(result_noparams))
  expect_true(is.list(result_params))
  expect_named(result_params, c("simulated_data", "parameters"))
  expect_true(is.list(result_params$simulated_data))
  expect_true(is.list(result_params$parameters))
  expect_named(result_params$simulated_data, c("SampleID", "ReplicateID", "MP_A", "MP_B"))
  expect_true(length(result_params$parameters) > 5) # Check if parameters list is populated
})

test_that("obs_tau overrides n and dist", {
  obs_tau_vals <- c(5, 10, 15, 20)
  n_obs <- length(obs_tau_vals)
  R_test <- 2
  params <- list(obs_tau = obs_tau_vals, R = R_test, n = 100, dist = "norm") # n/dist should be ignored
  result <- sim_eqa_data2(parameters = params, AR = TRUE, include_parameters = TRUE)
  
  expect_equal(length(result$simulated_data$SampleID), n_obs * R_test)
  expect_equal(length(unique(result$simulated_data$SampleID)), n_obs)
  expect_equal(result$parameters$n, n_obs)
  expect_equal(result$parameters$dist, "observed")
  expect_true(result$parameters$tau_provided)
})

test_that("Different distributions ('dist') run", {
  n_test <- 5
  R_test <- 1
  params_base <- list(n = n_test, R = R_test, cil = 1, ciu = 10)
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(dist = "unif"))))
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(dist = "norm"))))
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(dist = "lst", df_tau = 5))))
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(dist = "lnorm"))))
  
  # Rough check for uniform bounds
  res_unif <- sim_eqa_data2(parameters = c(params_base, list(dist = "unif")))
  # Need internal tau for exact check, but measurements should be mostly within wider bounds
  expect_true(all(res_unif$MP_B > 0 & res_unif$MP_B < 15)) # Allowing for some error variance
})

test_that("Different error distributions ('error_dist') run", {
  params <- list(n = 5, R = 2, error_dist = "norm")
  expect_no_error(sim_eqa_data2(parameters = params))
  params$error_dist <- "lt"
  params$dfx <- 4
  params$dfy <- 6
  expect_no_error(sim_eqa_data2(parameters = params))
})

test_that("Bias parameters (b0, b1) have an effect", {
  params_base <- list(n = 50, R = 1, cvx = 0.001, cvy = 0.001, cil = 10, ciu = 20) # Low noise
  res_nobias <- sim_eqa_data2(parameters = params_base)
  res_bias_intercept <- sim_eqa_data2(parameters = c(params_base, list(b0 = 5)))
  res_bias_slope <- sim_eqa_data2(parameters = c(params_base, list(b1 = 1.5)))
  
  # Use AR=FALSE for easier comparison
  res_nobias_mor <- sim_eqa_data2(parameters = params_base, AR = FALSE)
  res_bias_intercept_mor <- sim_eqa_data2(parameters = c(params_base, list(b0 = 5)), AR = FALSE)
  res_bias_slope_mor <- sim_eqa_data2(parameters = c(params_base, list(b1 = 1.5)), AR = FALSE)
  
  # Check if intercept bias shifts the mean difference
  mean_diff_nobias <- mean(res_nobias_mor$MP_A - res_nobias_mor$MP_B, na.rm = TRUE)
  mean_diff_intercept <- mean(res_bias_intercept_mor$MP_A - res_bias_intercept_mor$MP_B, na.rm = TRUE)
  expect_gt(mean_diff_intercept, mean_diff_nobias + 4) # Expect diff near b0=5
  
  # Check if slope bias changes the ratio
  mean_ratio_nobias <- mean(res_nobias_mor$MP_A / res_nobias_mor$MP_B, na.rm = TRUE)
  mean_ratio_slope <- mean(res_bias_slope_mor$MP_A / res_bias_slope_mor$MP_B, na.rm = TRUE)
  expect_gt(mean_ratio_slope, mean_ratio_nobias + 0.4) # Expect ratio near b1=1.5
})


test_that("Non-selectivity parameters run", {
  params_base <- list(n = 150, R = 1, cil = 10, ciu = 100)
  # Prop
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(prop = 0.2, mmax = 5))))
  # Qran lower
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(qpos = 0, qran = 0.3, mmax = 5))))
  # Qran upper
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(qpos = 1, qran = 0.2, mmax = 4))))
  
  # Check if prop introduces larger variability in differences
  res_nobias <- sim_eqa_data2(parameters = c(params_base, list(cvx=0.01, cvy=0.01)), AR=FALSE)
  res_prop <- sim_eqa_data2(parameters = c(params_base, list(prop = 0.2, mmax = 10, cvx=0.01, cvy=0.01)), AR=FALSE)
  sd_diff_nobias <- sd(res_nobias$MP_A - res_nobias$MP_B, na.rm = TRUE)
  sd_diff_prop <- sd(res_prop$MP_A - res_prop$MP_B, na.rm = TRUE)
  # Expect sd of differences to increase significantly when mmax is large
  expect_gt(sd_diff_prop, sd_diff_nobias * 2) 
})

test_that("Heteroscedasticity parameters run", {
  params_base <- list(n = 20, R = 1, cil = 10, ciu = 100)
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(eta = 0.5, eta0 = 0.5))))
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(eta = 1.0, eta0 = 0.1))))
  expect_no_error(sim_eqa_data2(parameters = c(params_base, list(eta = 0.0, eta0 = 2.0)))) # Only base scaling
})

test_that("Custom function 'g' runs and overrides type/bias", {
  my_g <- function(tau) { tau * 2 + 5 } # Define a simple R function
  params <- list(n = 10, R = 1, g = my_g, b0 = 100, b1 = 100, type = 3) # b0, b1, type should be ignored
  
  expect_no_error(sim_eqa_data2(parameters = params, include_parameters = TRUE))
  result <- sim_eqa_data2(parameters = params, include_parameters = TRUE, AR = FALSE)
  
  # Check if type/bias parameters were marked as ignored/NA in output
  expect_true(is.na(result$parameters$type))
  # expect_true(is.na(result$parameters$b0)) # b0/b1 are still stored, but not used by g() path
  # expect_true(is.na(result$parameters$b1))
  expect_true(result$parameters$g_provided)
  
  # Check if output roughly follows the custom function
  # MP_A should be approx 2*MP_B + 5 (before noise)
  mean_transform <- mean(result$simulated_data$MP_A - (2 * result$simulated_data$MP_B), na.rm = TRUE)
  expect_true(abs(mean_transform - 5) < 1.0) # Allow some deviation due to noise/c0/c1 effects
})

test_that("md_method = 'mar' introduces NAs", {
  params <- list(n = 500, R = 3, md_method = "mar", mar_prob = 0.4)
  result <- sim_eqa_data2(parameters = params, AR = TRUE)
  
  expect_true(any(is.na(result$MP_A)))
  expect_true(any(is.na(result$MP_B)))
  
  prop_na_a <- mean(is.na(result$MP_A))
  prop_na_b <- mean(is.na(result$MP_B))
  expect_true(abs(prop_na_a - 0.4) < 0.05) # Check proportion is roughly correct
  expect_true(abs(prop_na_b - 0.4) < 0.05)
  
  # Check AR=FALSE averages correctly over NAs
  result_mor <- sim_eqa_data2(parameters = params, AR = FALSE)
  expect_equal(length(result_mor$MP_A), params$n)
  expect_true(any(is.na(result_mor$MP_A))) # Mean can still be NA if all reps are NA
})

test_that("md_method = 'mnar' introduces NAs below threshold", {
  thresh <- 5
  # Ensure some values *will* be below threshold
  params <- list(n = 100, R = 1, cil = 0, ciu = 10, cvx = 0.3, cvy = 0.3, 
                 md_method = "mnar", mnar_threshold = thresh)
  result <- sim_eqa_data2(parameters = params, AR = TRUE)
  
  expect_true(any(is.na(result$MP_A)))
  expect_true(any(is.na(result$MP_B)))
  
  # Check that non-NA values are >= threshold
  expect_true(all(result$MP_A[!is.na(result$MP_A)] >= thresh - 1e-9)) # Allow for float precision
  expect_true(all(result$MP_B[!is.na(result$MP_B)] >= thresh - 1e-9))
})

test_that("md_method = 'mnar0' introduces NAs below threshold (default 0)", {
  # Ensure some negative values are possible
  params <- list(n = 100, R = 1, cil = -5, ciu = 5, cvx = 0.5, cvy = 0.5,
                 md_method = "mnar0") # Threshold defaults to 0
  result <- sim_eqa_data2(parameters = params, AR = TRUE)
  
  expect_true(any(is.na(result$MP_A)))
  expect_true(any(is.na(result$MP_B)))
  
  # Check that non-NA values are >= 0
  expect_true(all(result$MP_A[!is.na(result$MP_A)] >= 0 - 1e-9)) # Allow for float precision
  expect_true(all(result$MP_B[!is.na(result$MP_B)] >= 0 - 1e-9))
  
  # Check with explicit positive threshold
  params$mnar_threshold <- 2.0
  result2 <- sim_eqa_data2(parameters = params, AR = TRUE)
  expect_true(any(is.na(result2$MP_A)))
  expect_true(all(result2$MP_A[!is.na(result2$MP_A)] >= 2.0 - 1e-9))
})

test_that("md_method = 'marmnar' combines methods", {
  thresh <- 3
  params <- list(n = 100, R = 5, cil = 0, ciu = 10, cvx = 0.3, cvy = 0.3,
                 md_method = "marmnar", mar_prob = 0.05, mnar_threshold = thresh)
  result <- sim_eqa_data2(parameters = params, AR = TRUE)
  
  expect_true(any(is.na(result$MP_A)))
  expect_true(any(is.na(result$MP_B)))
  
  # Should have slightly more NAs than just MAR or just MNAR (probabilistically)
  result_mar_only <- sim_eqa_data2(parameters = list(n=100, R=5, md_method="mar", mar_prob=0.05))
  prop_na_marmnar <- mean(is.na(result$MP_A))
  prop_na_mar <- mean(is.na(result_mar_only$MP_A))
  expect_gt(prop_na_marmnar, prop_na_mar - 0.01) # Should be at least MAR level (minus noise)
})

test_that("Negative values are clamped at 0 if not using mnar0", {
  # Force negative values
  params <- list(n = 50, R = 1, cil = -10, ciu = -5, cvx = 0.1, cvy = 0.1, md_method = "none")
  result <- sim_eqa_data2(parameters = params)
  expect_true(all(result$MP_A >= 0))
  expect_true(all(result$MP_B >= 0))
  
  # Check MAR also clamps remaining values
  params$md_method <- "mar"
  params$mar_prob <- 0.1
  result_mar <- sim_eqa_data2(parameters = params)
  expect_true(all(result_mar$MP_A[!is.na(result_mar$MP_A)] >= 0))
  expect_true(all(result_mar$MP_B[!is.na(result_mar$MP_B)] >= 0))
})

test_that("Invalid inputs trigger errors", {
  expect_error(sim_eqa_data2(parameters = list(R = 0)), "must be >= 1")
  expect_error(sim_eqa_data2(parameters = list(cvx = -0.1)), "must be non-negative")
  expect_error(sim_eqa_data2(parameters = list(cil = 10, ciu = 5)), "must be greater than cil")
  expect_error(sim_eqa_data2(parameters = list(dist = "lst", df_tau = 0)), "must be positive")
  expect_error(sim_eqa_data2(parameters = list(error_dist = "lt", dfx = -1)), "must be positive")
  expect_error(sim_eqa_data2(parameters = list(prop = 1.1)), "must be between 0 and 1")
  expect_error(sim_eqa_data2(parameters = list(qran = -0.1)), "must be between 0 and 1")
  expect_error(sim_eqa_data2(parameters = list(md_method = "mar", mar_prob = 2)), "must be between 0 and 1")
  expect_error(sim_eqa_data2(parameters = list(dist = "bogus")), "Unknown distribution")
  expect_error(sim_eqa_data2(parameters = list(obs_tau = numeric(0))), "obs_tau is empty")
  expect_error(sim_eqa_data2(parameters = list(dist = "lnorm", cil = -1)), "must be > 0 for dist = 'lnorm'")
})

