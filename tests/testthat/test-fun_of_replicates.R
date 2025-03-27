# Required functions
cv <- function(x, na.rm = FALSE){
  return(sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
}

test_that(desc = "Check if function gives the expected results", code = {
  
  # Set seed for reproducibility
  set.seed(99)
  test_data <- sim_eqa_data(parameters = list(n = 10, R = 5, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10),
                            type = 0L,
                            AR = TRUE,
                            include_parameters = FALSE)
  test_data$MP_B[c(1:5, 45:49)] <- NA_real_
  
  true_mean_of_replicates_B <- unname(tapply(test_data$MP_B,
                                             test_data$SampleID,
                                             FUN = mean, na.rm = TRUE,
                                             simplify = TRUE))
  true_mean_of_replicates_A <- unname(tapply(test_data$MP_A,
                                             test_data$SampleID,
                                             FUN = mean,
                                             na.rm = TRUE,
                                             simplify = TRUE))
  true_median_of_replicates_B <- unname(tapply(test_data$MP_B,
                                               test_data$SampleID,
                                               FUN = median,
                                               na.rm = TRUE,
                                               simplify = TRUE))
  true_median_of_replicates_A <- unname(tapply(test_data$MP_A,
                                               test_data$SampleID,
                                               FUN = median,
                                               na.rm = TRUE,
                                               simplify = TRUE))
  true_var_of_replicates_B <- unname(tapply(test_data$MP_B,
                                            test_data$SampleID,
                                            FUN = var,
                                            na.rm = TRUE,
                                            simplify = TRUE))
  true_var_of_replicates_A <- unname(tapply(test_data$MP_A,
                                            test_data$SampleID,
                                            FUN = var,
                                            na.rm = TRUE,
                                            simplify = TRUE))
  true_sd_of_replicates_B <- unname(tapply(test_data$MP_B,
                                           test_data$SampleID,
                                           FUN = sd,
                                           na.rm = TRUE,
                                           simplify = TRUE))
  true_sd_of_replicates_A <- unname(tapply(test_data$MP_A,
                                           test_data$SampleID,
                                           FUN = sd,
                                           na.rm = TRUE,
                                           simplify = TRUE))
  true_cv_of_replicates_B <- unname(tapply(test_data$MP_B,
                                           test_data$SampleID,
                                           FUN = cv,
                                           na.rm = TRUE,
                                           simplify = TRUE))
  true_cv_of_replicates_A <- unname(tapply(test_data$MP_A,
                                           test_data$SampleID,
                                           FUN = cv,
                                           na.rm = TRUE,
                                           simplify = TRUE))
  true_min_of_replicates_A <- unname(tapply(test_data$MP_A,
                                            test_data$SampleID,
                                            FUN = function(x) if(length(x[!is.na(x)]) == 0){NA_real_}else{min(x, na.rm = TRUE)},
                                            simplify = TRUE))
  true_min_of_replicates_B <- unname(tapply(test_data$MP_B,
                                            test_data$SampleID,
                                            FUN = function(x) if(length(x[!is.na(x)]) == 0){NA_real_}else{min(x, na.rm = TRUE)},
                                            simplify = TRUE))
  true_max_of_replicates_A <- unname(tapply(test_data$MP_A,
                                            test_data$SampleID,
                                            FUN = function(x) if(length(x[!is.na(x)]) == 0){NA_real_}else{max(x, na.rm = TRUE)},
                                            simplify = TRUE))
  true_max_of_replicates_B <- unname(tapply(test_data$MP_B,
                                            test_data$SampleID,
                                            FUN = function(x) if(length(x[!is.na(x)]) == 0){NA_real_}else{max(x, na.rm = TRUE)},
                                            simplify = TRUE))
  
  expected_mean_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                      MP_A = true_mean_of_replicates_A,
                                      MP_B = true_mean_of_replicates_B)
  expected_median_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                        MP_A = true_median_of_replicates_A,
                                        MP_B = true_median_of_replicates_B)
  expected_var_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                     MP_A = true_var_of_replicates_A,
                                     MP_B = true_var_of_replicates_B)
  expected_sd_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                    MP_A = true_sd_of_replicates_A,
                                    MP_B = true_sd_of_replicates_B)
  expected_cv_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                    MP_A = true_cv_of_replicates_A,
                                    MP_B = true_cv_of_replicates_B)
  expected_min_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                     MP_A = true_min_of_replicates_A,
                                     MP_B = true_min_of_replicates_B)
  expected_max_of_replicates <- list(SampleID = as.character(unique(test_data$SampleID)),
                                     MP_A = true_max_of_replicates_A,
                                     MP_B = true_max_of_replicates_B)
  
  expect_equal(object = fun_of_replicates(test_data, fun = "mean"),
               expected = expected_mean_of_replicates,
               ignore_attr = TRUE)
  expect_equal(object = fun_of_replicates(test_data, fun = "median"),
               expected = expected_median_of_replicates,
               ignore_attr = TRUE)
  expect_equal(object = fun_of_replicates(test_data, fun = "var"),
               expected = expected_var_of_replicates,
               ignore_attr = TRUE)
  expect_equal(object = fun_of_replicates(test_data, fun = "sd"),
               expected = expected_sd_of_replicates,
               ignore_attr = TRUE)
  expect_equal(object = fun_of_replicates(test_data, fun = "cv"),
               expected = expected_cv_of_replicates,
               ignore_attr = TRUE)
  expect_equal(object = fun_of_replicates(test_data, fun = "min"),
               expected = expected_min_of_replicates,
               ignore_attr = TRUE)
  expect_equal(object = fun_of_replicates(test_data, fun = "max"),
               expected = expected_max_of_replicates,
               ignore_attr = TRUE)
  
})


