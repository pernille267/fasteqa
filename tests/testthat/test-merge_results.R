library(fasteqa)
library(testthat)
library(data.table)

test_that(desc = "Testing if merging is successful", code = {
  
  # Generate some fake components
  ce_data <- list("comparison" = rep(c("A - B", "A - C", "B - C"), each = 2),
                  "SampleID" = rep(c("EQAM 1", "EQAM 2"), times = 3),
                  "MP_A" = c(4.146, 7.931, 4.146, 7.931, 4.193, 7.957),
                  "MP_B" = c(4.193, 7.957, 4.313, 8.048, 4.313, 8.048),
                  "prediction" = c(4.176, 7.944, 4.214, 7.998, 4.213, 7.988),
                  "lwr" = c(4.092, 7.865, 4.151, 7.935, 3.792, 7.567),
                  "upr" = c(4.260, 8.028, 4.261, 8.077, 4.641, 8.399),
                  "inside" = c(1L, 1L, 0L, 0L, 1L, 1L))
  
  pb_data <- list("comparison" = rep(c("A - B", "A - C", "B - C"), each = 3),
                  "predictor" = c(2, 5, 8, 2, 5, 8, 2, 5, 8),
                  "prediction" = c(2.011, 5.045, 8.121, 2.081, 5.129, 8.229, 2.045, 5.093, 8.199),
                  "lwr" = c(1.981, 5.015, 8.091, 2.040, 5.088, 8.188, 1.837, 4.885, 7.991),
                  "upr" = c(2.041, 5.075, 8.151, 2.122, 5.170, 8.270, 2.253, 5.301, 8.407))
  
  zeta_data <- list("comparison" = c("A - B", "A - C", "B - C"),
                    "zeta" = c(1.2112, 1.3316, 2.9299),
                    "lwr" = c(0.8739, 0.9734, 1.7830),
                    "upr" = c(1.6737, 2.041, 5.5723),
                    "zeta_critical" = c(2.25, 2.25, 2.25),
                    "zeta_conclusion" = c(0L, 0L, 1L))
  
  impr_data <- list("comparison" = c("A - B", "A - C", "B - C"),
                    "CV_A" = c(0.0123, 0.0123, 0.0219),
                    "CV_A_lwr" = c(0.0102, 0.0102, 0.0096),
                    "CV_A_upr" = c(0.0160, 0.0160, 0.0285),
                    "CV_B" = c(0.0219, 0.0135, 0.0135),
                    "CV_B_lwr" = c(0.0096, 0.0078, 0.0078),
                    "CV_B_upr" = c(0.0285, 0.0491, 0.0491),
                    "lambda" = c(0.3154, 0.8301, 2.6316),
                    "lambda_lwr" = c(0.2639, 0.6513, 1.3403),
                    "lambda_upr" = c(1.2991, 1.5648, 4.0025))
  
  
  # Convert to data.table
  setDT(ce_data); setDT(pb_data); setDT(zeta_data); setDT(impr_data)
  
  # Try without including imprecision_data
  
  # Check if the function runs without error
  expect_no_error(object = merge_results(pb_data = pb_data,
                                         ce_data = ce_data,
                                         zeta_data = zeta_data,
                                         imprecision_data = impr_data,
                                         include_imprecision_estimates = FALSE,
                                         rounding = 2L))
  
  # Actual output
  actual <- merge_results(pb_data = pb_data,
                          ce_data = ce_data,
                          zeta_data = zeta_data,
                          imprecision_data = impr_data,
                          include_imprecision_estimates = FALSE,
                          rounding = 4L)
  
  actual_merged_ce_data <- copy(actual$merged_ce_data) |> setDT()
  actual_merged_pb_data <- copy(actual$merged_pb_data) |> setDT()
  
  # Check dimensions
  expect_length(object = actual_merged_ce_data$comparison, 6)
  expect_length(object = actual_merged_pb_data$comparison, 9)
  expect_length(object = names(actual_merged_ce_data), 13)
  expect_length(object = names(actual_merged_pb_data), 10)
  
  # Check names
  expect_named(object = actual_merged_ce_data,
               expected = c("comparison", "SampleID", "zeta", "zeta_ci_lwr",
                            "zeta_ci_upr", "zeta_upper", "dins_conclusion",
                            "MP_B", "MP_A", "prediction", "pi_lwr",
                            "pi_upr", "pi_inside"),
               ignore.order = TRUE,
               ignore.case = FALSE)
  expect_named(object = actual_merged_pb_data,
               expected = c("comparison", "zeta", "zeta_ci_lwr",
                            "zeta_ci_upr", "zeta_upper", "dins_conclusion",
                            "predictor", "prediction", "pi_lwr", "pi_upr"),
               ignore.order = TRUE,
               ignore.case = FALSE)
  
  # Check if all necessary information is transferred to output
  
  # Check ce_data
  
  # Check if information from ce_data is kept
  expect_contains(object = actual_merged_ce_data$SampleID,
                  expected = ce_data$SampleID)
  expect_contains(object = actual_merged_ce_data$comparison,
                  expected = ce_data$comparison)
  expect_contains(object = actual_merged_ce_data$MP_B,
                  expected = ce_data$MP_B)
  expect_contains(object = actual_merged_ce_data$MP_A,
                  expected = ce_data$MP_A)
  expect_contains(object = actual_merged_ce_data$prediction,
                  expected = ce_data$prediction)
  expect_contains(object = actual_merged_ce_data$pi_lwr,
                  expected = ce_data$lwr)
  expect_contains(object = actual_merged_ce_data$pi_upr,
                  expected = ce_data$upr)
  expect_contains(object = actual_merged_ce_data$pi_inside,
                  expected = ce_data$inside)
  
  # Check if information from zeta_data is transferred
  expect_contains(object = actual_merged_ce_data$zeta,
                  expected = zeta_data$zeta)
  expect_contains(object = actual_merged_ce_data$zeta_ci_lwr,
                  expected = zeta_data$lwr)
  expect_contains(object = actual_merged_ce_data$zeta_ci_upr,
                  expected = zeta_data$upr)
  expect_contains(object = actual_merged_ce_data$zeta_upper,
                  expected = zeta_data$zeta_critical)
  expect_contains(object = actual_merged_ce_data$dins_conclusion,
                  expected = zeta_data$zeta_conclusion)
  
  # Check pb_data
  
  # Check if information from ce_data is kept
  expect_contains(object = actual_merged_pb_data$comparison,
                  expected = pb_data$comparison)
  expect_contains(object = actual_merged_pb_data$predictor,
                  expected = pb_data$predictor)
  expect_contains(object = actual_merged_pb_data$prediction,
                  expected = pb_data$prediction)
  expect_contains(object = actual_merged_pb_data$pi_lwr,
                  expected = pb_data$lwr)
  expect_contains(object = actual_merged_pb_data$pi_upr,
                  expected = pb_data$upr)
  
  # Check if information from zeta_data is transferred
  expect_contains(object = actual_merged_pb_data$zeta,
                  expected = zeta_data$zeta)
  expect_contains(object = actual_merged_pb_data$zeta_ci_lwr,
                  expected = zeta_data$lwr)
  expect_contains(object = actual_merged_pb_data$zeta_ci_upr,
                  expected = zeta_data$upr)
  expect_contains(object = actual_merged_pb_data$zeta_upper,
                  expected = zeta_data$zeta_critical)
  expect_contains(object = actual_merged_pb_data$dins_conclusion,
                  expected = zeta_data$zeta_conclusion)
  
  
  
})


