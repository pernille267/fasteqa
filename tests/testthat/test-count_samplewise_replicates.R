library(fasteqa)
library(testthat)

test_that(desc = "Check summary = 'none'", code = {
  
  # Test 1: One replicate is missing for the response
  test_data_1 <- test_data
  test_data_1$MP_A[1] <- NA_real_
  expected_1 <- list(R_i = unname(tapply(X = !is.na(test_data_1$MP_A) & !is.na(test_data_1$MP_B),
                                         INDEX = test_data_1$SampleID,
                                         FUN = sum,
                                         simplify = TRUE)))
  expect_equal(object = count_samplewise_replicates(test_data_1, "none"),
               expected = expected_1,
               ignore_attr = TRUE)
  
  # Test 2: One replicate is missing from the response and
  #         One replicate missing for respones and one missing for predictor
  test_data_2 <- test_data_1
  test_data_1$MP_B[22] <- NA_real_
  test_data_1$MP_A[23] <- NA_real_
  expected_2 <- list(R_i = unname(tapply(X = !is.na(test_data_2$MP_A) & !is.na(test_data_2$MP_B),
                                         INDEX = test_data_2$SampleID,
                                         FUN = sum,
                                         simplify = TRUE)))
  expect_equal(object = count_samplewise_replicates(test_data_2, "none"),
               expected = expected_2,
               ignore_attr = TRUE)
  
  # Test 3: one sample have one missing replicate in the response &
  #         one sample have no valid replicates &
  #         one sample have only one valid replicate
  test_data_3 <- test_data_1
  test_data_3$MP_B[55:56] <- NA_real_
  test_data_3$MP_A[57] <- NA_real_
  test_data_3$MP_B[c(58, 60)] <- NA_real_
  test_data_3$MP_A[58] <- NA_real_
  expected_3 <- list(R_i = unname(tapply(X = !is.na(test_data_3$MP_A) & !is.na(test_data_3$MP_B),
                                         INDEX = test_data_3$SampleID,
                                         FUN = sum,
                                         simplify = TRUE)))
  expected_3$R_i <- expected_3$R_i[which(expected_3$R_i != 0)]
  expect_equal(object = count_samplewise_replicates(test_data_3, "none"),
               expected = expected_3,
               ignore_attr = TRUE)
  
})

test_that(desc = "Check summary = 'mean', 'floor', 'ceiling', 'round'", code = {
  
  # Test 1: 
  test_data_1 <- test_data
  test_data_1$MP_A[1] <- NA_real_
  test_data_1$MP_B[55:56] <- NA_real_
  test_data_1$MP_A[57] <- NA_real_
  test_data_1$MP_B[c(58, 60)] <- NA_real_
  test_data_1$MP_A[58] <- NA_real_
  
  expected <- tapply(X = !is.na(test_data_1$MP_A) & !is.na(test_data_1$MP_B),
                     INDEX = test_data_1$SampleID,
                     FUN = sum,
                     simplify = TRUE)
  expected <- expected[which(expected != 0)]
  expected_1 <- list(R_i = mean(expected))
  expected_2 <- list(R_i = floor(mean(expected)))
  expected_3 <- list(R_i = ceiling(mean(expected)))
  expected_4 <- list(R_i = round(mean(expected)))
  
  expect_equal(object = count_samplewise_replicates(test_data_1, "mean"),
               expected = expected_1,
               ignore_attr = TRUE)
  expect_equal(object = count_samplewise_replicates(test_data_1, "floor"),
               expected = expected_2,
               ignore_attr = TRUE)
  expect_equal(object = count_samplewise_replicates(test_data_1, "ceiling"),
               expected = expected_3,
               ignore_attr = TRUE)
  expect_equal(object = count_samplewise_replicates(test_data_1, "round"),
               expected = expected_4,
               ignore_attr = TRUE)
  
})

test_that(desc = "Check summary = 'median', 'mode'", code = {
  
  # Test 1: One sample with one replicate
  #         Nine samples with two replicates
  #         Ten samples with three replicates
  test_data_1 <- test_data
  test_data_1$MP_A[c(1, 23, 26, 29)] <- NA_real_
  test_data_1$MP_B[c(2, 4, 7, 10, 13, 16, 19)] <- NA_real_
  
  expected <- tapply(X = !is.na(test_data_1$MP_A) & !is.na(test_data_1$MP_B),
                     INDEX = test_data_1$SampleID,
                     FUN = sum,
                     simplify = TRUE)
  expected <- expected[which(expected != 0)]
  
  expected_1 <- list("R_i" = median(expected))
  expected_2 <- list("R_i" = as.numeric(names(table(expected)))[which.max(table(expected))])
  
  expect_equal(object = count_samplewise_replicates(test_data_1, "median"),
               expected = expected_1,
               ignore_attr = TRUE)
  expect_equal(object = count_samplewise_replicates(test_data_1, "mode"),
               expected = expected_2,
               ignore_attr = TRUE)
  
  # Test 2: One sample with one replicate
  #         Ten samples with two replicates
  #         Nine samples with three replicates
  test_data_2 <- test_data_1
  test_data_2$MP_B[32] <- NA_real_
  
  expected <- tapply(X = !is.na(test_data_2$MP_A) & !is.na(test_data_2$MP_B),
                     INDEX = test_data_2$SampleID,
                     FUN = sum,
                     simplify = TRUE)
  expected <- expected[which(expected != 0)]
  
  expected_1 <- list("R_i" = median(expected))
  expected_2 <- list("R_i" = as.numeric(names(table(expected)))[which.max(table(expected))])
  
  expect_equal(object = count_samplewise_replicates(test_data_2, "median"),
               expected = expected_1,
               ignore_attr = TRUE)
  expect_equal(object = count_samplewise_replicates(test_data_2, "mode"),
               expected = expected_2,
               ignore_attr = TRUE)
  
})

test_that(desc = "Test Error Handling", code = {
  
  # Check if the function returns list with NA if 'invalid_NA = TRUE'
  expect_equal(object = count_samplewise_replicates(data = "NA",
                                                    invalid_NA = TRUE),
               expected = list(R_i = NA_real_),
               ignore_attr = TRUE)
  expect_equal(object = count_samplewise_replicates(data = test_data,
                                                    summary = "something that is not allowed",
                                                    invalid_NA = TRUE),
               expected = list(R_i = NA_real_),
               ignore_attr = TRUE)
  # Check if the function throws error if 'invalid_NA = FALSE'
  expect_error(object = count_samplewise_replicates(data = "NA",
                                                    invalid_NA = FALSE),
               regexp = "data does not contain at least one of ")
  expect_error(object = count_samplewise_replicates(data = test_data,
                                                    summary = "something that is not allowed",
                                                    invalid_NA = FALSE),
               regexp = "'something that is not allowed' is not valid. Valid entries are:")
  
  
})

