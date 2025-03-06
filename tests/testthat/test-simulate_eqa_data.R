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

