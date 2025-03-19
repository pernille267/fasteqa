#include <Rcpp.h>
using namespace Rcpp;

//' @title Calculate Sample Skewness of a Random Sample
//' @name skewness
//'
//' @param x A \code{numeric} vector. A random sample.
//' @param na_rm A \code{logical} value. If \code{TRUE}, \code{NA}-values are
//'        removed prior to calculation of the sample skewness.
//'
//' @description
//' Calculates the sample skewness of a random sample \code{x}, using the
//' adjusted Fisher-Pearson standardized moment coefficient.
//' 
//' @details
//' Calculates the sample skewness of a random sample
//' \eqn{\lbrace x_i \rbrace_{i=1}^{N}}. The sample skewness, attempts to
//' estimate the theoretical skewness, \eqn{\gamma}, by using
//' the following estimator (adjusted Fisher-Pearson standardized moment
//' coefficient):
//'
//' \eqn{\hat{\gamma} = \frac{\sqrt{N(N-1)}}{N-2} \cdot
//' \frac{\frac{1}{N}\sum_{i=1}^{N}(x_i - \overline{x})^3}
//' {\Big[\frac{1}{N}\sum_{i=1}^{N}(x_i - \overline{x})^2\Big]^{1.5}}}.
//'
//' Note: This estimator biased if not \eqn{x_i \sim \mathrm{N}(\mu, \sigma^2)}.
//'
//' @return A \code{double}. The calculated sample skewness.
//'
//' @examples
//' 
//' # Reproducibility
//' set.seed(99)
//' 
//' # Simulate 1 million observations from log-normal(0, 0.25)
//' y <- rlnorm(n = 1e6, meanlog = 0, sdlog = 0.25)
//' sample_skew <- skewness(y)
//' 
//' # Theoretical skew for log-normal(0, 0.25)
//' theoretical_skew <- (exp(0.25^2) + 2) * sqrt(exp(0.25^2) - 1)
//' 
//' cat("Theoretical skewness:", round(theoretical_skew), "\n",
//'     " Empirical skewness:", round(sample_skew), "\n")
//' 
//' 

// [[Rcpp::export]]
double skewness(NumericVector x, bool na_rm = true) {
 int n = x.size();
 NumericVector x_clean;
 
 // Handle NA values
 if (na_rm) {
   x_clean = x[!is_na(x)];
 } else {
   // If any NA present and na_rm is false, return NA
   if (any(is_na(x))) {
     return NA_REAL;
   }
   x_clean = x;
 }
 
 // Get new size after NA removal
 n = x_clean.size();
 
 // Check if we have enough data points
 if (n < 3) {
   return NA_REAL;
 }
 
 double mean_x = mean(x_clean);
 double s2 = 0.0, s3 = 0.0;
 
 for(int i = 0; i < n; ++i) {
   double diff = x_clean[i] - mean_x;
   s2 += diff * diff;
   s3 += diff * diff * diff;
 }
 
 s2 = s2 / n;
 s3 = s3 / n;
 
 // Calculate skewness
 double skew = (s3 / pow(s2, 1.5)) * sqrt((double)n * (n - 1)) / (n - 2);
 
 return skew;
}

//' @title Calculate Sample Excess Kurtosis of a Random Sample
//' @name kurtosis
//'
//' @param x A \code{numeric} vector. A random sample.
//' @param na_rm A \code{logical} value. If \code{TRUE}, \code{NA}-values are
//'        removed prior to calculation of the sample excess kurtosis.
//'        
//' @description
//' Calculates the sample excess kurtosis of a random sample \code{x}, using
//' the adjusted Fisher-Pearson standardized moment coefficient estimator.
//' 
//' @details
//' Calculates the sample excess kurtosis of a random sample
//' \eqn{\lbrace x_i \rbrace_{i=1}^{N}}. The sample excess kurtosis, attempts
//' to estimate the theoretical excess kurtosis, \eqn{\Kappa}, by using
//' the following estimator (adjusted Fisher-Pearson standardized moment
//' coefficient):
//'
//' \eqn{\hat{\Kappa} = \frac{N(N+1)}{(N-1)(N-2)(N-3)} \cdot \frac{\sum_{i=1}^{N}(x_i - \overline{x})^4}{s^4} - 3 \cdot \frac{(N-1)^2}{(N-2)(N-3)}}.
//'
//' Here, \eqn{s^2} is the unbiased sample variance.
//' 
//' 
//' Note: this estimator will be biased if not
//' \eqn{x_i \sim \mathrm{N}(\mu, \sigma^2)}. Moreover, \eqn{\hat{\Kappa}}
//' estimates the theoretical excess kurtosis and not the raw kurtosis. Add
//' \eqn{3} to the sample excess kurtosis to get the sample raw kurtosis.
//'
//' @return A \code{double}. The calculated sample excess kurtosis.
//'
//' @examples
//' # Reproducibility
//' set.seed(99)
//' 
//' # Simulate from one million observations from N(0, 1)
//' y <- rnorm(n = 1e6)
//' sample_excess_kurt <- kurtosis(y)
//' 
//' # Compare with theoretical excess kurtosis of N(0, 1)
//' cat("Theoretical excess kurtosis:", round(0, 3), "\n",
//'     " Empirical excess kurtosis:", round(sample_excess_kurt, 3), "\n")
//' 

// [[Rcpp::export]]
double kurtosis(NumericVector x, bool na_rm = true) {
 int n = x.size();
 NumericVector x_clean;
 
 // Handle NA values
 if (na_rm) {
   x_clean = x[!is_na(x)];
 } else {
   // If any NA present and na_rm is false, return NA
   if (any(is_na(x))) {
     return NA_REAL;
   }
   x_clean = x;
 }
 
 // Get new size after NA removal
 n = x_clean.size();
 
 // Check if we have enough data points
 if (n < 3) {
   return NA_REAL;
 }
 
 double mean_x = mean(x_clean);
 double s2 = 0.0, s4 = 0.0;
 
 for(int i = 0; i < n; ++i) {
   double diff = x_clean[i] - mean_x;
   s2 += diff * diff;
   s4 += diff * diff * diff * diff;
 }
 
 s2 = s2 / (n - 1);
 
 double fact1 = (double)n * (n + 1) / (n - 1) / (n - 2) / (n - 3);
 double fact2 = 3.0 * pow(n - 1, 2) / (n - 2) / (n - 3);
 
 // Calculate excess kurtosis
 double kurt = fact1 * (s4 / pow(s2, 2.0)) - fact2;
 
 return kurt;
}

//' @title Fit a Ordinary Least Squares Model
//' @name ols_regression
//'
//' @param data A \code{list} or \code{data.table}. Must contain:
//'             \itemize{
//'                 \item \code{MP_A: } A \code{numeric} vector. The means of
//'                 replicated measurements from IVD-MD \code{MP_A} (response).
//'                 \item \code{MP_B: } A \code{numeric} vector. The means of
//'                 replicated measurements from IVD-MD \code{MP_B} (predictor).
//'             }
//'        
//' @description
//' Fits a Ordinary Least Squares (OLS) regression model to \code{data}. 
//' 
//' @details
//' No details are needed. It is OLS.
//'
//' @return A \code{list}. Contains information on the fitted model.
//'
//' @examples
//' # Reproducibility
//' set.seed(99)
//' 
//' print(rbinom(1, 5, 0.5))
//' 
// [[Rcpp::export]]
List ols_regression(List data) {
  
  // Extract traning data
  NumericVector x = data["MP_B"];
  NumericVector y = data["MP_A"];
  
  // Calculate length of traning data
  int n = x.length();
  
  // Converting integer types to double types when needed in division!
  double n_double = static_cast<double>(n);
  
  // Standard Statistics
  double sxx = var(x) * (n_double - 1);
  double mx = mean(x);
  double my = mean(y);
  double sxy = 0;
  double mse = 0;
  for(int i = 0; i < n; ++i){
    sxy = sxy + (x[i] - mx) * (y[i] - my);
  }
  
  double b1 = sxy / sxx;
  double b0 = my - b1 * mx;
  
  NumericVector residuals(n);
  NumericVector fitted(n);
  
  for (int i = 0; i < n; ++i) {
    fitted[i] = b0 + b1 * x[i];
    residuals[i] = y[i] - fitted[i];
    mse += pow(residuals[i], 2.0) / (n_double - 2.0);
  }
  
  
  List output = List::create(Named("b0") = b0,
                             Named("b1") = b1,
                             Named("fitted") = fitted,
                             Named("residuals") = residuals,
                             Named("mse") = mse);
  
  return output;
  
}

//' @title Perform a Breusch-Pagan Test
//' @name pb_test
//'
//' @param data A \code{list} or \code{data.table}. Must contain:
//'             \itemize{
//'                 \item \code{MP_A: } A \code{numeric} vector. The means of
//'                 replicated measurements from IVD-MD \code{MP_A} (response).
//'                 \item \code{MP_B: } A \code{numeric} vector. The means of
//'                 replicated measurements from IVD-MD \code{MP_B} (predictor).
//'             }
//' @param koenker A \code{logical} value. If \code{TRUE} (default), the Koenker
//'                modification is applied to the make the test more robust
//'                when data is non-normal.
//'        
//' @description
//' Performs a Breusch Pagan test on variance heterogeneity on \code{data}. 
//' 
//' @details
//' The Breusch-Pagan test and its Koenker modification are both used to detect
//' heteroskedasticity in linear regression models.
//' 
//' The Koenker version, proposed by Koenker in 1981, modifies the test to make
//' it more robust when the data is non-normal. Note however, that the Koenker
//' version will have poorer power than the unmodified Breusch-Pagan test if
//' the data is close to normal.
//' 
//' Note: \code{NA}-values are silently removed prior to calculating the 
//' Breusch-Pagan test statistic.
//'
//' @return A \code{double}. The calculated p-value of the test.
//'
//' @examples
//' 
//' # Required packages
//' library(fasteqa)
//' 
//' # Reproducibility
//' set.seed(99)
//' 
//' # Simulate some example data
//' x <- runif(n = 50, min = 0, max = 1)
//' y <- 0.1 + 0.9 * x + rnorm(n = 50, mean = 0, sd = 0.05 * (x + 1))
//' data <- list(MP_A = y,
//'              MP_B = x)
//' 
//' # The output (rounded to 3L)
//' print(round(bp_test(data), 3L))
//' 
//' 
// [[Rcpp::export]]
double bp_test(List data, bool koenker = true) {
 
 // Check if data contains required vectors
 if (!data.containsElementNamed("MP_B") || !data.containsElementNamed("MP_A")) {
   stop("Data must contain both 'MP_A' and 'MP_B' vectors");
 }
 
 // Extract predictor values data
 NumericVector x = data["MP_B"];
 NumericVector y = data["MP_A"];
 
 // Check for empty vectors
 if (x.size() == 0 || y.size() == 0) {
   stop("Input vectors cannot be empty");
 }
 
 // Check for size mismatch
 if (x.size() != y.size()) {
   stop("MP_A and MP_B must have the same length");
 }
 
 // Get number of pairs
 int n = x.size();
 
 // Check for NA values
 bool has_na = false;
 for (int i = 0; i < n; ++i) {
   if (ISNAN(x[i]) || ISNAN(y[i])) {
     has_na = true;
     break;
   }
 }
 
 // Handle NA values if present
 if (has_na) {
   
   // Create vectors without NA values
   NumericVector x_clean, y_clean;
   for (int i = 0; i < n; ++i) {
     if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(y[i])) {
       x_clean.push_back(x[i]);
       y_clean.push_back(y[i]);
     }
   }
   
   // Check if we have enough data after removing NAs
   if (x_clean.size() < 5) {
     stop("Not enough non-NA observations to perform the test");
   }
   
   // Create new data list with clean values
   List clean_data = List::create(
     Named("MP_A") = y_clean,
     Named("MP_B") = x_clean
   );
   
   // Get original OLS fit with clean data
   List ols_fit = ols_regression(clean_data);
   
   // Update n to the new size
   n = x_clean.size();
   
   // Extract residuals
   NumericVector residuals = ols_fit["residuals"];
  
  // Empty vectors
  NumericVector sq_residuals(n); // Sqaured residuals from orig OLS model
  NumericVector w_hat(n); // Modified residuals used in fitting aux. OLS model (response)
  
  // Get relevant measures
  double mle_mse = 0.0;
  double koenker_mse = 0.0;
  double denominator = 0.0;
  double numerator = 0.0;
  
  for (int i = 0; i < n; ++i) {
    sq_residuals[i] = pow(residuals[i], 2.0);
    mle_mse += sq_residuals[i] / (n + 0.0);
  }
  
  for (int i = 0; i < n; ++i) {
    w_hat[i] = sq_residuals[i] - mle_mse;
    koenker_mse += pow(sq_residuals[i] - mle_mse, 2.0) / (n + 0.0);
  }
  
  if (koenker) {
    denominator += koenker_mse;
  }
  
  else {
    denominator += 2.0 * pow(mle_mse, 2.0);
  }
  
  // Fit aux model
  List model_frame = List::create(Named("MP_A") = w_hat,
                                  Named("MP_B") = x_clean);
  List aux_ols_fit = ols_regression(model_frame);
  
  NumericVector w_hat_hat = aux_ols_fit["fitted"]; // Fitted aux. OLS model values
  
  for (int i = 0; i < n; ++i) {
    numerator += pow(w_hat_hat[i], 2.0);
  }
  
  double bp_test_stat = numerator / denominator;
  double pval = R::pchisq(bp_test_stat, 1.0, 0, 0);
  
  return pval;
  
  
 }
 
 // Get original OLS fit with data
 List ols_fit = ols_regression(data);
   
 // Extract residuals
 NumericVector residuals = ols_fit["residuals"];
 
 // Empty vectors
 NumericVector sq_residuals(n); // Sqaured residuals from orig OLS model
 NumericVector w_hat(n); // Modified residuals used in fitting aux. OLS model (response)
 
 // Get relevant measures
 double mle_mse = 0.0;
 double koenker_mse = 0.0;
 double denominator = 0.0;
 double numerator = 0.0;
 
 for (int i = 0; i < n; ++i) {
   sq_residuals[i] = pow(residuals[i], 2.0);
   mle_mse += sq_residuals[i] / (n + 0.0);
 }
 
 for (int i = 0; i < n; ++i) {
   w_hat[i] = sq_residuals[i] - mle_mse;
   koenker_mse += pow(sq_residuals[i] - mle_mse, 2.0) / (n + 0.0);
 }
 
 if (koenker) {
   denominator += koenker_mse;
 }
 
 else {
   denominator += 2.0 * pow(mle_mse, 2.0);
 }
 
 // Fit aux model
 List model_frame = List::create(Named("MP_A") = w_hat,
                                 Named("MP_B") = x);
 List aux_ols_fit = ols_regression(model_frame);
 
 NumericVector w_hat_hat = aux_ols_fit["fitted"]; // Fitted aux. OLS model values
 
 for (int i = 0; i < n; ++i) {
   numerator += pow(w_hat_hat[i], 2.0);
 }
 
 double bp_test_stat = numerator / denominator;
 double pval = R::pchisq(bp_test_stat, 1.0, 0, 0);
 
 return pval;
 
}


