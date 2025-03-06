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
