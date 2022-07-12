#include <Rcpp.h>
using namespace Rcpp;

//' Estimate bias and skewness corrected bootstrap confidence intervals
//' 
//' @title Estimate bias and skewness corrected bootstrap confidence intervals
//' @name BCa_bootstrap_ci
//' 
//' @param bootstrapped_parameter_estimates Numeric vector - Bootstrapped estimates of the parameter of interest
//' @param jackknife_parameter_estimates Numeric vector - Jack knife estimates of the parameter of interest
//' @param original_parameter_estimate Double - Original point estimate of the parameter of interest
//' @param level Double - Required confidence level for the estimated confidence interval
//' @param silence Integer - Controls the progress reports outputted for debugging and further examination of the command. \code{silence = -1} or \code{silence = 0} signify that progress reports should be printed. Default is \code{silence = 1} which suppresses all printing
//'
//' @description Estimation of bootstrap confidence intervals for a parameter of interest using the BCa-variant
//'
//' @details To get an adequate estimate of the theoretical confidence interval, it is essential that original data represents the underlying distribution for the population parameter. More data typically corresponds with larger probability of our data being representative. For bootstrap confidence interval estimation independent of jack knife estimates, \code{bootstrap_ci()} may be used instead
//'
//' @return A numeric vector. The first element of the vector is the lower bound of the confidence interval, and the second element of the vector is the upper bound of the confidence interval 
//'
//' @examples \dontrun{
//'   fictional_original_estimate <- 2.1
//'   fictional_bootstrap_estimates <- runif(n = 1e3, min = 0.5, max = 2.5)
//'   fictional_jackknife_estimates <- runif(n = 1e1, min = 1.5, max = 2.5)
//'   BCa_bootstrap_ci(fictional_bootstrap_estimates,
//'                    fictional_jackknife_estimates,
//'                    fictional_original_estimate,
//'                    level = 0.95,
//'                    silence = 1)
//' }





// [[Rcpp::export]]
NumericVector BCa_bootstrap_ci(NumericVector bootstrapped_parameter_estimates, NumericVector jackknife_parameter_estimates, float original_parameter_estimate = 0, float level = 0.95, int silence = 1){
  NumericVector sorted_bpe = sort_unique(bootstrapped_parameter_estimates);
  int n = sorted_bpe.size();
  int m = jackknife_parameter_estimates.size();
  if(silence == 0){
    Rcout << n << " bootstrap estimates are recorded" << "\n";
    Rcout << m << " Jack knife estimates are recorded" << "\n";
    Rcout << "-------------------------------------------" << "\n";
  }
  NumericVector diffs_cubed(m);
  NumericVector diffs_squared(m);
  float alpha = 1 - level;
  float jackknife_mean = mean(jackknife_parameter_estimates);
  float z0_kernel = 0;
  float a_jackknife = 0;
  // Calculation of bias-correction z0
  for(int i = 0; i < n; ++i){
    if(sorted_bpe[i] < original_parameter_estimate){
      z0_kernel = z0_kernel + 1;
    }
  }
  float z0 = R::qnorm(z0_kernel / n, 0, 1, 1, 0);
  // Calculation of the skewness-correction constant a_jackknife, via jack knife (leave-one-out cross validation)
  for(int j = 0; j < m; ++j){
    diffs_cubed[j] = pow(jackknife_parameter_estimates[j] - jackknife_mean, 3);
    diffs_squared[j] = pow(jackknife_parameter_estimates[j] - jackknife_mean, 2);
  }
  float num = sum(diffs_cubed);
  float den = 6 * pow(sum(diffs_squared), 3/2);
  a_jackknife = a_jackknife + (num / den);
  if(silence == 0){
    Rcout << "The bias-correction coefficient is estimated to be: " << z0 << "\n";
    Rcout << "The skewness-correction coefficient is estimated to be: " << a_jackknife << "\n";
    Rcout << "-------------------------------------------" << "\n";
  }
  float z_alpha_halves = R::qnorm(alpha / 2, 0, 1, 0, 0);
  float new_lower_quant_kernel = z0 + ((z0 + z_alpha_halves) / (1 - a_jackknife * (z0 + z_alpha_halves)));
  float new_upper_quant_kernel = z0 + ((z0 - z_alpha_halves) / (1 - a_jackknife * (z0 - z_alpha_halves)));
  float new_lower_quant = R::pnorm(new_lower_quant_kernel, 0, 1, 0, 0);
  float new_upper_quant = R::pnorm(new_upper_quant_kernel, 0, 1, 0, 0);
  int new_lower = ceil(new_lower_quant * n) - 1;
  int new_upper = ceil(new_upper_quant * n);
  if(silence == 0){
    Rcout << "The " << new_lower << "th bootstrap estimate is used as lower quantile" "\n";  
    Rcout << "The " << new_upper << "th bootstrap estimate is used as upper quantile" "\n"; 
    Rcout << "-------------------------------------------" << "\n";  
  }
  NumericVector bci(2);
  bci[0] = sorted_bpe[new_lower];
  bci[1] = sorted_bpe[new_upper];
  return bci;
}
