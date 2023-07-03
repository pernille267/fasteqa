#include <Rcpp.h>
using namespace Rcpp;

//' Estimate bias and skewness corrected bootstrap confidence intervals
//' 
//' @title Estimate bias and skewness corrected bootstrap confidence intervals
//' @name BCa_bootstrap_ci
//' 
//' @param bootstrapped_parameter_estimates A \code{numeric vector} storing bootstrapped estimates of the target parameter. A couple of missing values are allowed and will be removed if present.
//' @param jackknife_parameter_estimates A \code{numeric vector} containing Jackknife estimates of the target parameter. A couple of missing values are allowed and will be removed if present.
//' @param original_parameter_estimate A \code{double} representing the original point estimate of the target parameter. A missing value is not allowed.
//' @param level A \code{double} between 0 and 1, indicating the desired confidence level for the estimated BCa bootstrap confidence interval of the target parameter. Defaulting to \code{0.95}.
//' @param silence An \code{integer} controlling the level of progress reporting for debugging and analysis. Set to \code{-1} or \code{0} for detailed reports. Defaults to \code{1} to suppress output.
//'
//' @description This function estimates the bootstrap confidence intervals for a specified parameter using a Bias-Correction and Acceleration (BCa) approach. It requires precomputed Jackknife (leave-one-out) estimates of the parameter of interest. These can be generated by iterating over all unique SampleIDs in your dataset and employing the \code{leave_one_out()} function from this package for each unique ID.
//'
//' @details To accurately estimate the theoretical confidence interval for the parameter of interest, the original dataset must adequately represent the underlying distribution of the population parameter. Typically, larger datasets increase the likelihood of representativeness. To compute bootstrap confidence intervals independent of Jackknife estimates, consider using the \code{bootstrap_ci()} function instead.
//'
//' @return A \code{numeric vector} with two entries: the lower and upper bounds of the BCa bootstrap confidence interval, respectively.  
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
  
  NumericVector raw_bpe = clone(bootstrapped_parameter_estimates);
  NumericVector raw_jpe = clone(jackknife_parameter_estimates);
  NumericVector bpe;
  NumericVector jpe;
  
  int n = raw_bpe.size();
  int m = raw_jpe.size();
  int B = n;
  int J = m;
  
  for(int i = 0; i < n; ++i){
    bool is_na_bpe = ISNAN(raw_bpe[i]);
    if(is_na_bpe){
      B--;
      continue;
    }
    bpe.push_back(raw_bpe[i]);
  }
  
  NumericVector sorted_bpe = clone(bpe);
  std::sort(sorted_bpe.begin(), sorted_bpe.end());
  
  for(int i = 0; i < m; ++i){
    bool is_na_jpe = ISNAN(raw_jpe[i]);
    if(is_na_jpe){
      J--;
      continue;
    }
    jpe.push_front(raw_jpe[i]);
  }
  
  if(silence < 1){
    Rcout << B << " bootstrap estimates are non-missing out of " << n << " recorded values" << "\n";
    Rcout << J << " jackknife estimates are non-missing out of " << m << " recorded values" << "\n";
    Rcout << "-------------------------------------------" << "\n";
  }
  
  // Converting integer types to float types when needed in division!
  float B_float = static_cast<float>(B);
  
  NumericVector jackknife_bias_cubed(J);
  NumericVector jackknife_bias_squared(J);
  
  float alpha = 1 - level;
  float jackknife_mean = mean(jpe);
  float z0_kernel = 0;

  // Calculation of bias-correction z0
  for(int i = 0; i < n; ++i){
    if(sorted_bpe[i] < original_parameter_estimate){
      z0_kernel += 1.0 / B_float;
    }
  }
  
  float z0 = R::qnorm(z0_kernel, 0, 1, 1, 0);
  // Calculation of the skewness-correction constant a_jackknife, via jack knife (leave-one-out cross validation)
  for(int j = 0; j < J; ++j){
    jackknife_bias_cubed[j] = pow(jackknife_mean - jackknife_parameter_estimates[j], 3);
    jackknife_bias_squared[j] = pow(jackknife_mean - jackknife_parameter_estimates[j], 2);
  }
  float num = sum(jackknife_bias_cubed);
  float den = 6.0 * pow(sum(jackknife_bias_squared), 1.5);
  float a_jackknife = num / den;
  if(silence < 1){
    Rcout << "The bias-correction coefficient is estimated to be: " << z0 << "\n";
    Rcout << "The skewness-correction coefficient is estimated to be: " << a_jackknife << "\n";
    Rcout << "-------------------------------------------" << "\n";
  }
  float z_alpha_halves = R::qnorm(alpha * 0.5, 0, 1, 0, 0);
  float new_lower_quant_kernel = z0 + ((z0 + z_alpha_halves) / (1 - a_jackknife * (z0 + z_alpha_halves)));
  float new_upper_quant_kernel = z0 + ((z0 - z_alpha_halves) / (1 - a_jackknife * (z0 - z_alpha_halves)));
  float new_lower_quant = R::pnorm(new_lower_quant_kernel, 0, 1, 0, 0);
  float new_upper_quant = R::pnorm(new_upper_quant_kernel, 0, 1, 0, 0);
  int new_lower = std::min(std::max(static_cast<int>(ceil(new_lower_quant * B)) - 1, 0), B - 1);
  int new_upper = std::min(std::max(static_cast<int>(ceil(new_upper_quant * B)) - 1, 0), B - 1);
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
