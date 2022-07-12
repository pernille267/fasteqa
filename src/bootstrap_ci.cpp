#include <Rcpp.h>
using namespace Rcpp;

//' Estimate simple bootstrap confidence intervals
//' 
//' @title Estimate simple bootstrap confidence intervals
//' @name bootstrap_ci
//' 
//' @param bootstrapped_parameter_estimates Numeric vector - Bootstrapped estimates of the parameter of interest
//' @param original_parameter_estimate Double - Original point estimate of the parameter of interest
//' @param type Integer - Which type of bootstrap confidence interval is required
//' \itemize{
//'   \item{\code{1: }}{Standard normal bootstrap confidence interval}
//'   \item{\code{2: }}{Basic bootstrap confidence interval}
//'   \item{\code{3: }}{Percentile bootstrap confidence interval}
//' }
//' @param level Double - Required confidence level for the estimated confidence interval
//' @param silence Integer - Controls the progress reports outputted for debugging and further examination of the command. \code{silence = -1} or \code{silence = 0} signify that progress reports should be printed. Default is \code{silence = 1} which suppresses all printing
//'
//' @description Estimation of bootstrap confidence intervals for a parameter of interest
//'
//' @details To get an adequate estimate of the theoretical confidence interval, it is essential that original data represents the underlying distribution for the population parameter. More data typically corresponds with larger probability of our data being representative. For second-order accurate and transformation respective confidence intervals, you may employ \code{BCa_bootstrap_ci()} instead. Note that \code{BCa_bootstrap_ci()} require jack knife estimates of the parameter of interest as well as bootstrap estimates 
//'
//' @return A numeric vector. The first element of the vector is the lower bound of the confidence interval, and the second element of the vector is the upper bound of the confidence interval 
//'
//' @examples \dontrun{
//'   fictional_original_estimate <- 2.1
//'   fictional_bootstrap_estimates <- runif(n = 1e3, min = 0.5, max = 2.5)
//'   bootstrap_ci(fictional_bootstrap_estimates,
//'                fictional_original_estimate,
//'                type = 3,
//'                level = 0.95,
//'                silence = 1)
//' }



// [[Rcpp::export]]
NumericVector bootstrap_ci(NumericVector bootstrapped_parameter_estimates, float original_parameter_estimate = 0, int type = 3, float level = 0.95, int silence = 1){
  NumericVector sorted_bpe = sort_unique(bootstrapped_parameter_estimates);
  int n = sorted_bpe.size();
  if(silence == 0){
    Rcout << n << " bootstrap estimates are recorded" << "\n";
    Rcout << "-------------------------------------------" << "\n";
  }
  float alpha = 1 - level;
  int lwr = ceil(alpha * n / 2) - 1;
  int upr = n - lwr;
  float se_parameter = sd(sorted_bpe);
  if(silence == 0){
    Rcout << "The " << lwr << "th bootstrap estimate is used as lower quantile" "\n";  
    Rcout << "The " << upr << "th bootstrap estimate is used as upper quantile" "\n"; 
    Rcout << "-------------------------------------------" << "\n";
  }
  NumericVector bci(2);
  if(type == 1){
    float alpha_halved = alpha / 2;
    bci[0] = original_parameter_estimate - R::qnorm(alpha_halved, 0, 1, 0, 0) * se_parameter;
    bci[1] = original_parameter_estimate + R::qnorm(alpha_halved, 0, 1, 0, 0) * se_parameter;
    if(silence == 0){
      Rcout << "The standard normal bootstrap confidence interval method is used" << "\n";
      Rcout << "Output:" << "\n";
      Rcout << "-------------------------------------------" << "\n";
    }
  }
  else if(type == 2){
    bci[0] = 2 * original_parameter_estimate - sorted_bpe[upr];
    bci[1] = 2 * original_parameter_estimate - sorted_bpe[lwr];
    if(silence == 0){
      Rcout << "The basic bootstrap confidence interval method is used" << "\n";
      Rcout << "Output:" << "\n";
      Rcout << "-------------------------------------------" << "\n";
    }
  }
  else if(type == 3){
    bci[0] = sorted_bpe[lwr];
    bci[1] = sorted_bpe[upr];
    if(silence == 0){
      Rcout << "The percentile bootstrap confidence interval method is used" << "\n";
      Rcout << "Output:" << "\n";
      Rcout << " ------------------------------------------- " << "\n";
    }
  }
  return bci;
}
