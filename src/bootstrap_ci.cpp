#include <Rcpp.h>
using namespace Rcpp;

//' Estimate simple bootstrap confidence intervals
//' 
//' @title Estimate simple bootstrap confidence intervals
//' @name bootstrap_ci
//' 
//' @param bootstrapped_parameter_estimates A \code{numeric vector} storing bootstrapped estimates of the target parameter. A couple of missing values are allowed and will be removed if present.
//' @param original_parameter_estimate A \code{double} representing the original point estimate of the target parameter. A missing value is not allowed.
//' @param type An \code{integer} that refers to the type of bootstrap confidence interval to be estimated. Valid inputs include:
//' \itemize{
//'   \item{\code{1: }}{Standard normal bootstrap confidence interval}
//'   \item{\code{2: }}{Basic bootstrap confidence interval}
//'   \item{\code{3: }}{Percentile bootstrap confidence interval}
//' }
//' @param level A \code{double} between 0 and 1, indicating the desired confidence level for the estimated bootstrap confidence interval of the target parameter. Defaulting to \code{0.95}.
//' @param silence An \code{integer} controlling the level of progress reporting for debugging and analysis. Set to \code{-1} or \code{0} for detailed reports. Defaults to \code{1} to suppress output.
//'
//' @description Estimation of bootstrap confidence intervals for a parameter of interest based on \code{bootstrapped_parameter_estimates}, \code{original_parameter_estimate}, \code{type} and \code{level}. 
//'
//' @details To get an adequate estimate of the theoretical confidence interval for the target parameter, it is essential that the original data represents the underlying distribution for the population parameter. More data typically corresponds with larger probability of our data being representative. For second-order accurate and transformation respective confidence intervals, you may employ \code{BCa_bootstrap_ci()} instead. Note that \code{BCa_bootstrap_ci()} requires jack knife estimates of the parameter of interest as well as bootstrap estimates 
//'
//' @return A \code{numeric vector} with two entries: the lower and upper bounds of the bootstrap confidence interval, respectively.
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
  
  NumericVector raw_bpe = clone(bootstrapped_parameter_estimates);
  int n = raw_bpe.size();
  int B = n;
  NumericVector bpe;
  
  for(int i = 0; i < n; ++i){
    bool is_na_bpe = ISNAN(raw_bpe[i]);
    if(is_na_bpe){
      B--;
      continue;
    }
    bpe.push_back(raw_bpe[i]);
  }
  
  if(B <= 0 || B > n){
    stop("The effective number of bootstrap parameter estimates [%d] is either zero or negative. The effective number of bootstrap parameter estimates should an integer equal to or larger than 1.",
         B);
  }
  
  NumericVector sorted_bpe = clone(bpe);
  std::sort(sorted_bpe.begin(), sorted_bpe.end());
  
  if(silence < 1){
    Rcout << B << " bootstrap estimates are non-missing out of " << n << " recorded values" << "\n";
    Rcout << "-------------------------------------------" << "\n";
  }
  
  double lower_quant = (1.0 - level) * 0.5;
  double upper_quant = 1.0 - lower_quant;
  double ceil_rounding_tol = 1e-3;
  
  int lower = std::min(std::max(static_cast<int>(ceil(lower_quant * B - ceil_rounding_tol)) - 1, 0), B - 1);
  int upper = std::min(std::max(static_cast<int>(ceil(upper_quant * B - ceil_rounding_tol)) - 1, 0), B - 1);
  
  float se_parameter = sd(sorted_bpe);
  
  if(silence < 1){
    Rcout << "The " << lower + 1 << "th smallest bootstrap estimate is used as lower quantile" << "\n";  
    Rcout << "The " << upper + 1 << "th smallest bootstrap estimate is used as upper quantile" << "\n"; 
    Rcout << "-------------------------------------------" << "\n";
  }
  
  NumericVector bci(2);
  
  if(type == 1){
    bci[0] = original_parameter_estimate - R::qnorm(lower_quant, 0, 1, 0, 0) * se_parameter;
    bci[1] = original_parameter_estimate + R::qnorm(lower_quant, 0, 1, 0, 0) * se_parameter;
    if(silence < 1){
      Rcout << "The standard normal bootstrap confidence interval method is used" << "\n";
      Rcout << "Output:" << "\n";
      Rcout << "-------------------------------------------" << "\n";
    }
  }
  
  else if(type == 2){
    bci[0] = 2 * original_parameter_estimate - sorted_bpe[upper];
    bci[1] = 2 * original_parameter_estimate - sorted_bpe[lower];
    if(silence < 1){
      Rcout << "The basic bootstrap confidence interval method is used" << "\n";
      Rcout << "Output:" << "\n";
      Rcout << "-------------------------------------------" << "\n";
    }
  }
  
  else if(type == 3){
    bci[0] = sorted_bpe[lower];
    bci[1] = sorted_bpe[upper];
    if(silence < 1){
      Rcout << "The percentile bootstrap confidence interval method is used" << "\n";
      Rcout << "Output:" << "\n";
      Rcout << " ------------------------------------------- " << "\n";
    }
  }
  return bci;
}
