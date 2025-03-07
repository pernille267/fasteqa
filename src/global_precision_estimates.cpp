#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// (*) Helper function for estimating variance
double estimate_variance_gpe(const std::vector<double>& values) {
  double sum = 0.0, sq_sum = 0.0;
  const int n = values.size();
  
  for (double value : values) {
    sum += value;
    sq_sum += value * value;
  }
  
  double mean = sum / n;
  return (sq_sum - sum * mean) / (n - 1);
}

// (**) Helper function for obtaining unique elements maintaining their original order
CharacterVector unique_preserve_order_gpe(CharacterVector x) {
  std::unordered_set<std::string> seen;
  CharacterVector result;
  
  for(int i = 0; i < x.length(); i++) {
    std::string current = as<std::string>(x[i]);
    if(seen.find(current) == seen.end()) {
      seen.insert(current);
      result.push_back(x[i]);
    }
  }
  return result;
}

//' @title Calculate Repeatability Variability Components for a Particular IVD-MD Comparison
//' @name global_precision_estimates
//' 
//' @param data A \code{list} or \code{data.table}. Must contain the following
//'             variables:
//'             \itemize{
//'               \item \code{SampleID: } A \code{character} vector. The sample
//'                                       identifiers.
//'               \item \code{ReplicateID: } A \code{character} vector. The
//'                                          replicate measurement identifiers.
//'               \item \code{MP_A: } A \code{numeric} vector. The measurements
//'                                   from IVD-MD \code{MP_A} (response).
//'               \item \code{MP_B: } A \code{numeric} vector. The measurements
//'                                   from IVD-MD \code{MP_B} (predictor).
//'             }
//'             
//' @description
//' Calculate various repeatability components, such as variance, coefficient
//' of variability and the ratio of variances.
//' 
//' @details
//' Five statistics are estimated. The repeatability variances of \code{MP_A}
//' and \code{MP_B} are estimated by calculating the \eqn{n} (number of
//' unique values in \code{SampleID}) sample variances using the \code{numeric}
//' values in \code{MP_A} and \code{MP_B}. Then we take the mean of the sample
//' variances to obtain pooled estimated variances for the true repeatability
//' variances. We denote these estimated pooled variances by \code{Var_A} and
//' \code{Var_B}, respectively. Using the grand mean of \code{MP_A} and
//' \code{MP_B} we can then calculate three other statistics, \code{CV_A},
//' \code{CV_B} and \code{lambda}. Here is a summary the different estimated
//' statistics:
//' \itemize{
//'   \item \code{Var_A: } Pooled variance of all samplewise variances based on
//'                        IVD-MD \code{MP_A}. An estimator for
//'                        \eqn{\sigma_v^2}. Denote the estimator
//'                        \eqn{\hat{\sigma}_v^2}
//'   \item \code{Var_B: } Pooled variance of all samplewise variances based on
//'                        IVD-MD \code{MP_B}. An estimator for
//'                        \eqn{\sigma_h^2}. Denote the estimator
//'                        \eqn{\hat{\sigma}_h^2}
//'   \item \code{CV_A: } Estimated repeatbility coefficient of variation.
//'                       calculated by \eqn{\hat{\sigma}_v / \overline{y}},
//'                       where \eqn{\overline{y}} is the grand sample mean of
//'                       the values in \code{MP_A}.
//'   \item \code{CV_B: } Estimated repeatbility coefficient of variation.
//'                       calculated by \eqn{\hat{\sigma}_h / \overline{x}},
//'                       where \eqn{\overline{x}} is the grand sample mean of
//'                       the values in \code{MP_B}.
//'   \item \code{lambda: } Estimated repeatbility variance ratio. 
//'                         calculated by
//'                         \eqn{\hat{\sigma}_v^2 / \hat{\sigma}_h^2}.
//' }
//' 
//' By default, \code{CV_A} and \code{CV_B} are represented as a decimal
//' number. These values can also be represented as percentages, and to
//' covert to percentages, one should multiply their raw values with \code{100}.
//' 
//' Note: If one uses log-transformed \code{data}, the interpretation of
//' \code{CV_A} and \code{CV_B} may change, depending on the application.
//' In the log-transformation case, the square-root of \code{Var_A} and
//' \code{Var_B} have a similar interpretation as \code{CV_A} and \code{CV_B}
//' calculated using raw \code{data}.
//' 
//'
//' @return
//' A \code{list} of length five. Each element contains the an estimated
//' statistic. See details for information on each of them.
//' 
//' @examples
//' library(fasteqa)
//' # Calculate imprecision estimates
//' repeatability_statistics <- global_precision_estimates(test_data)
//' 
//' # Output
//' print(repeatability_statistics)
//' 
//' # Convert CV_A and CV_B to percentages
//' repeatability_statistics$CV_A <- repeatability_statistics$CV_A * 100
//' repeatability_statistics$CV_B <- repeatability_statistics$CV_B * 100
//' 
//' # Round results (two decimals)
//' repeatability_statistics <- lapply(X = repeatability_statistics,
//'                                    FUN = round,
//'                                    digits = 2L)
//' 
//' # Convert to data.frame
//' repeatability_statistics <- as.data.frame(repeatability_statistics)
//' 
//' print(repeatability_statistics)
//' 


// [[Rcpp::export]]
List global_precision_estimates(List data) {
  
  // Extract data columns
  CharacterVector SampleID = data["SampleID"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  
  // Get unique samples and initialize result vectors
  CharacterVector summary_SampleID = unique_preserve_order_gpe(SampleID);
  const int n = summary_SampleID.size();
  const int N = SampleID.size();
  NumericVector ith_var_MP_A(n, NA_REAL);
  NumericVector ith_var_MP_B(n, NA_REAL);
  const int replicate_number_requirement = 2;
  
  // Create hash map for faster sample lookup
  std::unordered_map<String, std::vector<int>> sample_indices;
  for (int i = 0; i < N; ++i) {
    sample_indices[String(SampleID[i])].push_back(i);
  }
  
  // Calculate variances for each sample
  for (int i = 0; i < n; ++i) {
    String current_sample = String(summary_SampleID[i]);
    const std::vector<int>& indices = sample_indices[current_sample];
    
    std::vector<double> valid_measurements_A;
    std::vector<double> valid_measurements_B;
    valid_measurements_A.reserve(indices.size());
    valid_measurements_B.reserve(indices.size());
    
    // Collect valid measurements
    for (int idx : indices) {
      if (!ISNAN(MP_A[idx])) {
        valid_measurements_A.push_back(MP_A[idx]);
      }
      if (!ISNAN(MP_B[idx])) {
        valid_measurements_B.push_back(MP_B[idx]);
      }
    }
    
    // Calculate variances if enough replicates
    if (valid_measurements_A.size() >= replicate_number_requirement) {
      ith_var_MP_A[i] = estimate_variance_gpe(valid_measurements_A);
    }
    if (valid_measurements_B.size() >= replicate_number_requirement) {
      ith_var_MP_B[i] = estimate_variance_gpe(valid_measurements_B);
    }
  }
  
  // Calculate global statistics
  double var_MP_A = 0, var_MP_B = 0;
  int effective_n_A = 0, effective_n_B = 0;
  
  for (int i = 0; i < n; ++i) {
    if (!ISNAN(ith_var_MP_A[i])) {
      var_MP_A += ith_var_MP_A[i];
      effective_n_A++;
    }
    if (!ISNAN(ith_var_MP_B[i])) {
      var_MP_B += ith_var_MP_B[i];
      effective_n_B++;
    }
  }
  
  var_MP_A = effective_n_A > 0 ? var_MP_A / effective_n_A : NA_REAL;
  var_MP_B = effective_n_B > 0 ? var_MP_B / effective_n_B : NA_REAL;
  var_MP_A = var_MP_A > 0 ? var_MP_A : NA_REAL;
  var_MP_B = var_MP_B > 0 ? var_MP_B : NA_REAL;
  
  
  // Calculate means
  double mean_MP_A = 0, mean_MP_B = 0;
  int valid_count_A = 0, valid_count_B = 0;
  
  for (int i = 0; i < N; ++i) {
    if (!ISNAN(MP_A[i])) {
      mean_MP_A += MP_A[i];
      valid_count_A++;
    }
    if (!ISNAN(MP_B[i])) {
      mean_MP_B += MP_B[i];
      valid_count_B++;
    }
  }
  
  mean_MP_A = valid_count_A > 0 ? mean_MP_A / valid_count_A : NA_REAL;
  mean_MP_B = valid_count_B > 0 ? mean_MP_B / valid_count_B : NA_REAL;
  
  // Calculate final statistics
  double cv_MP_A = (!ISNAN(var_MP_A) && !ISNAN(mean_MP_A) && mean_MP_A != 0) ?
  sqrt(var_MP_A) / mean_MP_A : NA_REAL;
  double cv_MP_B = (!ISNAN(var_MP_B) && !ISNAN(mean_MP_B) && mean_MP_B != 0) ?
  sqrt(var_MP_B) / mean_MP_B : NA_REAL;
  double lambda = (!ISNAN(var_MP_A) && !ISNAN(var_MP_B) && var_MP_B != 0) ?
  var_MP_A / var_MP_B : NA_REAL;
  
  return List::create(
    Named("Var_A") = var_MP_A,
    Named("Var_B") = var_MP_B,
    Named("CV_A") = cv_MP_A,
    Named("CV_B") = cv_MP_B,
    Named("lambda") = lambda
  );
}

