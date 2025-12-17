#include <Rcpp.h>
#include <vector>        // Use std::vector for temporary storage
#include <numeric>       // For std::accumulate
#include <cmath>         // For std::sqrt, std::pow, ISNAN
#include <unordered_map> // For efficient grouping
#include <limits>        // For infinity or quiet NaN if needed
using namespace Rcpp;

// (*) Helper function for estimating variance
double estimate_variance(const std::vector<double>& values) {
  const int n = values.size();
  if (n < 2) {
    return NA_REAL;
  }
  double sum = 0.0, sq_sum = 0.0;
  
  for (double value : values) {
    sum += value;
    sq_sum += value * value;
  }
  
  double mean = sum / n;
  return (sq_sum - sum * mean) / (n - 1);
}

// (**) Helper function for obtaining unique elements maintaining their original order
CharacterVector unique_preserve_order(CharacterVector x) {
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

// (***) Helper function for obtaining estimated imprecision statistics
List global_precision_estimates_local(List data) {
  
  // Extract data columns
  CharacterVector SampleID = data["SampleID"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  
  // Get unique samples and initialize result vectors
  CharacterVector summary_SampleID = unique_preserve_order(SampleID);
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
      ith_var_MP_A[i] = estimate_variance(valid_measurements_A);
    }
    if (valid_measurements_B.size() >= replicate_number_requirement) {
      ith_var_MP_B[i] = estimate_variance(valid_measurements_B);
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

// Helper to calculate sum and count for mean
struct SumCount {
  double sum = 0.0;
  int count = 0;
};

// (****) Helper function for obtaining estimated zeta values
double estimate_zeta_local(List data) {
  
  // --- 1. Data Extraction and Preparation ---
  CharacterVector SampleID = data["SampleID"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  const int N = SampleID.size(); // Total number of rows
  
  // Basic input checks
  if (N == 0) {
    Rcpp::warning("Input data has zero rows.");
    return NA_REAL;
  }
  if (MP_A.size() != N || MP_B.size() != N) {
    stop("Input vectors (SampleID, MP_A, MP_B) must have the same length.");
  }
  
  // Build index map for efficient grouping by SampleID
  std::unordered_map<String, std::vector<int>> sample_indices;
  for (int i = 0; i < N; ++i) {
    sample_indices[String(SampleID[i])].push_back(i);
  }
  int n = sample_indices.size(); // Number of unique samples
  
  // --- 2. Calculate Sample Variances (Robustly) ---
  std::vector<double> sample_vars_A; // Store only valid variances
  std::vector<double> sample_vars_B;
  sample_vars_A.reserve(n); // Pre-allocate potential space
  sample_vars_B.reserve(n);
  const int replicate_number_requirement = 2;
  
  for (auto const& [sample_str, indices] : sample_indices) {
    std::vector<double> valid_measurements_A;
    std::vector<double> valid_measurements_B;
    valid_measurements_A.reserve(indices.size());
    valid_measurements_B.reserve(indices.size());
    
    // Collect valid (non-NA) measurements for this sample
    for (int idx : indices) {
      if (!ISNAN(MP_A[idx])) {
        valid_measurements_A.push_back(MP_A[idx]);
      }
      if (!ISNAN(MP_B[idx])) {
        valid_measurements_B.push_back(MP_B[idx]);
      }
    }
    
    // Calculate variance if enough valid replicates exist
    if (valid_measurements_A.size() >= replicate_number_requirement) {
      sample_vars_A.push_back(estimate_variance(valid_measurements_A));
    }
    if (valid_measurements_B.size() >= replicate_number_requirement) {
      sample_vars_B.push_back(estimate_variance(valid_measurements_B));
    }
  }
  
  // --- 3. Pool Variances & Calculate Lambda ---
  double var_MP_A = NA_REAL;
  double var_MP_B = NA_REAL;
  int effective_n_A = sample_vars_A.size();
  int effective_n_B = sample_vars_B.size();
  
  if (effective_n_A > 0) {
    // Remove potential NAs returned by estimate_variance if n<2 (if applicable)
    sample_vars_A.erase(std::remove_if(sample_vars_A.begin(), sample_vars_A.end(),
                                       [](double v){ return ISNAN(v); }),
                                       sample_vars_A.end());
    effective_n_A = sample_vars_A.size(); // Update count after potential removal
    if (effective_n_A > 0) {
      var_MP_A = std::accumulate(sample_vars_A.begin(), sample_vars_A.end(), 0.0) / effective_n_A;
    }
  }
  if (effective_n_B > 0) {
    sample_vars_B.erase(std::remove_if(sample_vars_B.begin(), sample_vars_B.end(),
                                       [](double v){ return ISNAN(v); }),
                                       sample_vars_B.end());
    effective_n_B = sample_vars_B.size();
    if (effective_n_B > 0) {
      var_MP_B = std::accumulate(sample_vars_B.begin(), sample_vars_B.end(), 0.0) / effective_n_B;
    }
  }
  
  if (ISNAN(var_MP_A) || ISNAN(var_MP_B)) {
    Rcpp::warning("Could not estimate pooled variance for MP_A and/or MP_B (insufficient replicates or all NA).");
    return NA_REAL;
  }
  
  // Handle near-zero variance for B
  const double tolerance = 1e-9; // Use a small tolerance
  if (std::abs(var_MP_B) < tolerance) {
    if (var_MP_B == 0 && var_MP_A == 0) {
      // Special case: If both are zero, lambda is undefined but maybe treat as 1? Or NA?
      // Current code implies lambda -> Inf or 0. Let's return NA as safest.
      Rcpp::warning("Estimated variance for MP_B is zero (or near zero). Cannot calculate lambda reliably.");
      return NA_REAL;
    }
    if (var_MP_B == 0 && var_MP_A != 0) {
      Rcpp::warning("Estimated variance for MP_B is zero, but MP_A variance is non-zero. Lambda -> Inf.");
      // Depending on context, maybe return Inf or NA? Returning NA is safer.
      return NA_REAL;
    }
    // If var_MP_B is just very small but non-zero:
    Rcpp::warning("The estimate of Var_MP_B is very small. Consider results carefully.");
    // Proceed with calculation, lambda might be huge. Original code replaced with 0.0000001, let's avoid modification.
  }
  if (var_MP_A < 0 || var_MP_B < 0) {
    // This shouldn't happen with correct variance calculation but good check
    Rcpp::warning("Estimated variance negative (%f, %f), check variance calculation.", var_MP_A, var_MP_B);
    return NA_REAL;
  }
  
  // Check for division by zero before calculating lambda
  if (var_MP_B == 0.0) {
    Rcpp::warning("Division by zero: var_MP_B is exactly zero.");
    return NA_REAL;
  }
  double lambda = var_MP_A / var_MP_B;
  
  
  // --- 4. Calculate Overall Means and Prepare for Regression ---
  std::vector<double> x_clean, y_clean;
  x_clean.reserve(N);
  y_clean.reserve(N);
  SumCount mean_calc_A = {0.0, 0};
  SumCount mean_calc_B = {0.0, 0};
  
  bool swap_xy = (lambda >= 0.5); // Determine which variable is X
  
  for (int i = 0; i < N; ++i) {
    bool na_A = ISNAN(MP_A[i]);
    bool na_B = ISNAN(MP_B[i]);
    
    // Only consider pairs where *both* are non-NA for regression
    if (!na_A && !na_B) {
      if (swap_xy) { // B is X, A is Y
        x_clean.push_back(MP_B[i]);
        y_clean.push_back(MP_A[i]);
        mean_calc_B.sum += MP_B[i]; mean_calc_B.count++;
        mean_calc_A.sum += MP_A[i]; mean_calc_A.count++; // Still need overall mean of A
      } else { // A is X, B is Y
        x_clean.push_back(MP_A[i]);
        y_clean.push_back(MP_B[i]);
        mean_calc_A.sum += MP_A[i]; mean_calc_A.count++;
        mean_calc_B.sum += MP_B[i]; mean_calc_B.count++; // Still need overall mean of B
      }
    }
    // Accumulate for overall means even if one is NA (matches original logic, but could be debated)
    // Update: Sticking to only non-NA pairs for means used in regression for consistency.
    // If truly *overall* means were needed including single NAs, calculate separately.
    // The original code accumulated means *only* when both were non-NA in the lambda>=0.5 branch,
    // but seemed to accumulate overall means before that. Let's be consistent and use
    // means *of the data used in the regression*.
  }
  
  int effective_N_reg = x_clean.size(); // Number of valid pairs for regression
  
  if (effective_N_reg < 3) { // Need at least 3 points for OLS (2 for slope, 1 for error df)
    Rcpp::warning("Insufficient valid data pairs (< 3) for regression.");
    return NA_REAL;
  }
  
  // Calculate means of the data used for regression
  double mx = (mean_calc_B.count > 0) ? mean_calc_B.sum / mean_calc_B.count : NA_REAL; // Mean of X (if B is X)
  double my = (mean_calc_A.count > 0) ? mean_calc_A.sum / mean_calc_A.count : NA_REAL; // Mean of Y (if A is Y)
  if (!swap_xy) { // If A is X, swap the calculated means
    std::swap(mx, my);
  }
  
  if (ISNAN(mx) || ISNAN(my)) {
    Rcpp::warning("Could not calculate means for regression data.");
    return NA_REAL;
  }
  
  // --- 5. Perform OLS Regression ---
  double sxx = 0.0, sxy = 0.0, sse = 0.0;
  for (int i = 0; i < effective_N_reg; ++i) {
    sxx += std::pow(x_clean[i] - mx, 2);
    sxy += (x_clean[i] - mx) * (y_clean[i] - my);
  }
  
  if (std::abs(sxx) < tolerance) {
    Rcpp::warning("Sum of squares for X (sxx) is zero or near zero. Cannot calculate regression slope.");
    return NA_REAL;
  }
  
  double b1 = sxy / sxx;
  double b0 = my - b1 * mx;
  
  for (int i = 0; i < effective_N_reg; ++i) {
    double yhat = b0 + b1 * x_clean[i];
    sse += std::pow(y_clean[i] - yhat, 2);
  }
  
  double mse = sse / (effective_N_reg - 2);
  
  // --- 6. Calculate Zeta ---
  // Original formula: varpar = mse * (effective_N + 2) / effective_N;
  // This seems unusual. Often mse is used directly or var(residuals). Let's stick to it.
  double varpar = mse * (static_cast<double>(effective_N_reg) + 2.0) / static_cast<double>(effective_N_reg);
  
  double zeta_denominator = 0.0;
  if (swap_xy) { // B was X, A was Y
    zeta_denominator = (var_MP_B * std::pow(b1, 2) + var_MP_A);
  } else { // A was X, B was Y
    zeta_denominator = (var_MP_A * std::pow(b1, 2) + var_MP_B);
  }
  
  if (std::abs(zeta_denominator) < tolerance) {
    Rcpp::warning("Denominator for zeta calculation is zero or near zero.");
    return NA_REAL;
  }
  
  double zeta = varpar / zeta_denominator;
  
  return zeta;
}

//' @title Resample Clustered Data
//' @name resample_samples
//'
//' @param data A \code{list} or a \code{data.table}. Must contain \code{SampleID},
//'        \code{ReplicateID}, \code{MP_A} and \code{MP_B}. The ID variables \code{SampleID}
//'        and \code{ReplicateID} must be of character type for the function to operate correctly.
//'
//' @details
//' This function is a very efficient method to resample clinical sample data on sample-level.
//' It is convenient to combine this function with \code{fasteqa} functions such as
//' \itemize{
//'   \item \code{global_precision_estimates()} to estimate bootstrap imprecision confidence intervals
//'   \item \code{estimate_zeta()} to estimate bootstrap zeta confidence intervals
//' }
//'
//' @return
//' A \code{list} containing the resampled data.
//'
//' @examples
//' print(1)

// [[Rcpp::export]]
List resample_samples(List data) {
  CharacterVector samples = data["SampleID"];
  CharacterVector replicates = data["ReplicateID"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  CharacterVector unique_samples = unique_preserve_order(samples);
  int n = unique_samples.size();
  int N = samples.size();
  
  CharacterVector resampled_samples = sample(unique_samples, n, true);
  
  // Create a map for faster lookup
  std::unordered_map<String, std::vector<int>> sample_indices;
  for (int i = 0; i < N; ++i) {
    sample_indices[samples[i]].push_back(i);
  }
  
  // Calculate output size and prepare output vectors
  int output_size = 0;
  for (int i = 0; i < n; ++i) {
    output_size += sample_indices[resampled_samples[i]].size();
  }
  
  CharacterVector new_samples(output_size);
  CharacterVector new_replicates(output_size);
  NumericVector new_MP_A(output_size);
  NumericVector new_MP_B(output_size);
  
  // Fill output vectors
  int counter = 0;
  for (int i = 0; i < n; ++i) {
    String current_sample = resampled_samples[i];
    const std::vector<int>& indices = sample_indices[current_sample];
    for (int idx : indices) {
      new_samples[counter] = unique_samples[i];
      new_replicates[counter] = replicates[idx];
      new_MP_A[counter] = MP_A[idx];
      new_MP_B[counter] = MP_B[idx];
      ++counter;
    }
  }
  
  return List::create(
    Named("SampleID") = new_samples,
    Named("ReplicateID") = new_replicates,
    Named("MP_A") = new_MP_A,
    Named("MP_B") = new_MP_B
  );
}

//' @title Resample Cluster Statistics Based on External Quality Assessment (EQA) Clinical Sample Data
//' @name resample_fun_of_samples
//'
//' @param data A \code{list} or \code{data.table}. The mean of replicates (MOR)
//'        clinical sample data. Must contain:
//'        \itemize{
//'           \item \code{SampleID}: A \code{character} vector of sample identifiers.
//'           \item \code{MP_A}: A \code{numeric} vector of measurements from method A (response).
//'           \item \code{MP_B}: A \code{numeric} vector of measurements from method B (predictor).
//'        }
//'
//' @description
//' This function resamples resamples sample statistics based on EQA clinical sample data
//' with replacement. This function can for example be used to create bootstrap-based estimations
//' of predictions and conclusion strengths.
//' 
//' @details
//' This function is an optimized method to resample clinical sample statistics
//' based on clinical sample data.
//' 
//' Note, if \code{ReplicateID} is part of \code{data}, this method may yield
//' unexpected results. For methods that require resampling of all replicated
//' measurements, please use \code{resample_samples()} instead.
//' 
//' Note also that the \code{SampleID} of the returned object does not correspond
//' with the actual values in the input \code{SampleID}. Thus, this function
//' should not be used to do inference on individual clinical samples!
//'
//' @return
//' A \code{list} that contains the resampled \code{data}.
//'
//' @examples
//' library(fasteqa)
//' mor_test_data <- fun_of_replicates(test_data)
//' mor_test_data$MP_A <- round(mor_test_data$MP_A, 2)
//' mor_test_data$MP_B <- round(mor_test_data$MP_B, 2)
//' resampled_data <- resample_fun_of_samples(mor_test_data)
//' print(as.data.frame(resampled_data))
 
 // [[Rcpp::export]]
 List resample_fun_of_samples(List data) {
   
   // Error handling
   if (data.length() == 0) {
     stop("Input data is empty");
   }
   if (!data.containsElementNamed("SampleID") ||
       !data.containsElementNamed("MP_A") ||
       !data.containsElementNamed("MP_B")) {
       stop("Missing required columns in input data");
   }
   
   // Initialization
   CharacterVector samples = data["SampleID"];
   NumericVector MP_A = data["MP_A"];
   NumericVector MP_B = data["MP_B"];
   int n = samples.size();
   
   // Resample samples
   CharacterVector resampled_samples = sample(samples, n, true);
   
   // Create a map for faster lookup
   std::unordered_map<String, std::vector<int>> sample_indices;
   for (int i = 0; i < n; ++i) {
     sample_indices[samples[i]].push_back(i);
   }
   
   // Calculate output size and prepare output vectors
   int output_size = 0;
   for (int i = 0; i < n; ++i) {
     output_size += sample_indices[resampled_samples[i]].size();
   }
   
   // Set up new empty vectors
   CharacterVector new_samples(output_size);
   NumericVector new_MP_A(output_size);
   NumericVector new_MP_B(output_size);
   
   // Fill output vectors
   int counter = 0;
   for (int i = 0; i < n; ++i) {
     String current_sample = resampled_samples[i];
     const std::vector<int>& indices = sample_indices[current_sample];
     for (int idx : indices) {
       new_samples[counter] = samples[i];
       new_MP_A[counter] = MP_A[idx];
       new_MP_B[counter] = MP_B[idx];
       ++counter;
     }
   }
   
   // Structure of output
   List resampled_cs_data = List::create(Named("SampleID") = new_samples,
                                         Named("MP_A") = new_MP_A,
                                         Named("MP_B") = new_MP_B);
   
   
   return resampled_cs_data;
 }

//' @title Resample imprecision estimates based on clustered EQA clinical sample data 
//' @name resample_global_precision_estimates
//' 
//' @param data A \code{list} or \code{data.table}. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//'
//' @description
//' Resample \code{Var_B}, \code{Var_A}, \code{CV_B}, \code{CV_A} and \code{lambda} one time.
//'
//' @details
//' Uses \code{resample_samples()} to resample \code{data} on sample-level.
//' Then \code{global_imprecision_estimates()} is used to calculate the estimated values
//' of \code{Var_B}, \code{Var_A}, \code{CV_B}, \code{CV_A} and \code{lambda},
//' based on the resampled \code{data}. 
//'
//' @return
//' A \code{list} contaning the resampled imprecision estimates:
//' \itemize{
//'   \item \code{Var_A:} Pooled variance of all sample-variances based on \code{MP_A}
//'   \item \code{Var_B:} Pooled variance of all sample-variances based on \code{MP_B}
//'   \item \code{CV_A:} CV estimate based on Var_A and the grand mean of all measurements from \code{MP_A}
//'   \item \code{CV_B:} CV estimate based on Var_B and the grand mean of all measurements from \code{MP_B}
//'   \item \code{lambda:} Ratio of pooled variances \code{Var_A} and \code{Var_B}
//' }
//'
//' @examples 
//' library(fasteqa)
//' 
//' # Resample global imprecision estimates 3 times
//' rie <- replicate(n = 3,
//'                  expr = resample_global_precision_estimates(test_data),
//'                  simplify = FALSE)
//' rie <- lapply(X = rie,
//'               FUN = as.data.frame)
//' resampled_impr_estimates <- rbind(rie[[1]],
//'                                   rie[[2]],
//'                                   rie[[3]])
//' print(resampled_impr_estimates)
//' 

// [[Rcpp::export]]
List resample_global_precision_estimates(List data) {
  List resampled_samples = resample_samples(data);
  List output = global_precision_estimates_local(resampled_samples);
  return output;
} 

//' @title Resample \eqn{\hat{\zeta}} values based on clustered EQA clinical sample data 
//' @name resample_zeta_estimates
//' 
//' @param data A \code{list} or \code{data.table}. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//'
//' @description
//' Resample \code{zeta} one time.
//'
//' @details
//' Uses \code{resample_samples()} to resample \code{data} on sample-level.
//' Then \code{estimate_zeta()} is used to calculate \eqn{\hat{\zeta}},
//' based on the resampled \code{data}. 
//'
//' @return
//' A \code{list} contaning the resampled \eqn{\hat{\zeta}}.
//'
//' @examples
//' # Resample global imprecision estimates 1000 times
//' rze <- replicate(n = 1000,
//'                  expr = resample_zeta_estimates(test_data)$zeta,
//'                  simplify = TRUE)
//' # Get summary statistics of resampled zeta estimates
//' summary(rze)
//' 

// [[Rcpp::export]]
List resample_zeta_estimates(List data) {
 List resampled_samples = resample_samples(data);
 List output = List::create(Named("zeta") = estimate_zeta_local(resampled_samples));
 return output;
}

//' @title Leave-one-out on clustered EQA clinical sample data
//' @name leave_one_out
//' 
//' @param data A \code{list} or \code{data.table}. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//' @param loo_id An \code{integer}. The \code{SampleID} that is left out.
//'        Default value is \code{1}.
//'
//' @description Needed to calculate jack knife estimates of a parameter, that is required when using \code{BCa_bootstrap_ci()}. Alternatively one could use the Jack knife estimates to calculate standard error or bias of the estimator of relevance
//'
//' @details loo_ids can not be vectorized directly in R. Use \code{sapply()} or \code{replicate()} to leave out sample IDs one by one 
//'
//' @return A \code{list} containing the original data, but without the sample id corresponding to the given \code{loo_id} 
//'
//' @examples
//' # Get data without the fifth sample
//' print(as.data.frame(leave_one_out(test_data, 5L)))
//' 


// [[Rcpp::export]]
List leave_one_out(List data, int loo_id = 1){
 
 // Definitions
 CharacterVector samples = data["SampleID"];
 CharacterVector replicates = data["ReplicateID"];
 NumericVector MP_A = data["MP_A"];
 NumericVector MP_B = data["MP_B"];
 CharacterVector unique_samples = unique(samples);
 int n = unique_samples.size();
 int N = samples.size();
 int output_size = 0;
 int counter = 0;
 
 // How large should the output be?
 for(int i = 0; i < n; ++i){
   CharacterVector to_match(1);
   to_match[0] = unique_samples[i];
   LogicalVector matching_resamples = in(samples, to_match);
   int index_vector_length = sum(matching_resamples);
   if((i + 1) != loo_id){
     output_size += index_vector_length;  
   }
 }
 IntegerVector new_indicies(output_size);
 CharacterVector new_samples(output_size);
 CharacterVector new_replicates(output_size);
 NumericVector new_MP_A(output_size);
 NumericVector new_MP_B(output_size);
 
 
 for(int i = 0; i < n; ++i){
   // Empty vector
   CharacterVector to_match(1);
   // Fill empty vector with unique sample id number i
   to_match[0] = unique_samples[i];
   // Unique sample id in vector of all replicates
   LogicalVector matching_resamples = in(samples, to_match);
   // Count the number of replicates for sample id i
   int index_vector_length = sum(matching_resamples);
   // Create an empty index vector with same length as number of replicates for sample id i
   NumericVector index_vector(index_vector_length);
   int relevant_ind = 0;
   for(int j = 0; j < N; ++j){
     // if jth element is 1 and we are not at the sample ID that is left out, we record the index number
     if(matching_resamples[j] == 1 and (i + 1) != loo_id){
       index_vector[relevant_ind] = j;
       ++relevant_ind;
     }
   }
   for(int k = 0; k < index_vector_length; ++k){
     if((i + 1) != loo_id){
       new_indicies[counter] = index_vector[k];
       new_samples[counter] = unique_samples[i];
       ++counter;  
     }
   }
 }
 // Filling in new stuff
 for(int j = 0; j < output_size; ++j){
   new_replicates[j] = replicates[new_indicies[j]];
   new_MP_A[j] = MP_A[new_indicies[j]];
   new_MP_B[j] = MP_B[new_indicies[j]];
 }
 List output = List::create(Named("SampleID") = new_samples, Named("ReplicateID") = new_replicates, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A);
 return output;
}


//' @title Leave-One-Out (LOO) Imprecision Estimates Based on Clustered EQA Clinical Sample Data 
//' @name loo_global_precision_estimates
//' 
//' @param data A \code{list} or \code{data.table}. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//' @param loo_id A \code{integer}. The ID of the \code{SampleID} to be
//'        excluded before calculating the imprecision estimates.
//'
//' @description
//' Calculate \code{Var_B}, \code{Var_A}, \code{CV_B}, \code{CV_A} and \code{lambda} without
//' the \code{loo_id}-th SampleID.
//'
//' @details
//' Uses \code{leave_one_out()} to remove one unique value of \code{SampleID} from \code{data}
//' Then \code{global_imprecision_estimates()} is used to calculate the estimated values
//' of \code{Var_B}, \code{Var_A}, \code{CV_B}, \code{CV_A} and \code{lambda},
//' based on the \code{data} without the data from the \code{loo_id}-th \code{SampleID} value. 
//'
//' @return
//' A \code{list} contaning the LOO imprecision estimates:
//' \itemize{
//'   \item \code{Var_A:} Pooled variance of all sample-variances based on \code{MP_A}
//'   \item \code{Var_B:} Pooled variance of all sample-variances based on \code{MP_B}
//'   \item \code{CV_A:} CV estimate based on Var_A and the grand mean of all measurements from \code{MP_A}
//'   \item \code{CV_B:} CV estimate based on Var_B and the grand mean of all measurements from \code{MP_B}
//'   \item \code{lambda:} Ratio of pooled variances \code{Var_A} and \code{Var_B}
//' }
//'
//' @examples
//' library(fasteqa)
//' # Estimate imprecision estimates without the second clinical sample
//' gpe <- loo_global_precision_estimates(test_data, 2L)
//' 
//' # Output
//' print(as.data.frame(gpe))
//'

// [[Rcpp::export]]
List loo_global_precision_estimates(List data, int loo_id = 1) {
  List loo_data = leave_one_out(data, loo_id);
  List output = global_precision_estimates_local(loo_data);
  return output;
  
}


//' @title Calculate the \eqn{\hat{\zeta}} Value When One Observation is Removed
//' @name loo_zeta_estimates
//' 
//' @param data A \code{list} or \code{data.table}. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//'
//' @description
//' Resample \code{zeta} one time.
//'
//' @details
//' One unique value of \code{SampleID} is removed from \code{data}.
//' Then \code{estimate_zeta()} is used to calculate \eqn{\hat{\zeta}},
//' based on \code{data} without the data from the \code{loo_id}-th \code{SampleID} value.
//'
//' @return
//' A \code{list} contaning the resampled \eqn{\hat{\zeta}}.
//'
//' @examples \dontrun{
//'   # Estimate zeta without the second clinical sample
//'   print(loo_zeta_estimatees(test_data, 2L)$zeta)
//' }

// [[Rcpp::export]]
List loo_zeta_estimates(List data, int loo_id = 1) {
 List loo_data = leave_one_out(data, loo_id);
 List output = List::create(Named("zeta") = estimate_zeta_local(loo_data));
 return output;
}

// Helper function (*****) 
NumericVector BCa_bootstrap_ci(NumericVector bootstrapped_parameter_estimates,
                               NumericVector jackknife_parameter_estimates,
                               double original_parameter_estimate = 0,
                               double level = 0.95, int silence = 1){
 
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
   jpe.push_back(raw_jpe[i]);
 }
 
 if(silence < 1){
   Rcout << B << " bootstrap estimates are non-missing out of " << n << " recorded values" << "\n";
   Rcout << J << " jackknife estimates are non-missing out of " << m << " recorded values" << "\n";
   Rcout << "-------------------------------------------" << "\n";
 }
 
 // Converting integer types to double types when needed in division!
 double B_double = static_cast<double>(B);
 
 NumericVector jackknife_bias_cubed(J);
 NumericVector jackknife_bias_squared(J);
 
 double alpha = 1 - level;
 double jackknife_mean = mean(jpe);
 double z0_kernel = 0;
 
 // Calculation of bias-correction z0
 for(int i = 0; i < n; ++i){
   if(sorted_bpe[i] < original_parameter_estimate){
     z0_kernel += 1.0 / B_double;
   }
 }
 
 double z0 = R::qnorm5(z0_kernel, 0, 1, 1, 0);
 // Calculation of the skewness-correction constant a_jackknife, via jack knife (leave-one-out cross validation)
 for(int j = 0; j < J; ++j){
   jackknife_bias_cubed[j] = pow(jackknife_mean - jackknife_parameter_estimates[j], 3);
   jackknife_bias_squared[j] = pow(jackknife_mean - jackknife_parameter_estimates[j], 2);
 }
 double num = mean(jackknife_bias_cubed);
 double den = 6.0 * pow(sum(jackknife_bias_squared), 1.5);
 double a_jackknife = num / den;
 if(silence < 1){
   Rcout << "The bias-correction coefficient is estimated to be: " << z0 << "\n";
   Rcout << "The skewness-correction coefficient is estimated to be: " << a_jackknife << "\n";
   Rcout << "-------------------------------------------" << "\n";
 }
 double z_alpha_halves = R::qnorm5(1.0 - alpha * 0.5, 0, 1, 1, 0);
 double new_lower_quant_kernel = z0 + ((z0 - z_alpha_halves) / (1.0 - a_jackknife * (z0 - z_alpha_halves)));
 double new_upper_quant_kernel = z0 + ((z0 + z_alpha_halves) / (1.0 - a_jackknife * (z0 + z_alpha_halves)));
 double new_lower_quant = R::pnorm5(new_lower_quant_kernel, 0, 1, 1, 0);
 double new_upper_quant = R::pnorm5(new_upper_quant_kernel, 0, 1, 1, 0);
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

//' @title Simple Bootstrap Confidence Interval Estimation
//' @name bootstrap_ci
//' 
//' @param bootstrapped_parameter_estimates A \code{numeric} vector containing 
//'        resampled parameter estimates. Missing values (\code{NA}) are automatically 
//'        removed without a warning. Must contain at least \eqn{2} valid observations.
//' @param jackknife_parameter_estimates A \code{numeric} vector containing Jackknife
//'        estimates of the target parameter.
//'        Missing values (\code{NA}) will be removed automatically.
//' @param original_parameter_estimate A \code{double} representing the parameter 
//'        estimate from the original sample. Must be a single non-missing value.
//' @param type An `integer` specifying the confidence interval method:
//' \itemize{
//'   \item \code{1} - Standard normal bootstrap confidence interval
//'   \item \code{2} - Basic bootstrap confidence interval
//'   \item \code{3} - Percentile bootstrap confidence interval
//'   \item \code{4} - Bias and skewness corrected bootstrap confidence interval
//' }
//' @param level A \code{double} between 0.50 and 0.99 indicating the confidence level. 
//'        Default: 0.95 (95\% interval)
//' @param silence An \code{integer} controlling verbosity:
//' \itemize{
//'   \item \code{-1} - Detailed debug output
//'   \item \code{0} - Basic progress reporting
//'   \item \code{1} - No output (default)
//' }
//'
//' @description
//' Calculates nonparametric bootstrap confidence intervals using four common methods.
//' Suitable for situations where only bootstrap replicates and the original estimate 
//' are available.
//'
//' @details
//' Methodological Notes
//' \itemize{
//'   \item Standard Normal (Type 1): Assumes normal distribution of bootstrap estimates
//'   \item Basic (Type 2): Computes interval using bootstrap distribution's quantiles
//'   \item Percentile (Type 3): Directly uses empirical quantiles of bootstrap estimates
//'   \item BCa (Type 4): Uses bias- and skewness corrected bootstrap quantiles
//' }
//' 
//' The function automatically removes \code{NA} values from bootstrap estimates without a warning.
//' For reliable results, use \eqn{\geq 1000} bootstrap replicates.
//' BCa intervals require jackknife estimates. 
//'
//' @return A \code{numeric} vector \code{x} of length 2 containing:
//' \itemize{
//'   \item \code{x[1]} - Lower confidence bound
//'   \item \code{x[2]} - Upper confidence bound
//' }
//'
//' @examples
//' # Simulated dataset
//' true_mean <- 5.2
//' bootstrap_means <- rnorm(1000, mean = 5.0, sd = 0.5)
//' jackknife_means <- rnorm(25, mean = 5.1, sd = 0.75)
//' 
//' # Calculate 90% percentile interval
//' bootstrap_ci(
//'   bootstrapped_parameter_estimates = bootstrap_means,
//'   jackknife_parameter_estimates = bootstrap_means,
//'   original_parameter_estimate = true_mean,
//'   type = 3,
//'   level = 0.90
//' )
//' 
//' # Calculate 95% BCa interval
//' bootstrap_ci(
//'   bootstrapped_parameter_estimates = bootstrap_means,
//'   jackknife_parameter_estimates = jackknife_means,
//'   original_parameter_estimate = true_mean,
//'   type = 4,
//'   level = 0.95
//' )
//'
//' @references
//' Bootstrap methods and their application (Davison & Hinkley, 1997)

// [[Rcpp::export]]
NumericVector bootstrap_ci(NumericVector bootstrapped_parameter_estimates,
                           NumericVector jackknife_parameter_estimates,
                           double original_parameter_estimate = 0,
                           int type = 3,
                           double level = 0.95,
                           int silence = 1){
  
  if(type == 4){
    return BCa_bootstrap_ci(bootstrapped_parameter_estimates,
                            jackknife_parameter_estimates,
                            original_parameter_estimate,
                            level,
                            silence);
  }
  
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
  
  double se_parameter = sd(sorted_bpe);
  double bias_parameter = original_parameter_estimate - mean(sorted_bpe);
  
  if(silence < 1){
    Rcout << "The " << lower + 1 << "th smallest bootstrap estimate is used as lower quantile" << "\n";  
    Rcout << "The " << upper + 1 << "th smallest bootstrap estimate is used as upper quantile" << "\n"; 
    Rcout << "-------------------------------------------" << "\n";
  }
  
  NumericVector bci(2);
  
  if(type == 1){
    bci[0] = original_parameter_estimate + bias_parameter - R::qnorm(lower_quant, 0, 1, 0, 0) * se_parameter;
    bci[1] = original_parameter_estimate + bias_parameter + R::qnorm(lower_quant, 0, 1, 0, 0) * se_parameter;
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





