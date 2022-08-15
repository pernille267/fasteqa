#include <Rcpp.h>
using namespace Rcpp;

//' Estimate differences in non-selectivity with zeta using Deming regression
//' 
//' @title Estimate differences in non-selectivity with zeta using Deming regression
//' @name estimate_zeta_deming
//' @param data \code{list} or \code{data table} - Data with elements/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
//' @param silence \code{integer} - How much progress reports should be returned. Note that returning progress reports will slow down the performance drastically. There are three valid inputs:
//' \itemize{
//'   \item{\code{1: }}{All progress reports are silenced and is the default}
//'   \item{\code{0: }}{Some progress reports are delivered, but debugging reports are suppressed}
//'   \item{\code{-1: }}{All prorgress reports are delivered}
//' }
//' 
//' @description Estimate the degree of differences in non-selectivity with zeta. Zeta is is the ratio of the pooled average prediction error variance and the sum of analytical variances.  
//' 
//' @details Differences in non-selectivity between measurement systems may cause problems in e.g., evaluation of commutability. A large value of zeta indicates that we have have large differences in non-selectivity between compared measurement systems. An upper limit of acceptable zeta may be determined based on the allowable increase in prediction interval width and analyte of relevance
//' 
//' @return A list with the point estimate of zeta. The zeta value is a float value, meaning that the precision is 1e-6 (six decimals precision).
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   data <- simulate_data_eqa(list(n = 25, R = 3, qran = 0.20, qpos = 1, mmax = 8))
//'   estimate_zeta_deming(data)
//' }


// [[Rcpp::export]]
List estimate_zeta_deming(List data, int silence = 1) {
  
  // Extract components from data
  CharacterVector SampleID = data["SampleID"];
  CharacterVector ReplicateID = data["ReplicateID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  CharacterVector summary_SampleID = unique(SampleID);
  int n = summary_SampleID.size();
  int N = SampleID.size();
  NumericVector ith_var_MS_A(n);
  NumericVector ith_var_MS_B(n);
  int c = 0;
  int replicate_number_requirement = 2;
  
  for(int i = 0; i < n; ++i){
    // Extract necessary information on the ith sample (indices and number of replicated meas)
    CharacterVector ith_sample(1);
    ith_sample[0] = summary_SampleID[i];
    LogicalVector matches = in(SampleID, ith_sample);
    int number_of_matches = sum(matches);
    NumericVector indices(number_of_matches);
    int k = 0;
    for(int j = 0; j < N; ++j){
      if(matches[j] == 1){
        indices[k] = j;
        ++k;
      }
    }
    
    // Create empty vectors to be filled with measurements of the ith sample
    NumericVector ith_sample_measurements_A(number_of_matches);
    NumericVector ith_sample_measurements_B(number_of_matches);
    
    // Checks for NA-values
    // If R_IsNA results in 1 (i.e., NA-value) let the kth measurement be zero
    // Otherwise, the kth measurement will be the kth measurement of MS_A / MS_B
    for(int k = 0; k < number_of_matches; ++k){
      int na_check_A = R_IsNA(MS_A[indices[k]]);
      if(na_check_A == 0){
        ith_sample_measurements_A[k] = MS_A[indices[k]];
      }
      int na_check_B = R_IsNA(MS_B[indices[k]]);
      if(na_check_B == 0){
        ith_sample_measurements_B[k] = MS_B[indices[k]];  
      }
    }
    
    // Create NA-search vectors with zeros and ones matching with ith sample of MS_A and MS_B
    IntegerVector NA_search_A(number_of_matches);
    IntegerVector NA_search_B(number_of_matches);
    
    // recall: ith_sample_measurements_*[k] = 0 signify that the value is a NA-value
    // NA_search_* will return 0 if ith_sample_measurements_*[k] = 0, and 1 otherwise 
    for(int k = 0; k < number_of_matches; ++k){
      if(ith_sample_measurements_A[k] <= 0){
        NA_search_A[k] = 0;
      }
      if(ith_sample_measurements_A[k] > 0){
        NA_search_A[k] = 1;
      }
      if(ith_sample_measurements_B[k] <= 0){
        NA_search_B[k] = 0;
      }
      if(ith_sample_measurements_B[k] > 0){
        NA_search_B[k] = 1;
      }
    }
    
    // Aligning vectors making it easier to delete NA-values
    ith_sample_measurements_A = ith_sample_measurements_A.sort();
    ith_sample_measurements_B = ith_sample_measurements_B.sort();
    NA_search_A = NA_search_A.sort();
    NA_search_B = NA_search_B.sort();
    
    // Step-wise check for NA values at vector start and delete if it is a NA-value for MS_A
    for(int k = 0; k < number_of_matches; ++k){
      if(NA_search_A[k] == 0){
        ith_sample_measurements_A.erase(0); 
      }
    }
    
    // Step-wise check for NA values at vector start and delete if it is a NA-value for MS_B
    for(int k = 0; k < number_of_matches; ++k){
      if(NA_search_B[k] == 0){
        ith_sample_measurements_B.erase(0);
      }
    }
    
    // Checks if the number of replicates meets the conditions of the relevant summary function
    // If met, the summary function is applied. Otherwise, the result will be a NA-value
    if(ith_sample_measurements_A.size() < replicate_number_requirement){
      ith_var_MS_A[c] = NA_REAL;  
    }
    if(ith_sample_measurements_A.size() >= replicate_number_requirement){
      ith_var_MS_A[c] = var(ith_sample_measurements_A);
    }
    if(ith_sample_measurements_B.size() < replicate_number_requirement){
      ith_var_MS_B[c] = NA_REAL;
    }
    if(ith_sample_measurements_B.size() >= replicate_number_requirement){
      ith_var_MS_B[c] = var(ith_sample_measurements_B);  
    }
    ++c;
  }
  
  
  int effective_N_A = N;
  int effective_n_A = n;
  int effective_N_B = N;
  int effective_n_B = n;
  float var_MS_A = 0;
  float var_MS_B = 0;
  
  for(int j = 0; j < n; ++j){
    int na_check_A = ISNAN(ith_var_MS_A[j]);
    int na_check_B = ISNAN(ith_var_MS_B[j]);
    if(na_check_A <= 0){
      var_MS_A += ith_var_MS_A[j];
    }
    if(na_check_A > 0){
      effective_n_A = effective_n_A - 1;
    }
    if(na_check_B <= 0){
      var_MS_B += ith_var_MS_B[j];
    }
    if(na_check_B > 0){
      effective_n_B = effective_n_B - 1;
    }
  }
  
  if(effective_n_A > 0){
    var_MS_A = var_MS_A / effective_n_A; 
  }
  if(effective_n_A <= 0){
    var_MS_A = NA_REAL;
  }
  if(effective_n_B > 0){
    var_MS_B = var_MS_B / effective_n_B;
  }
  if(effective_n_B <= 0){
    var_MS_B = NA_REAL;
  }
  float lambda = 0;
  int na_check_A = ISNAN(var_MS_A);
  int na_check_B = ISNAN(var_MS_B);
  if(na_check_A <= 0 and na_check_B <= 0){
    lambda = var_MS_A / var_MS_B;  
  }
  if(na_check_A > 0 or na_check_B > 0){
    lambda = NA_REAL;
  }
  
  float mean_MS_A = 0;
  float mean_MS_B = 0;
  
  for(int i = 0; i < N; ++i){
    int na_check_A = ISNAN(MS_A[i]);
    int na_check_B = ISNAN(MS_B[i]);
    if(na_check_A <= 0){
      mean_MS_A += MS_A[i];
    }
    if(na_check_A > 0){
      effective_N_A = effective_N_A - 1;
    }
    if(na_check_B <= 0){
      mean_MS_B += MS_B[i];
    }
    if(na_check_B > 0){
      effective_N_B = effective_N_B - 1;
    }
  }
  if(effective_N_A >= 1){
    mean_MS_A = mean_MS_A / effective_N_A;  
  }
  if(effective_N_A < 1){
    mean_MS_A = NA_REAL;
  }
  if(effective_N_B >= 1){
    mean_MS_B = mean_MS_B / effective_N_B;
  }
  if(effective_N_B < 1){
    mean_MS_B = NA_REAL;  
  }
  
  if(silence == 0 or silence == -1){
    Rcout << "----- Progress (1 / 3) Calculations of variances and means -------" << "\n";
    Rcout << "## Pooled variance of MS_A is : " << var_MS_A << "\n";
    Rcout << "## Pooled variance of MS_B is : " << var_MS_B << "\n";
    Rcout << "## Mean of MS_A is : " << mean_MS_A << "\n";
    Rcout << "## Mean of MS_B is : " << mean_MS_B << "\n";
    Rcout << "## lambda is : " << lambda << "\n";
    Rcout << "## N is : " << N << "\n";
    Rcout << "------------------------------------------------------------------" << "\n";  
  }
  
  NumericVector x = MS_B;
  NumericVector y = MS_A;
  float mx = mean(x);
  float my = mean(y);
  float msxx = var(x);
  float msyy = var(y);
  float msxy = 0;
  for(int i = 0; i < N; ++i){
    msxy = msxy + (x[i] - mx) * (y[i] - my);
  }
  msxy = msxy / (N - 1);
  float sub_expression_1 = msyy - lambda * msxx;
  float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
  float sub_expression_3 = 2 * msxy;
  float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
  float varb1 = (pow(b1, 2) / (N * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
  float hvar = (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda);
  float varpar = varb1 * msxx + varb1 * hvar + (1 + 1 / N) * (pow(b1, 2) + lambda) * hvar;
  float zeta = varpar / (var_MS_A + var_MS_B * pow(b1, 2));
  
  if(silence == 0 or silence == -1){
    Rcout << "----- Progress (2 / 3) Calculations of variances and means -------" << "\n";
    Rcout << "## var[y0 - y0^] : " << varpar << "\n";
    Rcout << "## b1 : " << b1 << "\n";
    Rcout << "## Pooled variance of MS_B is : " << var_MS_B << "\n";
    Rcout << "## Estimated variance of MS_B is : " << hvar << "\n";
    Rcout << "------------------------------------------------------------------" << "\n";  
  }
  
  List out = List::create(Named("zeta") = zeta);
  return out;
}
