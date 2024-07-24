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
    // If ISNAN results in 1 (i.e., NA-value) let the kth measurement be zero
    // Otherwise, the kth measurement will be the kth measurement of MS_A / MS_B
    for(int k = 0; k < number_of_matches; ++k){
      bool is_na_MS_A = ISNAN(MS_A[indices[k]]);
      bool is_na_MS_B = ISNAN(MS_B[indices[k]]);
      if(!is_na_MS_A){
        ith_sample_measurements_A[k] = MS_A[indices[k]];
      }
      else if(is_na_MS_A){
        ith_sample_measurements_A[k] = 0;
      }
      if(!is_na_MS_B){
        ith_sample_measurements_B[k] = MS_B[indices[k]];  
      }
      else if(is_na_MS_B){
        ith_sample_measurements_B[k] = 0;
      }
    }
    
    // Create NA-search vectors with zeros and ones matching with ith sample of MS_A and MS_B
    IntegerVector NA_search_A(number_of_matches);
    IntegerVector NA_search_B(number_of_matches);
    
    // recall: ith_sample_measurements_*[k] = 0 signify that the value is a NA-value
    // NA_search_* will return 0 if ith_sample_measurements_*[k] = 0, and 1 otherwise 
    for(int k = 0; k < number_of_matches; ++k){
      if(ith_sample_measurements_A[k] == 0){
        NA_search_A[k] = 0;
      }
      else if(ith_sample_measurements_A[k] > 0){
        NA_search_A[k] = 1;
      }
      else if(ith_sample_measurements_A[k] < 0){
        NA_search_A[k] = 0;
      }
      if(ith_sample_measurements_B[k] == 0){
        NA_search_B[k] = 0;
      }
      else if(ith_sample_measurements_B[k] > 0){
        NA_search_B[k] = 1;
      }
      else if(ith_sample_measurements_B[k] < 0){
        NA_search_B[k] = 0;
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
    bool is_na_var_A_j = ISNAN(ith_var_MS_A[j]);
    bool is_na_var_B_j = ISNAN(ith_var_MS_B[j]);
    if(!is_na_var_A_j){
      var_MS_A += ith_var_MS_A[j];
    }
    else if(is_na_var_A_j){
      effective_n_A = effective_n_A - 1;
    }
    if(!is_na_var_B_j){
      var_MS_B += ith_var_MS_B[j];
    }
    else if(is_na_var_B_j){
      effective_n_B = effective_n_B - 1;
    }
  }
  
  if(effective_n_A >= 1){
    var_MS_A = var_MS_A / effective_n_A; 
  }
  if(effective_n_A < 1){
    var_MS_A = NA_REAL;
  }
  if(effective_n_B >= 1){
    var_MS_B = var_MS_B / effective_n_B;
  }
  if(effective_n_B < 1){
    var_MS_B = NA_REAL;
  }
  float lambda = 0;
  bool is_na_pooled_var_A = ISNAN(var_MS_A);
  bool is_na_pooled_var_B = ISNAN(var_MS_B);
  if(is_na_pooled_var_A | is_na_pooled_var_B){
    List out = List::create(Named("zeta") = NA_REAL);
    return out;
  }
  else if((!is_na_pooled_var_A) & (!is_na_pooled_var_B)){
    lambda = var_MS_A / var_MS_B;  
  }
  float mean_MS_A = 0;
  float mean_MS_B = 0;
  
  for(int i = 0; i < N; ++i){
    bool is_na_mean_A = ISNAN(MS_A[i]);
    bool is_na_mean_B = ISNAN(MS_B[i]);
    if(!is_na_mean_A){
      mean_MS_A += MS_A[i];
    }
    else if(is_na_mean_A){
      effective_N_A = effective_N_A - 1;
    }
    if(!is_na_mean_B){
      mean_MS_B += MS_B[i];
    }
    else if(is_na_mean_B){
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
  float mx = mean_MS_B;
  float my = mean_MS_A;
  
  float msxx = 0;
  float msyy = 0;
  float msxy = 0;
  
  for(int i = 0; i < N; ++i){
    bool na_check_x = ISNAN(x[i]);
    bool na_check_y = ISNAN(y[i]);
    if((!na_check_x) & (!na_check_y)){
      msxx = msxx + pow(x[i] - mx, 2);
      msyy = msyy + pow(y[i] - my, 2);
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
  }
  
  msxx = msxx / (effective_N_B - 1.0);
  msyy = msyy / (effective_N_A - 1.0);
  msxy = msxy / (effective_N_A - 1.0);
  
  float sub_expression_1 = msyy - lambda * msxx;
  float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
  float sub_expression_3 = 2 * msxy;
  float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
  float varb1 = (pow(b1, 2) / ((effective_N_A) * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
  float hvar = (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda);
  float varpar = varb1 * msxx + varb1 * hvar + (1 + 1.0 / effective_N_A) * (pow(b1, 2) + lambda) * hvar;
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
