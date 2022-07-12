#include <Rcpp.h>
using namespace Rcpp;

//' Calculate imprecision point estimates of measurements in a given MS comparison 
//' 
//' @title Calculate imprecision point estimates of measurements in a given MS comparison
//' @name global_precision_estimates
//' 
//' @param data \code{list} or \code{data table} - Data with elements/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
//' @param silence \code{integer} - Should debugging messages be printed. Default is no
//'
//' @description Structure-requirement of \code{data}:
//' \itemize{
//'   \item{\code{SampleID: }}{Must be a character vector}
//'   \item{\code{ReplicateID: }}{Must be a character vector}
//'   \item{\code{MP_A: }}{Must be a numeric vector}
//'   \item{\code{MP_B: }}{Must be a numeric vector}
//' }
//'
//' @details Calculates various global imprecision estimates. To get CVs in percent you need only to multiply the raw CV estimates with 100. Here is a rough explaination of the output list:
//' \itemize{
//'   \item{\code{Var_A: }}{Pooled variance of all sample-variances based on MS_A}
//'   \item{\code{Var_B: }}{Pooled variance of all sample-variances based on MS_B}
//'   \item{\code{CV_A: }}{Global CV estimate based on Var_A and the grand mean of all measurements from MS_A}
//'   \item{\code{CV_B: }}{Global CV estimate based on Var_B and the grand mean of all measurements from MS_B}
//'   \item{\code{lambda: }}{Ratio of pooled variances Var_A and Var_B}
//' }
//' 
//'
//' @return \code{list} - with point imprecision estimates \code{Var_A}, \code{Var_B}, \code{CV_A}, \code{CV_B} and \code{lambda}
//' @examples \dontrun{
//'   library(fasteqa)
//'   data <- simulate_eqa_data(list(n = 25, R = 3, cvx = 0.02, cvy = 0.3))
//'   data$SampleID <- as.character(data$SampleID)
//'   data$ReplicateID <- as.character(data$ReplicateID)
//'   global_prcision_estimates(data = data)
//' }


// [[Rcpp::export]]
List global_precision_estimates(List data, int silence = 1) {
  CharacterVector SampleID = data["SampleID"];
  CharacterVector ReplicateID = data["ReplicateID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  CharacterVector summary_SampleID = unique(SampleID);
  int n = summary_SampleID.size();
  int N = SampleID.size();
  NumericVector ith_var_MS_A(n);
  NumericVector ith_var_MS_B(n);
  int counter = 0;
  int replicate_number_requirement = 2;
  
  for(int i = 0; i < n; ++i){
    // Extract necessary information on the ith sample (indices and number of replicated meas)
    CharacterVector ith_sample(1);
    ith_sample[0] = summary_SampleID[i];
    LogicalVector sample_match = in(SampleID, ith_sample);
    int index_vector_length = sum(sample_match);
    NumericVector index_vector(index_vector_length);
    int relevant_ind = 0;
    for(int j = 0; j < N; ++j){
      if(sample_match[j] == 1){
        index_vector[relevant_ind] = j;
        ++relevant_ind;
      }
    }
    
    // Create empty vectors to be filled with measurements of the ith sample
    NumericVector ith_sample_measurements_A(index_vector_length);
    NumericVector ith_sample_measurements_B(index_vector_length);
    
    // Checks for NA-values
    // If R_IsNA results in 1 (i.e., NA-value) let the kth measurement be zero
    // Otherwise, the kth measurement will be the kth measurement of MS_A / MS_B
    for(int k = 0; k < index_vector_length; ++k){
      int na_check_A = R_IsNA(MS_A[index_vector[k]]);
      if(na_check_A == 0){
        ith_sample_measurements_A[k] = MS_A[index_vector[k]];
      }
      int na_check_B = R_IsNA(MS_B[index_vector[k]]);
      if(na_check_B == 0){
        ith_sample_measurements_B[k] = MS_B[index_vector[k]];  
      }
    }
    
    // Create NA-search vectors with zeros and ones matching with ith sample of MS_A and MS_B
    IntegerVector NA_search_A(index_vector_length);
    IntegerVector NA_search_B(index_vector_length);
    
    // recall: ith_sample_measurements_*[k] = 0 signify that the value is a NA-value
    // NA_search_* will return 0 if ith_sample_measurements_*[k] = 0, and 1 otherwise 
    for(int k = 0; k < index_vector_length; ++k){
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
    for(int k = 0; k < index_vector_length; ++k){
      if(NA_search_A[k] == 0){
        ith_sample_measurements_A.erase(0); 
      }
    }
    
    // Step-wise check for NA values at vector start and delete if it is a NA-value for MS_A
    for(int k = 0; k < index_vector_length; ++k){
      if(NA_search_B[k] == 0){
        ith_sample_measurements_B.erase(0);
      }
    }
    
    // Checks if the number of replicates meets the conditions of the relevant summary function
    // If met, the summary function is applied. Otherwise, the result will be a NA-value
    if(ith_sample_measurements_A.size() < replicate_number_requirement){
      ith_var_MS_A[counter] = NA_REAL;  
    }
    if(ith_sample_measurements_A.size() >= replicate_number_requirement){
      ith_var_MS_A[counter] = var(ith_sample_measurements_A);
    }
    if(ith_sample_measurements_B.size() < replicate_number_requirement){
      ith_var_MS_B[counter] = NA_REAL;
    }
    if(ith_sample_measurements_B.size() >= replicate_number_requirement){
      ith_var_MS_B[counter] = var(ith_sample_measurements_B);  
    }
    ++counter;
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
  
  float cv_MS_A = 0;
  float cv_MS_B = 0;
  
  if(effective_N_A >= 1 and effective_n_A >= 1){
    cv_MS_A = sqrt(var_MS_A) / mean_MS_A;
  }
  if(effective_N_A < 1 or effective_n_A < 1){
    cv_MS_A = NA_REAL;
  }
  if(effective_N_B >= 1 and effective_n_B >= 1){
    cv_MS_B = sqrt(var_MS_B) / mean_MS_B;
  }
  if(effective_N_B < 1 or effective_n_B < 1){
    cv_MS_B = NA_REAL;
  }
  
  List out = List::create(Named("Var_A") = var_MS_A, Named("Var_B") = var_MS_B, Named("CV_A") = cv_MS_A, Named("CV_B") = cv_MS_B, Named("lambda") = lambda);
  return out; 
}
  


