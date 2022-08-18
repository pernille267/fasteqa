#include <Rcpp.h>
using namespace Rcpp;

//' Estimate differences in non-selectivity with zeta
//' 
//' @title Estimate differences in non-selectivity with zeta
//' @name estimate_zeta
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
//'   data <- simulate_data_eqa(list(n = 25, R = 3, cvx = 0.06, cvy = 0.04))
//'   estimate_zeta(data)
//' }
//'

// [[Rcpp::export]]
List estimate_zeta(List data, int silence = 1) {
  
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
    // If ISNAN results in TRUE (i.e., NA-value) let the kth measurement be zero
    // Otherwise, the kth measurement will be the kth measurement of the corresponding MS_A or MS_B
    for(int k = 0; k < number_of_matches; ++k){
      bool is_na_A = ISNAN(MS_A[indices[k]]);
      bool is_na_B = ISNAN(MS_B[indices[k]]);
      if(!is_na_A){
        ith_sample_measurements_A[k] = MS_A[indices[k]];
      }
      else{
        ith_sample_measurements_A[k] = 0;
      }
      if(!is_na_B){
        ith_sample_measurements_B[k] = MS_B[indices[k]];  
      }
      else{
        ith_sample_measurements_B[k] = 0;
      }
    }
    
    // Create NA-search vectors with zeros and ones matching with ith sample of MS_A and MS_B
    IntegerVector NA_search_A(number_of_matches);
    IntegerVector NA_search_B(number_of_matches);
    
    // Recall: ith_sample_measurements_*[k] = 0 signify that the value is a NA-value
    // NA_search_* will return 0 if ith_sample_measurements_*[k] = 0, and 1 otherwise 
    for(int k = 0; k < number_of_matches; ++k){
      if(ith_sample_measurements_A[k] == 0){
        NA_search_A[k] = 0;
      }
      if(ith_sample_measurements_A[k] > 0){
        NA_search_A[k] = 1;
      }
      if(ith_sample_measurements_B[k] == 0){
        NA_search_B[k] = 0;
      }
      if(ith_sample_measurements_B[k] > 0){
        NA_search_B[k] = 1;
      }
      if((ith_sample_measurements_A[k] < 0) | (ith_sample_measurements_B[k] < 0)){
        if((ith_sample_measurements_A[k] < 0) & (silence == -2)){
          Rcout << "SampleID " << i << " and replicate measurement " << k << " is registered as a negative values: "<< ith_sample_measurements_A[k] << "\n";  
        }
        if((ith_sample_measurements_B[k] < 0) & (silence == -2)){
          Rcout << "SampleID " << i << " and replicate measurement " << k << " is registered as a negative values: "<< ith_sample_measurements_B[k] << "\n";  
        }
        stop("Negative measurements are recorded");
      }
    }
    
    // Aligning vectors NA search vectors and measurement vectors making it easier to delete NA-values
    ith_sample_measurements_A = ith_sample_measurements_A.sort();
    ith_sample_measurements_B = ith_sample_measurements_B.sort();
    NA_search_A = NA_search_A.sort();
    NA_search_B = NA_search_B.sort();
    
    // Rcout << "ith sample measurement for A: " << ith_sample_measurements_A << "\n";
    // Rcout << "ith NA values: " << NA_search_A << "\n";
    
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
    bool is_na_var_A = ISNAN(ith_var_MS_A[j]);
    bool is_na_var_B = ISNAN(ith_var_MS_B[j]);
    if(!is_na_var_A){
      var_MS_A += ith_var_MS_A[j];
    }
    if(is_na_var_A){
      effective_n_A = effective_n_A - 1;
    }
    if(!is_na_var_B){
      var_MS_B += ith_var_MS_B[j];
    }
    if(is_na_var_A){
      effective_n_B = effective_n_B - 1;
    }
  }
  
  if(effective_n_A >= 1){
    var_MS_A = var_MS_A / effective_n_A; 
  }
  else if(effective_n_A < 1){
    if(silence == -2){
      Rcout << "effective sample size for MS_A was smaller than 1..." << "\n";
    }
    var_MS_A = NA_REAL;
  }
  if(effective_n_B >= 1){
    var_MS_B = var_MS_B / effective_n_B;
  }
  else if(effective_n_B < 1){
    if(silence == -2){
      Rcout << "effective sample size for MS_B was smaller than 1..." << "\n";
    }
    var_MS_B = NA_REAL;
  }
  float lambda = 0;
  bool is_na_pooled_var_MS_A = ISNAN(var_MS_A);
  bool is_na_pooled_var_MS_B = ISNAN(var_MS_B);
  if((!is_na_pooled_var_MS_A) & (!is_na_pooled_var_MS_B)){
    lambda = var_MS_A / var_MS_B;  
  }
  if(is_na_pooled_var_MS_A | is_na_pooled_var_MS_B){
    lambda = 1;
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
  
  // Calculate varpar and then zeta
  
  if((lambda < 1) & (!is_na_pooled_var_MS_A) & (!is_na_pooled_var_MS_B)){
    
    NumericVector x = MS_A;
    NumericVector y = MS_B;
    
    float mx = mean_MS_A;
    float my = mean_MS_B;
    
    float sxx = 0;
    float sxy = 0;
    float sse = 0;
    for(int i = 0; i < N; ++i){
      bool na_check_x = ISNAN(x[i]);
      bool na_check_y = ISNAN(y[i]);
      if((!na_check_x) & (!na_check_y)){
        sxx = sxx + pow(x[i] - mx, 2);
        sxy = sxy + (x[i] - mx) * (y[i] - my);
      }
    }
    if(silence == -1){
      Rcout << "----- Progress (2 / 3) Calculations of sum of squares ------------" << "\n";
      Rcout << "Sum of squares of x is : " << sxx << "\n";
      Rcout << "Sum of squares of x and y is : " << sxy << "\n";
      Rcout << "------------------------------------------------------------------" << "\n";  
    }
    
    float b1 = sxy / sxx;
    float b0 = my - b1 * mx;
    
    int effective_N = N;
    for(int i = 0; i < N; ++i){
      bool na_check_x = ISNAN(x[i]);
      bool na_check_y = ISNAN(y[i]);
      if((!na_check_x) & (!na_check_y)){
        float yhat = b0 + b1 * x[i];  
        sse = sse + pow(y[i] - yhat, 2);
      }
      else{
        effective_N = effective_N - 1;
        continue;
      }
    }
    float mse = sse / (effective_N - 2);
    float varpar = mse * (effective_N + 2);
    varpar = varpar / effective_N;
    
    if(silence == 0 or silence == -1){
      Rcout << "----- Progress (3 / 3) remaining calculations --------------------" << "\n";
      Rcout << "## Estimtated slope estimator : " << b1 << "\n";
      Rcout << "## Estimtated intercept estimator : " << b0 << "\n";
      Rcout << "## Estimated sum of squares error : " << sse << "\n";
      Rcout << "## Estimated mean sum of squares error : " << mse << "\n";
      Rcout << "## Effective N (after removing NA values) : " << effective_N << "\n";
      Rcout << "------------------------------------------------------------------" << "\n";  
    }
    
    
    float zeta = varpar / (var_MS_A * pow(b1, 2) + var_MS_B);
    
    List out = List::create(Named("zeta") = zeta);
    return out;
  }
  
  else if((lambda >= 1) & (!is_na_pooled_var_MS_A) & (!is_na_pooled_var_MS_B)){
    
    NumericVector x = MS_B;
    NumericVector y = MS_A;
    
    float mx = mean_MS_B;
    float my = mean_MS_A;
    
    float sxx = 0;
    float sxy = 0;
    float sse = 0;
    for(int i = 0; i < N; ++i){
      bool na_check_x = ISNAN(x[i]);
      bool na_check_y = ISNAN(y[i]);
      if((!na_check_x) & (!na_check_y)){
        sxx = sxx + pow(x[i] - mx, 2);
        sxy = sxy + (x[i] - mx) * (y[i] - my);
      }
    
    }
    
    if(silence == -1){
      Rcout << "----- Progress (2 / 3) Calculations of sum of squares ------------" << "\n";
      Rcout << "## Sum of squares of x is : " << sxx << "\n";
      Rcout << "## Sum of squares of x and y is : " << sxy << "\n";
      Rcout << "------------------------------------------------------------------" << "\n";  
    }
    
    
    float b1 = sxy / sxx;
    float b0 = my - b1 * mx;
    
    int effective_N = N;
    for(int i = 0; i < N; ++i){
      bool na_check_x = ISNAN(x[i]);
      bool na_check_y = ISNAN(y[i]);
      if((!na_check_x) & (!na_check_y)){
        float yhat = b0 + b1 * x[i];  
        sse = sse + pow(y[i] - yhat, 2);
      }
      else{
        effective_N = effective_N - 1;
        continue;
      }
    }
    float mse = sse / (effective_N - 2);
    float varpar = mse * (effective_N + 2);
    varpar = varpar / effective_N;
    
    if(silence == 0 or silence == -1){
      Rcout << "----- Progress (3 / 3) remaining calculations --------------------" << "\n";
      Rcout << "## Estimtated slope estimator : " << b1 << "\n";
      Rcout << "## Estimtated intercept estimator : " << b0 << "\n";
      Rcout << "## Estimated sum of squares error : " << sse << "\n";
      Rcout << "## Estimated mean sum of squares error : " << mse << "\n";
      Rcout << "## Estimated SD_R squared : " << varpar << "\n";
      Rcout << "## Effective N (after removing NA values) : " << effective_N << "\n";
      Rcout << "------------------------------------------------------------------" << "\n";  
    }
    
    float zeta = varpar / (var_MS_B * pow(b1, 2) + var_MS_A);
    List out = List::create(Named("zeta") = zeta);
    return out;
  }
  
  else if(is_na_pooled_var_MS_A & is_na_pooled_var_MS_B){
    if(silence == -2){
      warning("Pooled variance of both MS_A and MS_B are NA values. Zeta could not be calculated!");  
    }
    List out = List::create(Named("zeta") = NA_REAL);
    return out;  
  }
  
  else if(is_na_pooled_var_MS_A | is_na_pooled_var_MS_A){
    if(silence == -2){
      warning("Pooled variance of both MS_A and MS_B are NA values. Zeta could not be calculated!");  
    }
    List out = List::create(Named("zeta") = NA_REAL);
    return out;  
  }
  List out = List::create(Named("zeta") = NA_REAL);
  return out;
  
}
