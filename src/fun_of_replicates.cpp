#include <Rcpp.h>
using namespace Rcpp;

//' Apply a mathematical summary function on every SampleID
//' 
//' @title Apply a mathematical summary function on every SampleID
//' @name fun_of_replicates
//' 
//' @param data \code{list} or \code{data table} - Data with elements/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
//' @param fun \code{string} - Used to choose the summary function one want to apply on each samples' measurements. Here is the list of possible valid summary functions:
//' \itemize{
//'   \item{\code{mean: }}{Taking the average of replicated measurments for samples with at least one measurement. Default}
//'   \item{\code{median: }}{Taking the median of replicated measurments for samples with at least one measurement}
//'   \item{\code{min: }}{Taking the minimum of replicated measurments for samples with at least one measurement}
//'   \item{\code{max: }}{Taking the maximum of replicated measurements for samples with at least one measurement}
//'   \item{\code{var: }}{Taking the variance of replicated measurements for samples with at least two measurements}
//'   \item{\code{sd: }}{Taking the standard deviation of replicated measurements for samples with at least two measurements}
//'   \item{\code{cv: }}{Taking the coefficient of variation of replicated measurements for samples with at least two measurements}
//' }
//' @param silence \code{integer} - Should debugging messages be printed. Default is no
//'
//' @description A practical function to evaluate summary typical summary functions of samples' replicated measurements. Taking mean of replicates, sd of replicates or cv of replicates are typical when analyzing EQA data so these are very useful. The remaning functions are not used as much, but may be needed in some cases so they are included based on this fact. In order for this function to work propery you need to ensure that:
//' \itemize{
//'   \item{\code{SampleID: }}{Must be a character vector}
//'   \item{\code{ReplicateID: }}{Must be a character vector}
//'   \item{\code{MP_A: }}{Must be a numeric vector}
//'   \item{\code{MP_B: }}{Must be a numeric vector}
//' }
//'
//' @details The difference between this function and \code{mean_of_replicates()} method in the \code{commutability.selectivity}, is that this is more than ten times faster.
//'
//' @return list - data containing the elements that is needed to build the resulting data of results. Use \code{setDT()} for maximum efficiency when converting the list into a data table 
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   data <- simulate_eqa_data()
//'   mean_of_replicates_data <- fun_of_replicates(data)
//'   var_of_replicates_data <- fun_of_replicates(data, "var")
//' }


// [[Rcpp::export]]
List fun_of_replicates(List data, String fun = "mean", int silence = 1){
  CharacterVector SampleID = data["SampleID"];
  CharacterVector ReplicateID = data["ReplicateID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  CharacterVector summary_SampleID = unique(SampleID);
  int n = summary_SampleID.size();
  int N = SampleID.size();
  NumericVector new_MS_A(n);
  NumericVector new_MS_B(n);
  int counter = 0;
  int replicate_number_requirement = 0;
  CharacterVector valid_funs = CharacterVector::create("mean", "median", "min", "max", "var", "sd", "cv");
  CharacterVector fun_check(1);
  fun_check[0] = fun;
  int valid_fun_input = sum(in(fun_check, valid_funs));
  if(valid_fun_input < 1){
    stop("Given input for fun is not recongnized. See ?fun_of_replicates() for valid fun inputs");
  }
  if(valid_fun_input > 1){
    stop("Vector of functions are not allowed. See ?fun_of_replicates() for valid fun inputs");
  }
  
  if(fun == "mean" or fun == "median" or fun == "min" or fun == "max"){
    replicate_number_requirement = 1;
  }
  if(fun == "sd" or fun == "var" or fun == "cv"){
    replicate_number_requirement = 2;
  }
  
  if(fun == "mean"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = mean(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = mean(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  // Median
  //
  // // // // // // //
  
  if(fun == "median"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = median(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = median(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  // MIN
  //
  // // // //
  
  if(fun == "min"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = min(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = min(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  // MAX
  //
  // // // //
  
  if(fun == "max"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = max(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = max(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  // VAR
  //
  // // // //
  
  if(fun == "var"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = var(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = var(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  // SD
  //
  // // // //
  
  if(fun == "sd"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = sd(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = sd(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  // CV
  //
  // // // //
  
  if(fun == "cv"){
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
        new_MS_A[counter] = NA_REAL;  
      }
      if(ith_sample_measurements_A.size() >= replicate_number_requirement){
        new_MS_A[counter] = sd(ith_sample_measurements_A) / mean(ith_sample_measurements_A);
      }
      if(ith_sample_measurements_B.size() < replicate_number_requirement){
        new_MS_B[counter] = NA_REAL;
      }
      if(ith_sample_measurements_B.size() >= replicate_number_requirement){
        new_MS_B[counter] = sd(ith_sample_measurements_B) / mean(ith_sample_measurements_B);  
      }
      ++counter;
    }
  }
  
  List out = List::create(Named("SampleID") = summary_SampleID, Named("MP_A") = new_MS_A, Named("MP_B") = new_MS_B);
  return out;
}


