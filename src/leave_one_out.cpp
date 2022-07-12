#include <Rcpp.h>
using namespace Rcpp;

//' Leave-one-out on clustered EQA clinical sample data
//' 
//' @title Leave-one-out on clustered EQA clinical sample data
//' @name leave_one_out
//' 
//' @param data \code{List} or \code{data table} - Data with list elements or data table columns with names \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}. \code{SampleID} and \code{ReplicateID} must be of character type for the function to operate correctly
//' @param loo_id Integer - Which of the samples of \code{SampleID} should be left out. Default value is 1
//'
//' @description Needed to calculate jack knife estimates of a parameter, that is required when using \code{BCa_bootstrap_ci()}. Alternatively one could use the Jack knife estimates to calculate standard error or bias of the estimator of relevance
//'
//' @details loo_ids can not be vectorized directly in R. Use \code{sapply()} or \code{replicate()} to leave out sample IDs one by one 
//'
//' @return A \code{list} containing the original data, but without the sample id corresponding to the given \code{loo_id} 
//'
//' @examples \dontrun{
//'   library(commutability.selectivity)
//'   data <- sdwdnsp2()
//'   loo_data <- leave_one_out(data, loo_id = 5)
//' }


// [[Rcpp::export]]
List leave_one_out(List data, int loo_id = 1){
  
  // Definitions
  CharacterVector samples = data["SampleID"];
  CharacterVector replicates = data["ReplicateID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
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
  NumericVector new_MS_A(output_size);
  NumericVector new_MS_B(output_size);
  
  
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
    new_MS_A[j] = MS_A[new_indicies[j]];
    new_MS_B[j] = MS_B[new_indicies[j]];
  }
  List output = List::create(Named("SampleID") = new_samples, Named("ReplicateID") = new_replicates, Named("MP_B") = new_MS_B, Named("MP_A") = new_MS_A);
  return output;
}

