#include <Rcpp.h>
using namespace Rcpp;

//' Resample clustered EQA clinical sample data
//' 
//' @title Resample clustered EQA clinical sample data
//' @name resample_samples
//' 
//' @param data A list or a data table with elements/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}. \code{SampleID} and \code{ReplicateID} must be of character type for the function to operate correctly
//' @param silence An integer that controls the progress reports outputted for debugging and further examination of the command. \code{silence = -1} or \code{silence = 0} signify that progress reports should be printed. Default is \code{silence = 1} which suppresses all printing
//'
//' @description In order to construct bootstrap confidence intervals and do inference on a set of population parameters where the underlying distribution is complex, will require resample of clustered data. This function is both efficient and does its job, but at a cost of a strict input requirement
//'
//' @details \code{resample_samples()} is a very efficient algorithm to resample EQA. Combine with \code{Estimatek()} to resample k or combine with \code{CharacterEstimatePrecision()} to resample variability measures such as CVs and variances. May also be combined with other functions
//'
//' @return A list containing the resampled EQA clinical sample data. Use \code{setDT()} for maximum efficiency if needed to convert the resampled data to a data table
//'
//' @examples \dontrun{
//'   library(commutability.selectivity)
//'   data <- sdwdnsp2()
//'   data$SampleID <- as.character(data$SampleID)
//'   data$ReplicateID <- as.character(data$ReplicateID)
//'   resampled_data <- resample_samples(data)
//' }


// [[Rcpp::export]]
List resample_samples(List data, int silence = 1){
  CharacterVector samples = data["SampleID"];
  CharacterVector replicates = data["ReplicateID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  CharacterVector unique_samples = unique(samples);
  int n = unique_samples.size();
  int N = samples.size();
  CharacterVector resampled_samples = sample(unique_samples,
                                             unique_samples.size(),
                                             true);
  int output_size = 0;
  if(silence == 0){
    Rcout << "there are recorded " << n << "unique clinical samples" << "\n";
    Rcout << "there are recorded " << N << "measurement results" << "\n";
  }
  
  // How large should the output be?
  for(int i = 0; i < n; ++i){
    CharacterVector to_match(1);
    to_match[0] = resampled_samples[i];
    LogicalVector matching_resamples = in(samples, to_match);
    int index_vector_length = sum(matching_resamples);
    output_size += index_vector_length;
  }
  
  if(silence == 0){
    Rcout << "the size of the output is expected to be : " << output_size << "\n";
  }
  
  // Defining output-helping vectors to be filled in the loop below
  IntegerVector new_indicies(output_size);
  CharacterVector new_samples(output_size);
  CharacterVector new_replicates(output_size);
  NumericVector new_MS_A(output_size);
  NumericVector new_MS_B(output_size);
  int counter = 0;
  
  // Fill output-helping vectors
  for(int i = 0; i < n; ++i){
    CharacterVector to_match(1);
    to_match[0] = resampled_samples[i];
    LogicalVector matching_resamples = in(samples, to_match);
    int index_vector_length = sum(matching_resamples);
    NumericVector index_vector(index_vector_length);
    int relevant_ind = 0;
    for(int j = 0; j < N; ++j){
      if(matching_resamples[j] == 1){
        index_vector[relevant_ind] = j;
        ++relevant_ind;
      }
    }
    for(int k = 0; k < index_vector_length; ++k){
      new_indicies[counter] = index_vector[k];
      new_samples[counter] = unique_samples[i];
      ++counter;
    }
  }
  
  if(silence == 0){
    Rcout << "the number of newly selected values are : " << counter << "\n";
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

