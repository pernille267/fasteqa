#include <Rcpp.h>
using namespace Rcpp;

//' Resample fun-of-replicates data
//' 
//' @title Resample fun-of-replicates data
//' @name resample_fun_data
//' 
//' @param data \code{list} or \code{data table} - Must contain \code{SampleID}, \code{MP_A} and \code{MP_B}. \code{SampleID} and \code{ReplicateID} must be of character type, or else an error is thrown
//' @param make_unique - \code{integer} that controls the output SampleID. If \code{make_unique = 1}, the default, new SampleIDs will be made. Conversely \code{make_unique = 0} will use the original SampleIDs. The latter is not recommended because potential calculations are affected by this 
//'
//' @description In order to construct bootstrap confidence intervals and do inference on a set of population parameters where the underlying distribution is complex, will require bootstrap replicates. This function is both efficient and does its job, but at a cost of a strict input requirement
//'
//' @details Combine with e.g., \code{predict_eqa()}, to estimate classification rates and more
//'
//' @return A list containing the resampled fun-of-replicates data. Use \code{setDT()} for maximum efficiency if you desire to convert the resampled data to a data table
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   data <- simulate_data_eqa(list(n = 25, R = 3, cvx = 0.06, cvy = 0.04))
//'   data$SampleID <- as.character(data$SampleID)
//'   data$ReplicateID <- as.character(data$ReplicateID)
//'   data <- fun_of_replicates(data, fun = "median")
//'   resampled_data <- resample_fun_data(data)
//' }


// [[Rcpp::export]]
List resample_fun_data(List data, int make_unique = 1) {
  CharacterVector SampleID = data["SampleID"];
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  int n = MS_A.size();
  CharacterVector resampled_SampleID = sample(SampleID, n, true);
  IntegerVector new_indices(n);
  CharacterVector new_SampleID(n);
  CharacterVector uni_SampleID(n);
  NumericVector new_MS_A(n);
  NumericVector new_MS_B(n);
  
  int c = 0;
  for(int i = 0; i < n; ++i){
    CharacterVector to_match(1);
    to_match[0] = resampled_SampleID[i];
    LogicalVector matches = in(SampleID, to_match);
    int number_of_matches = sum(matches);
    IntegerVector indices(number_of_matches);
    int k = 0;
    for(int j = 0; j < n; ++j){
      if(matches[j] == 1){
          indices[k] = j;
          ++k;
      }
    }
    for(int j = 0; j < number_of_matches; ++j){
      new_indices[c] =  indices[j];
      new_SampleID[c] = SampleID[indices[j]];
      String uni(c + 1);
      uni_SampleID[c] = uni;
      ++c;
    }
    
  }
  for(int i = 0; i < n; ++i){
    new_MS_A[i] = MS_A[new_indices[i]];
    new_MS_B[i] = MS_B[new_indices[i]];
  }
  
  if(make_unique == 1){
    List output = List::create(Named("SampleID") = uni_SampleID, Named("MP_A") = new_MS_A, Named("MP_B") = new_MS_B);
    return output;
  }
  else if(make_unique == 0){
    List output = List::create(Named("SampleID") = new_SampleID, Named("MP_A") = new_MS_A, Named("MP_B") = new_MS_B);
    return output;
  }
  List output = List::create(Named("SampleID2") = new_SampleID, Named("SampleID") = uni_SampleID, Named("MP_A") = new_MS_A, Named("MP_B") = new_MS_B);
  return output;
}


