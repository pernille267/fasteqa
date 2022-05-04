#include <Rcpp.h>
using namespace Rcpp;

//' Calculates means of replicate
//' 
//' @title Calculate means of replicates
//' @name MeanOfReplicates
//' 
//' @param data A data frame or data table converted to a list so that every column is unique list elements
//' @param silence should progress reports be outputted? For debugging and additional information, where you set \code{silence} to either -1 or 0. Default is 1, which means that no progress reports should be outputted
//'
//' @description Constructing prediction intervals require EQA data to be on mean-of-replicates-form. So this function should be applied to the EQAM data before it is included into the \code{PredictDeming()} or \code{PredictOLS()}
//'
//' @details The difference between this function and \code{mean_of_replicates()} method in the \code{commutability.selectivity}, is that this is more than ten times faster.
//'
//' @return A list containing the elements that is needed to build the resulting data of results. Use \code{setDT()} for maximum efficiency when converting the list into a data table 
//'
//' @examples \dontrun{
//'   library(commutability.selectivity)
//'   library(data.table)
//'   data <- sampled_cs_measurements
//'   data <- MS_wise(data = data)[Comparison=="MP1 - MP2",]
//'   data$SampleID <- as.character(data$SampleID)
//'   data$ReplicateID <- as.character(data$ReplicateID)
//'   input <- as.list(data)
//'   output <- setDT(MeanOfReplicates(data = input, silence = 1))
//' }

// [[Rcpp::export]]
List MeanOfReplicates(List data, int silence = 1) {
  
  // Extracting vectors from list (IDs)
  CharacterVector SampleID = data["SampleID"];
  CharacterVector ReplicateID = data["ReplicateID"];
  int nR = SampleID.size();
  CharacterMatrix IDs ( nR, 2 );
  // Extracting vectors from list (Measurements)
  NumericVector MS_A = data["MP_A"];
  NumericVector MS_B = data["MP_B"];
  NumericMatrix Measurements ( nR, 2 );
  // Filling in matrices
  for(int i = 0; i < nR; ++i){
    IDs ( i, 0 ) = SampleID[i];
    IDs ( i, 1 ) = ReplicateID[i];
    Measurements ( i, 0) = MS_A[i];
    Measurements ( i, 1) = MS_B[i];
  }
  // Summary over SampleID
  int n = Rcpp::unique(SampleID).size();
  int k = 1;
  int idx = 0;
  IntegerVector Ri (n);
  CharacterVector Si (n);
  for(int i = 1; i <= nR; ++i){
    if(i == nR){
      Ri[idx] = k;
      Si[idx] = SampleID[i - 1];
      ++idx;
      break;
    }
    if(SampleID[i] == SampleID[i - 1]){
      ++k;
    }
    else{
      Ri[idx] = k;
      Si[idx] = SampleID[i-1];
      k = 1;
      ++idx;
    }
  }
  NumericVector mean_A(n);
  NumericVector mean_B(n);
  IntegerVector bounds = cumsum(Ri);
  for(int i = 0; i < n; ++i){
    if(i == 0){
      NumericMatrix sub = Measurements( Range(0, bounds[i] - 1), Range(0,1) );
      NumericVector sub_A = sub( _, 0);
      NumericVector sub_B = sub( _, 1);
      mean_A[i] = mean(sub_A);
      mean_B[i] = mean(sub_B);
      continue;
    }
    NumericMatrix sub = Measurements( Range( bounds[i-1], bounds[i] - 1), Range(0,1) );
    NumericVector sub_A = sub( _, 0);
    NumericVector sub_B = sub( _, 1);
    mean_A[i] = mean(sub_A);
    mean_B[i] = mean(sub_B);
  }
  List out = List::create(Named("SampleID") = Si, Named("MP_A") = mean_A, Named("MP_B") = mean_B);
  return out;
}
