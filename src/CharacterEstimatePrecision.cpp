#include <Rcpp.h>
using namespace Rcpp;

//' Estimate MS precision measures (variances, CVs and lambda) for character IDs
//' 
//' @title Estimate MS precision measures
//' @name CharacterEstimatePrecision
//' @param data A data frame or table with format LFDT converted to list so that each column is a list element
//' @param silence - Integer, hide helpful (maybe unhelpful for most people) progress reports. Set to 0 for these progress report to appear. For additional debugging information, set this to \code{-1}.
//' 
//' @description Estimate precision of EQA data when the ID columns are of string type. All ID columns must be character vectors for this function to work
//' 
//' @details These estimates are used as part of calculation of relative selectivity, and construction of prediction intervals. They can also be used to visualize the the MS precisions using the much used CVs in percent
//' 
//' @return A list containing estimated MS precision estimates for one group. Can be melted together by using \code{setDT()} from the \code{data.table} package
//' 
//' @examples \dontrun{
//'   print(1)
//' }
//' 

// [[Rcpp::export]]
List CharacterEstimatePrecision(List data, int silence = 1) {
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
  for(int i = 1; i <= nR; ++i){
    if(i == nR){
      Ri[idx] = k;
      ++idx;
      break;
    }
    if(SampleID[i] == SampleID[i - 1]){
      ++k;
    }
    else{
      Ri[idx] = k;
      k = 1;
      ++idx;
    }
  }
  NumericVector Var_A(n);
  NumericVector Var_B(n);
  NumericVector mean_A(n);
  NumericVector mean_B(n);
  IntegerVector bounds = cumsum(Ri);
  for(int i = 0; i < n; ++i){
    if(i == 0){
      NumericMatrix sub = Measurements( Range(0, bounds[i] - 1), Range(0,1) );
      NumericVector sub_A = sub( _, 0);
      NumericVector sub_B = sub( _, 1);
      Var_A[i] = var(sub_A);
      Var_B[i] = var(sub_B);
      mean_A[i] = mean(sub_A);
      mean_B[i] = mean(sub_B);
      continue;
    }
    NumericMatrix sub = Measurements( Range( bounds[i-1], bounds[i] - 1), Range(0,1) );
    NumericVector sub_A = sub( _, 0);
    NumericVector sub_B = sub( _, 1);
    Var_A[i] = var(sub_A);
    Var_B[i] = var(sub_B);
    mean_A[i] = mean(sub_A);
    mean_B[i] = mean(sub_B);
  }
  
  float out_Var_A = mean(Var_A);
  float out_Var_B = mean(Var_B);
  float out_CV_A = 100 * sqrt(out_Var_A) / mean(mean_A);
  float out_CV_B = 100 * sqrt(out_Var_B) / mean(mean_B);
  float lambda = out_Var_A / out_Var_B;
  List out = List::create(Named("Var_A") = out_Var_A, Named("Var_B") = out_Var_B, Named("CV_A (%)") = out_CV_A, Named("CV_B (%)") =  out_CV_B, Named("lambda") = lambda);
  return out;
}
