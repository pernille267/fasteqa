#include <Rcpp.h>
using namespace Rcpp;

//' Simulation of EQAM data based on study design, that is, n and R
//'
//' @title Simulate EQAM data for any combination of n and R
//' @name SimulateEqaData
//' 
//' @param n The number of clinical sampels to generate. Must be an integer
//' @param R The number of replicates on each clinical sample. Must be an integer
//' @param silence Should progress reports be outputted. Defaults to no, which is \code{silence = 1}
//' 
//' @description Simulates a data set wit n x R rows, and 4 columns, where the two first is the base ID columns (SampleID and ReplicateID), and the last it the base Measurement columns (MP_A and MP_B). Enjoy
//' @details In order to glue the elements of the output you may use \code{as.data.frame()}, \code{as.data.table()}, \code{as.tibble()}, etc., but the fastest is \code{setDT()} from the \code{data.table} package.
//'
//' @return A list where each list element is a column of the generated data
//'
//' @examples \dontrun{
//'   library(data.table)
//'   n <- 25
//'   R <- 3
//'   simulated_data <- setDT(SimulateEqaData(n = n, R = R))
//' }
//'

// [[Rcpp::export]]
List SimulateEqaData(int n = 25, int R = 3, int silence = 1){
  NumericVector tau(n);
  int nR = n * R;
  IntegerVector SampleID(nR);
  IntegerVector ReplicateID(nR);
  NumericVector MP_A(nR);
  NumericVector MP_B(nR);
  NumericVector Ci(2);
  float CVX = R::rbeta(2, 5) / 10;
  float CVY = R::rbeta(2, 5) / 10;
  
  Ci[0] = R::runif(5.0, 10.0);
  Ci[1] = Ci[0] + R::runif(10.0, 20.0);
  float Mean = mean(Ci);
  float SDX = Mean * CVX;
  float SDY = Mean * CVY;
  if(silence != 1){
    Rcout << "The value of CVX is: " << CVX << "\n";
    Rcout << "The value of CVY is: " << CVY << "\n";
    Rcout << "The value of Ci is: " << Ci << "\n";
    Rcout << "The value of SDX is: " << SDX << "\n";
    Rcout << "The value of SDY is: " << SDY << "\n";
  }
  for(int i = 0; i < n; ++i){
    tau[i] = R::runif(Ci[0], Ci[1]);
    int lower = i * R;
    int upper = R * (1 + i) - 1;
    IntegerVector idx = Rcpp::Range(lower, upper);
    for(int r = 0; r < R; ++r){
      SampleID[idx[r]] = i + 1;
      ReplicateID[idx[r]] = r + 1;
      MP_A[idx[r]] = tau[i] + R::rnorm(0, SDY);
      MP_B[idx[r]] = tau[i] + R::rnorm(0, SDX);
    }
  }
  List out = List::create(Named("SampleID") = SampleID, Named("ReplicateID") = ReplicateID, Named("MP_A") = MP_A, Named("MP_B") = MP_B);
  return out;
}
  
