#include <Rcpp.h>
using namespace Rcpp;

//' Estimate relative non-selectivity measure k based on data
//' 
//' @title Estimate MS precision measures
//' @name Estimatek
//' @param data A data table or data frame converted to a list so that each column is a unique list element. The ID columns must be of character type and contain SampleID, ReplicateID, MP_A and MP_B 
//' @param silence An integer (-1, 0, 1) that controls how much is printed. 0 shows some messages regarding temporary output, whereas -1 shows even more messages for debugging. -1 is used for superusers.
//' 
//' @description Estimate the relative non-selectivity measure k, that is the ratio of the pooled prediction error variance and the sum of repeatability variances.  
//' 
//' @details Differences in non-selectivity between measurement systems may cause problems in e.g., evaluation of commutability. A large value of k indicate that we have have large differences in selectivity
//' 
//' @return A single estimate of k (in a list) based on the specific grouped data
//'
//' @examples \dontrun{
//'   print("Estimate k")
//' }
//'

// [[Rcpp::export]]
List Estimatek(List data, int silence = 1) {
  // Extracting vectors from list
  CharacterVector SampleID = data["SampleID"];
  CharacterVector ReplicateID = data["ReplicateID"];
  int nR = SampleID.size();
  CharacterMatrix IDs ( nR, 2 );
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
  NumericVector VA(n);
  NumericVector VB(n);
  NumericVector MA(n);
  NumericVector MB(n);
  IntegerVector bounds = cumsum(Ri);
  if(silence == -1){
    Rcout << "The value of bounds is: " << bounds << "\n";  
  }
  for(int i = 0; i < n; ++i){
    if(i == 0){
      NumericMatrix sub = Measurements( Range(0, bounds[i] - 1), Range(0, 1) );
      NumericVector sub_A = sub( _, 0);
      NumericVector sub_B = sub( _, 1);
      VA[i] = var(sub_A);
      VB[i] = var(sub_B);
      MA[i] = mean(sub_A);
      MB[i] = mean(sub_B);
      continue;
    }
    NumericMatrix sub = Measurements( Range( bounds[i-1], bounds[i] - 1), Range(0,1) );
    NumericVector sub_A = sub( _, 0);
    NumericVector sub_B = sub( _, 1);
    VA[i] = var(sub_A);
    VB[i] = var(sub_B);
    MA[i] = mean(sub_A);
    MB[i] = mean(sub_B);
  }
  float Var_A = mean(VA);
  float Var_B = mean(VB);
  float mean_A = mean(MA);
  float mean_B = mean(MB);
  float lambda = Var_A / Var_B;
  if(silence == 0){
    Rcout << "The value of Var_A is: " << Var_A << "\n";
    Rcout << "The value of Var_B is: " << Var_B << "\n";
    Rcout << "The value of mean_B is: " << mean_B << "\n";
    Rcout << "The value of mean_A is: " << mean_A << "\n";
    Rcout << "The value of lambda is: " << lambda << "\n";  
  }
  // Calculate remaining components
  if(lambda < 1){
    NumericVector x = MS_A;
    NumericVector y = MS_B;
    
    float mean_x = mean_A;
    float mean_y = mean_B;
    
    float sxx = 0;
    float sxy = 0;
    float sse = 0;
    
    for(int i = 0; i < nR; ++i){
      sxx += (x[i] - mean_x) * (x[i] - mean_x);
      sxy += (x[i] - mean_x) * (y[i] - mean_y);
    }
    
    float b1 = sxy / sxx;
    float b0 = mean_y - b1 * mean_x;
    if(silence == 0){
      Rcout << "The value of b1 is: " << b1 << "\n";
      Rcout << "The value of b0 is: " << b0 << "\n";  
    }
    NumericVector yhat(nR);
    for(int i = 0; i < nR; ++i){
      yhat[i] = b0 + b1 * x[i];
      sse += (y[i] - yhat[i]) * (y[i] - yhat[i]);
    }
    
    float varpar = (sse / (nR - 2)) * ((nR + 2) / nR);
    float k = varpar / (Var_A * b1 * b1 + Var_B);
    List out = List::create(Named("k") = k);
    return out;
  }
  else{
    NumericVector x = Measurements( _, 0 );
    NumericVector y = Measurements( _, 1 );
    
    float mean_x = mean_B;
    float mean_y = mean_A;
    float sxx = 0;
    float sxy = 0;
    float sse = 0;
    for(int i = 0; i < nR; ++i){
      sxx += (x[i] - mean_x) * (x[i] - mean_x);
      sxy += (x[i] - mean_x) * (y[i] - mean_y);
    }
    float b1 = sxy / sxx;
    float b0 = mean_y - b1 * mean_x;
    if(silence == 0){
      Rcout << "The value of b1 is: " << b1 << "\n";
      Rcout << "The value of b0 is: " << b0 << "\n"; 
    }
    NumericVector yhat(nR);
    for(int i = 0; i < nR; ++i){
      yhat[i] = b0 + b1 * x[i];
      sse += (y[i] - yhat[i]) * (y[i] - yhat[i]);
    }
    float varpar = (sse / (n - 2)) * ((n + 2) / n);
    float k = varpar / (Var_B * b1 * b1 + Var_A);
    List out = List::create(Named("k") = k);
    return out;
  }
}
