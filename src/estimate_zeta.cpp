#include <Rcpp.h>
using namespace Rcpp;

//' @title Quantify Differences in Non-selectivity Using \eqn{\hat{\zeta}}
//' @name estimate_zeta_ols
//' @param data A \code{list} or \code{data.table}. Must contain:
//'             \itemize{
//'               \item \code{SampleID: } A \code{character} vector. The
//'                     clinical sample identifiers.
//'               \item \code{ReplicateID: } A \code{character} vector. The
//'                     replicate measurement identifiers.
//'               \item \code{MP_A: } A \code{numeric} vector. The measurement
//'                     results from IVD-MD \code{MP_A} (response).
//'               \item \code{MP_B: } A \code{numeric} vector. The measurement
//'                     results from IVD-MD \code{MP_B} (predictor).
//'               
//'             }
//' @param silence An \code{integer}. Controls verbose. Verbose can be either
//'                informative or just noise. Keep in mind that printing to
//'                the console will slow down calcuation drastically.
//'                \itemize{
//'                   \item \code{1: } All progress reports are silenced. This is the default.
//'                   \item \code{0: } Estimation steps and temporary results are printed to the console
//'                   \item \code{< 0: } Debugging. Expert use only.
//'                }
//'                
//' @description
//' Estimates the magnitude of differences in non-selectivity (DINS) in a
//' IVD-MD comparison between \code{MP_A} and \code{MP_B}. The parameter
//' \eqn{\zeta} measures DINS and is estimated using \eqn{\hat{\zeta}}. See
//' details.
//' 
//' @details
//' Differences in non-selectivity (DINS) between in vitro diagnostic medical
//' devices (IVD-MDs) may cause problems in e.g., evaluation of commutability
//' of external quality assessment materials or certified reference materials. 
//' A large value of \eqn{\hat{\zeta}} may indicate an unacceptable magnitude of
//' DINS between compared IVD-MDs. The estimator is defined by
//' 
//' \eqn{\hat{\zeta} = \frac{S^2 \cdot (nR + 2) / nR}
//' {\hat{\sigma}_v^2 + b_1^2 \cdot \hat{\sigma}_h^2}},
//' 
//' where:
//' \itemize{
//'   \item \eqn{S^2}: The estimated residual variance of the regression model
//'   \item \eqn{b_1^2}: The square of the slope estimator of the regression
//'         model.
//'   \item \eqn{\hat{\sigma}_v^2}: The estimated pooled repeatability variance
//'         of \code{MP_A}.
//'   \item \eqn{\hat{\sigma}_h^2}: The estimated pooled repeatability variance
//'         of \code{MP_B}.
//'   \item \eqn{n}: The number of unique clinical samples in \code{data}
//'   \item \eqn{R}: A suitable statistic (e.g., mean) of the number of
//'         replicated measurements performed on each clinical sample in
//'         \code{data}.
//' }
//' 
//' See references for finer details. This estimator of \eqn{\zeta} utilizes
//' the ordinary least squares regression model. Thus, it is sensitive to
//' non-linearity, which may or may not be due to systematic differences in
//' non-selectivity.
//' 
//' To mitigate model endogeneity effects from measurement error in the
//' predictor, the roles of \code{MP_A} and \code{MP_B} shifts if
//' 
//' \eqn{\lambda = \hat{\sigma}_v^2 / \hat{\sigma}_h^2 < 0.5}.
//' 
//' Be careful using this estimator of \eqn{\zeta} if both
//' \eqn{\hat{\sigma}_v^2} and \eqn{\hat{\sigma}_h^2} are large and the domain of
//' the measurement results is very narrow.
//' 
//' Note: If the relationship between \code{MP_A} and \code{MP_B} is
//' non-linear, it is advisable to use \code{estimate_zeta_ss()} from the
//' \code{smooth.commutable} package instead.
//' 
//' @return
//' A \code{list} of length one. Contains \code{zeta}, which is a
//' \code{double}. This is the calculated \eqn{\hat{\zeta}}.
//'
//' @examples
//' # Required packages
//' library(fasteqa)
//' 
//' # Estimate zeta based on the raw data
//' zeta <- estimate_zeta_ols(test_data)$zeta
//' 
//' # The output
//' print(round(zeta, 2L))
//' 
//' # Log-transformed data
//' log_test_data <- test_data
//' log_test_data$MP_A <- log(log_test_data$MP_A)
//' log_test_data$MP_B <- log(log_test_data$MP_B)
//' 
//' # Estimate zeta based on the log-transformed data
//' zeta_log <- estimate_zeta_ols(log_test_data)$zeta
//' 
//' # The output
//' print(round(zeta_log, 2L))
//'
//' @references
//' Fauskanger P.K., et al. (2025) Quantification of Difference in
//' Nonselectivity Between In Vitro Diagnostic Medical Devices.
//' \emph{Biometrical Journal}. 67: e70032.
//' \url{https://doi.org/10.1002/bimj.70032}
//'

// [[Rcpp::export]]
List estimate_zeta_ols(List data, int silence = 1) {
 
 // Extract components from data
 CharacterVector SampleID = data["SampleID"];
 CharacterVector ReplicateID = data["ReplicateID"];
 NumericVector MP_A = data["MP_A"];
 NumericVector MP_B = data["MP_B"];
 
 // Unique SampleIDs
 CharacterVector summary_SampleID = unique(SampleID);
 
 // Sample sizes
 int n = summary_SampleID.size();
 int N = SampleID.size();
 
 // Allocate sample variance vector
 NumericVector ith_var_MP_A(n);
 NumericVector ith_var_MP_B(n);
 
 // Requirements
 int replicate_number_requirement = 2;
 
 // Counter variable
 int c = 0;
 
 // Look at Sample by Sample
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
   // Otherwise, the kth measurement will be the kth measurement of the corresponding MP_A or MP_B
   for(int k = 0; k < number_of_matches; ++k){
     bool is_na_A = ISNAN(MP_A[indices[k]]);
     bool is_na_B = ISNAN(MP_B[indices[k]]);
     if(!is_na_A){
       ith_sample_measurements_A[k] = MP_A[indices[k]];
     }
     else{
       ith_sample_measurements_A[k] = -100.0;
     }
     if(!is_na_B){
       ith_sample_measurements_B[k] = MP_B[indices[k]];  
     }
     else{
       ith_sample_measurements_B[k] = -100.0;
     }
   }
   
   // Create NA-search vectors with zeros and ones matching with ith sample of MP_A and MP_B
   IntegerVector NA_search_A(number_of_matches);
   IntegerVector NA_search_B(number_of_matches);
   
   // Recall: ith_sample_measurements_*[k] = -100.0 signify that the value is a NA-value
   // NA_search_* will return 0 if ith_sample_measurements_*[k] = -100.0, and 1 otherwise 
   for(int k = 0; k < number_of_matches; ++k){
     if(ith_sample_measurements_A[k] == -100.0){
       NA_search_A[k] = 0;
     }
     if(ith_sample_measurements_A[k] > -100.0){
       NA_search_A[k] = 1;
     }
     if(ith_sample_measurements_B[k] == -100.0){
       NA_search_B[k] = 0;
     }
     if(ith_sample_measurements_B[k] > -100.0){
       NA_search_B[k] = 1;
     }
     if((ith_sample_measurements_A[k] < -100.0) | (ith_sample_measurements_B[k] < -100.0)){
       stop("Unrealistic measurements are recorded!");
     }
   }
   
   // Aligning vectors NA search vectors and measurement vectors making it easier to delete NA-values
   ith_sample_measurements_A = ith_sample_measurements_A.sort();
   ith_sample_measurements_B = ith_sample_measurements_B.sort();
   NA_search_A = NA_search_A.sort();
   NA_search_B = NA_search_B.sort();
  
   // Step-wise check for NA values at vector start and delete if it is a NA-value for MP_A
   for(int k = 0; k < number_of_matches; ++k){
     if(NA_search_A[k] == 0){
       ith_sample_measurements_A.erase(0); 
     }
   }
   
   // Step-wise check for NA values at vector start and delete if it is a NA-value for MP_B
   for(int k = 0; k < number_of_matches; ++k){
     if(NA_search_B[k] == 0){
       ith_sample_measurements_B.erase(0);
     }
   }
   
   // Checks if the number of replicates meets the conditions of the relevant summary function
   // If met, the summary function is applied. Otherwise, the result will be a NA-value
   
   if(ith_sample_measurements_A.size() >= replicate_number_requirement){
     ith_var_MP_A[c] = var(ith_sample_measurements_A);
   }
   else if(ith_sample_measurements_A.size() < replicate_number_requirement){
     ith_var_MP_A[c] = NA_REAL;  
   }
   if(ith_sample_measurements_B.size() >= replicate_number_requirement){
     ith_var_MP_B[c] = var(ith_sample_measurements_B);  
   }
   else if(ith_sample_measurements_B.size() < replicate_number_requirement){
     ith_var_MP_B[c] = NA_REAL;
   }
   ++c;
 }
 
 int effective_N_A = N;
 int effective_n_A = n;
 int effective_N_B = N;
 int effective_n_B = n;
 double var_MP_A = 0;
 double var_MP_B = 0;
 
 for(int j = 0; j < n; ++j){
   bool is_na_var_A = ISNAN(ith_var_MP_A[j]);
   bool is_na_var_B = ISNAN(ith_var_MP_B[j]);
   if(!is_na_var_A){
     var_MP_A += ith_var_MP_A[j];
   }
   else if(is_na_var_A){
     effective_n_A = effective_n_A - 1;
   }
   if(!is_na_var_B){
     var_MP_B += ith_var_MP_B[j];
   }
   else if(is_na_var_B){
     effective_n_B = effective_n_B - 1;
   }
 }
 
 if(effective_n_A >= 1){
   var_MP_A = var_MP_A / effective_n_A; 
 }
 else if(effective_n_A < 1){
   var_MP_A = NA_REAL;
 }
 if(effective_n_B >= 1){
   var_MP_B = var_MP_B / effective_n_B;
 }
 else if(effective_n_B < 1){
   var_MP_B = NA_REAL;
 }
 double lambda = 0;
 bool is_na_pooled_var_MP_A = ISNAN(var_MP_A);
 bool is_na_pooled_var_MP_B = ISNAN(var_MP_B);
 bool can_calculate_zeta = (!is_na_pooled_var_MP_A) & (!is_na_pooled_var_MP_B);
 if(can_calculate_zeta){
   if(var_MP_A < 0){
     List out = List::create(Named("zeta") = NA_REAL);
     return out;
   }
   lambda = var_MP_A / var_MP_B;  
 }
 else{
   List out = List::create(Named("zeta") = NA_REAL);
   return out;
 }
 
 double mean_MP_A = 0.0;
 double mean_MP_B = 0.0;
 
 for(int i = 0; i < N; ++i){
   bool is_na_mean_A = ISNAN(MP_A[i]);
   bool is_na_mean_B = ISNAN(MP_B[i]);
   if((!is_na_mean_A) & (!is_na_mean_B)){
     mean_MP_A += MP_A[i];
     mean_MP_B += MP_B[i];
   }
   else{
     effective_N_A = effective_N_A - 1;
     effective_N_B = effective_N_B - 1;
   }
 }
 if(effective_N_A >= 1){
   mean_MP_A = mean_MP_A / effective_N_A;  
 }
 if(effective_N_A < 1){
   mean_MP_A = NA_REAL;
 }
 if(effective_N_B >= 1){
   mean_MP_B = mean_MP_B / effective_N_B;
 }
 if(effective_N_B < 1){
   mean_MP_B = NA_REAL;  
 }
 
 if(silence == 0 or silence == -1){
   Rcout << "----- Progress (1 / 3) Calculations of variances and means -------" << "\n";
   Rcout << "## Pooled variance of MP_A is : " << var_MP_A << "\n";
   Rcout << "## Pooled variance of MP_B is : " << var_MP_B << "\n";
   Rcout << "## Mean of MP_A is : " << mean_MP_A << "\n";
   Rcout << "## Mean of MP_B is : " << mean_MP_B << "\n";
   Rcout << "## lambda is : " << lambda << "\n";
   Rcout << "## N is : " << N << "\n";
   Rcout << "------------------------------------------------------------------" << "\n";  
 }
 
 // Calculate varpar and then zeta
 
 if(lambda < 0.5){
   
   NumericVector x = MP_A;
   NumericVector y = MP_B;
   
   double mx = mean_MP_A;
   double my = mean_MP_B;
   
   double sxx = 0;
   double sxy = 0;
   double sse = 0;
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
   
   double b1 = sxy / sxx;
   double b0 = my - b1 * mx;
   
   int effective_N = N;
   for(int i = 0; i < N; ++i){
     bool na_check_x = ISNAN(x[i]);
     bool na_check_y = ISNAN(y[i]);
     if((!na_check_x) & (!na_check_y)){
       double yhat = b0 + b1 * x[i];  
       sse = sse + pow(y[i] - yhat, 2);
     }
     else{
       effective_N = effective_N - 1;
       continue;
     }
   }
   double mse = sse / (effective_N - 2);
   double varpar = mse * (effective_N + 2);
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
   
   
   double zeta = varpar / (var_MP_A * pow(b1, 2) + var_MP_B);
   
   List out = List::create(Named("zeta") = zeta);
   return out;
 }
 
 else if(lambda >= 0.5){
   
   NumericVector x = MP_B;
   NumericVector y = MP_A;
   
   double mx = mean_MP_B;
   double my = mean_MP_A;
   
   double sxx = 0;
   double sxy = 0;
   double sse = 0;
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
   
   
   double b1 = sxy / sxx;
   double b0 = my - b1 * mx;
   
   int effective_N = N;
   for(int i = 0; i < N; ++i){
     bool na_check_x = ISNAN(x[i]);
     bool na_check_y = ISNAN(y[i]);
     if((!na_check_x) & (!na_check_y)){
       double yhat = b0 + b1 * x[i];  
       sse = sse + pow(y[i] - yhat, 2);
     }
     else{
       effective_N = effective_N - 1;
       continue;
     }
   }
   double mse = sse / (effective_N - 2);
   double varpar = mse * (effective_N + 2);
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
   
   double zeta = varpar / (var_MP_B * pow(b1, 2) + var_MP_A);
   List out = List::create(Named("zeta") = zeta);
   return out;
 }
 
 List out = List::create(Named("zeta") = NA_REAL);
 return out;
 
}

//' Estimate differences in non-selectivity with zeta using Deming regression
//' 
//' @title Estimate differences in non-selectivity with zeta using Deming regression
//' @name estimate_zeta_deming
//' @param data \code{list} or \code{data table} - Data with elements/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
//' @param silence \code{integer} - How much progress reports should be returned. Note that returning progress reports will slow down the performance drastically. There are three valid inputs:
//' \itemize{
//'   \item \code{1: } All progress reports are silenced. This is the default.
//'   \item \code{0: } Estimation steps and temporary results are printed to the console.
//'   \item \code{< 0: } Debugging. Expert use only.
//' }
//' 
//' @description Estimate the degree of differences in non-selectivity with zeta. Zeta is is the ratio of the pooled average prediction error variance and the sum of analytical variances.  
//' 
//' @details Differences in non-selectivity between measurement systems may cause problems in e.g., evaluation of commutability. A large value of zeta indicates that we have have large differences in non-selectivity between compared measurement systems. An upper limit of acceptable zeta may be determined based on the allowable increase in prediction interval width and analyte of relevance
//' 
//' @return A list with the point estimate of zeta. The zeta value is a double value, meaning that the precision is 1e-6 (six decimals precision).
//'
//' @examples
//' 
//' # Estimate zeta for a particular dataset
//' zeta <- estimate_zeta_deming(test_data)$zeta
//' print(round(zeta, 2L))
//' 
//' # Estimate zeta based on log-transformed data
//' log_test_data <- test_data
//' log_test_data$MP_A <- log(log_test_data$MP_A)
//' log_test_data$MP_B <- log(log_test_data$MP_B)
//' zeta_log <- estimate_zeta_deming(log_test_data)$zeta
//' print(round(zeta_log, 2L))
//' 
 
 
 // [[Rcpp::export]]
 List estimate_zeta_deming(List data, int silence = 1) {
   
   // Extract components from data
   CharacterVector SampleID = data["SampleID"];
   CharacterVector ReplicateID = data["ReplicateID"];
   NumericVector MP_A = data["MP_A"];
   NumericVector MP_B = data["MP_B"];
   CharacterVector summary_SampleID = unique(SampleID);
   int n = summary_SampleID.size();
   int N = SampleID.size();
   NumericVector ith_var_MP_A(n);
   NumericVector ith_var_MP_B(n);
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
     // If ISNAN results in 1 (i.e., NA-value) let the kth measurement be zero
     // Otherwise, the kth measurement will be the kth measurement of MP_A / MP_B
     for(int k = 0; k < number_of_matches; ++k){
       bool is_na_MP_A = ISNAN(MP_A[indices[k]]);
       bool is_na_MP_B = ISNAN(MP_B[indices[k]]);
       if(!is_na_MP_A){
         ith_sample_measurements_A[k] = MP_A[indices[k]];
       }
       else if(is_na_MP_A){
         ith_sample_measurements_A[k] = -100.0;
       }
       if(!is_na_MP_B){
         ith_sample_measurements_B[k] = MP_B[indices[k]];  
       }
       else if(is_na_MP_B){
         ith_sample_measurements_B[k] = -100.0;
       }
     }
     
     // Create NA-search vectors with zeros and ones matching with ith sample of MP_A and MP_B
     IntegerVector NA_search_A(number_of_matches);
     IntegerVector NA_search_B(number_of_matches);
     
     // recall: ith_sample_measurements_*[k] = 0 signify that the value is a NA-value
     // NA_search_* will return 0 if ith_sample_measurements_*[k] = 0, and 1 otherwise 
     for(int k = 0; k < number_of_matches; ++k){
       if(ith_sample_measurements_A[k] >= -100.0){
         NA_search_A[k] = 1;
       }
       else if(ith_sample_measurements_A[k] < -100.0){
         NA_search_A[k] = 0;
       }
       if(ith_sample_measurements_B[k] >= -100.0){
         NA_search_B[k] = 1;
       }
       else if(ith_sample_measurements_B[k] < -100.0){
         NA_search_B[k] = 0;
       }
     }
     
     // Aligning vectors making it easier to delete NA-values
     ith_sample_measurements_A = ith_sample_measurements_A.sort();
     ith_sample_measurements_B = ith_sample_measurements_B.sort();
     NA_search_A = NA_search_A.sort();
     NA_search_B = NA_search_B.sort();
     
     // Step-wise check for NA values at vector start and delete if it is a NA-value for MP_A
     for(int k = 0; k < number_of_matches; ++k){
       if(NA_search_A[k] == 0){
         ith_sample_measurements_A.erase(0); 
       }
     }
     
     // Step-wise check for NA values at vector start and delete if it is a NA-value for MP_B
     for(int k = 0; k < number_of_matches; ++k){
       if(NA_search_B[k] == 0){
         ith_sample_measurements_B.erase(0);
       }
     }
     
     // Checks if the number of replicates meets the conditions of the relevant summary function
     // If met, the summary function is applied. Otherwise, the result will be a NA-value
     if(ith_sample_measurements_A.size() < replicate_number_requirement){
       ith_var_MP_A[c] = NA_REAL;  
     }
     if(ith_sample_measurements_A.size() >= replicate_number_requirement){
       ith_var_MP_A[c] = var(ith_sample_measurements_A);
     }
     if(ith_sample_measurements_B.size() < replicate_number_requirement){
       ith_var_MP_B[c] = NA_REAL;
     }
     if(ith_sample_measurements_B.size() >= replicate_number_requirement){
       ith_var_MP_B[c] = var(ith_sample_measurements_B);  
     }
     ++c;
   }
   
   
   int effective_N_A = N;
   int effective_n_A = n;
   int effective_N_B = N;
   int effective_n_B = n;
   double var_MP_A = 0;
   double var_MP_B = 0;
   
   for(int j = 0; j < n; ++j){
     bool is_na_var_A_j = ISNAN(ith_var_MP_A[j]);
     bool is_na_var_B_j = ISNAN(ith_var_MP_B[j]);
     if(!is_na_var_A_j){
       var_MP_A += ith_var_MP_A[j];
     }
     else if(is_na_var_A_j){
       effective_n_A = effective_n_A - 1;
     }
     if(!is_na_var_B_j){
       var_MP_B += ith_var_MP_B[j];
     }
     else if(is_na_var_B_j){
       effective_n_B = effective_n_B - 1;
     }
   }
   
   if(effective_n_A >= 1){
     var_MP_A = var_MP_A / effective_n_A; 
   }
   if(effective_n_A < 1){
     var_MP_A = NA_REAL;
   }
   if(effective_n_B >= 1){
     var_MP_B = var_MP_B / effective_n_B;
   }
   if(effective_n_B < 1){
     var_MP_B = NA_REAL;
   }
   double lambda = 0;
   bool is_na_pooled_var_A = ISNAN(var_MP_A);
   bool is_na_pooled_var_B = ISNAN(var_MP_B);
   if(is_na_pooled_var_A | is_na_pooled_var_B){
     List out = List::create(Named("zeta") = NA_REAL);
     return out;
   }
   else if((!is_na_pooled_var_A) & (!is_na_pooled_var_B)){
     lambda = var_MP_A / var_MP_B;  
   }
   double mean_MP_A = 0;
   double mean_MP_B = 0;
   
   for(int i = 0; i < N; ++i){
     bool is_na_mean_A = ISNAN(MP_A[i]);
     bool is_na_mean_B = ISNAN(MP_B[i]);
     if(!is_na_mean_A){
       mean_MP_A += MP_A[i];
     }
     else if(is_na_mean_A){
       effective_N_A = effective_N_A - 1;
     }
     if(!is_na_mean_B){
       mean_MP_B += MP_B[i];
     }
     else if(is_na_mean_B){
       effective_N_B = effective_N_B - 1;
     }
   }
   if(effective_N_A >= 1){
     mean_MP_A = mean_MP_A / effective_N_A;  
   }
   if(effective_N_A < 1){
     mean_MP_A = NA_REAL;
   }
   if(effective_N_B >= 1){
     mean_MP_B = mean_MP_B / effective_N_B;
   }
   if(effective_N_B < 1){
     mean_MP_B = NA_REAL;  
   }
   
   if(silence == 0 or silence == -1){
     Rcout << "----- Progress (1 / 3) Calculations of variances and means -------" << "\n";
     Rcout << "## Pooled variance of MP_A is : " << var_MP_A << "\n";
     Rcout << "## Pooled variance of MP_B is : " << var_MP_B << "\n";
     Rcout << "## Mean of MP_A is : " << mean_MP_A << "\n";
     Rcout << "## Mean of MP_B is : " << mean_MP_B << "\n";
     Rcout << "## lambda is : " << lambda << "\n";
     Rcout << "## N is : " << N << "\n";
     Rcout << "------------------------------------------------------------------" << "\n";  
   }
   
   NumericVector x = MP_B;
   NumericVector y = MP_A;
   double mx = mean_MP_B;
   double my = mean_MP_A;
   
   double msxx = 0;
   double msyy = 0;
   double msxy = 0;
   
   for(int i = 0; i < N; ++i){
     bool na_check_x = ISNAN(x[i]);
     bool na_check_y = ISNAN(y[i]);
     if((!na_check_x) & (!na_check_y)){
       msxx = msxx + pow(x[i] - mx, 2);
       msyy = msyy + pow(y[i] - my, 2);
       msxy = msxy + (x[i] - mx) * (y[i] - my);
     }
   }
   
   msxx = msxx / (effective_N_B - 1.0);
   msyy = msyy / (effective_N_A - 1.0);
   msxy = msxy / (effective_N_A - 1.0);
   
   double sub_expression_1 = msyy - lambda * msxx;
   double sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
   double sub_expression_3 = 2 * msxy;
   double b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
   double varb1 = (pow(b1, 2) / ((effective_N_A) * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
   double hvar = (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda);
   double varpar = varb1 * msxx + varb1 * hvar + (1 + 1.0 / effective_N_A) * (pow(b1, 2) + lambda) * hvar;
   double zeta = varpar / (var_MP_A + var_MP_B * pow(b1, 2));
   
   if(silence == 0 or silence == -1){
     Rcout << "----- Progress (2 / 3) Calculations of variances and means -------" << "\n";
     Rcout << "## var[y0 - y0^] : " << varpar << "\n";
     Rcout << "## b1 : " << b1 << "\n";
     Rcout << "## Pooled variance of MP_B is : " << var_MP_B << "\n";
     Rcout << "## Estimated variance of MP_B is : " << hvar << "\n";
     Rcout << "------------------------------------------------------------------" << "\n";  
   }
   
   List out = List::create(Named("zeta") = zeta);
   return out;
 }



