#include <Rcpp.h>
using namespace Rcpp;

//' Obtain the residuals based on data and imprecision estimates
//' 
//' @title Obtain the residuals based on data and imprecision estimates
//' @name residuals_eqa
//' 
//' @param data \code{list} or \code{data table} - Mean-of-replicates clinical sample data with elements/columns \code{SampleID}, \code{MP_A} and \code{MP_B}
//' @param imprecision_estimates \code{list} or \code{data table} - Imprecision estimates e.g., that which are outputted by \code{global_precision_estimates()}. Minimum requirement: \code{lambda} and \code{Var_B} must be part of the  
//' @param method character - Which method should be used to obtain the residuals. Default is \code{fg}. At the moment, possible values for method is
//' \itemize{
//'   \item{\code{fg: }}{Standard Deming regression, using components from W. Fuller and J. Gillards work}
//'   \item{\code{clsi: }}{Standard Deming regression based the EP14 standard, which is derived by J. Vaks}
//'   \item{\code{ols: }}{Ordinary least squares regression}
//' }
//' @param studentize integer - \code{1} (default) for yes, and \code{0} for no.
//'
//' @description A simple function extracting the residuals based on \code{data}
//' 
//'
//' @details The practical differences of outcome between the three available methods may be small. \code{fg} estimates latent values, and uses n - 1 degrees of freedom instead of n which is used by CLSI. CLSI does not estimate latent variables. OLS does of course ignore imprecision in x. 
//'
//' @return \code{list} containing two numeric vectors, that are residuals and fitted values
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   parameters <- list(n = 20, R = 3, cvx = 0.03, cvy = 0.02, cil = 10, ciu = 15)
//'   data <- simulate_eqa_data(training_parameters)
//'   data$SampleID <- as.character(data$SampleID)
//'   data$ReplicateID <- as.character(data$ReplicateID)
//'   imprecision <- global_precision_estimates(data)
//'   mean_of_replicates_data <- fun_of_replicates(data)
//'   # Extracting raw residuals based on the clsi Deming regression approach
//'   residuals_eqa(mean_of_replicates_data, imprecision, "clsi", 0)
//' }

// [[Rcpp::export]]
List residuals_eqa(List data, List imprecision_estimates, String method = "fg", int studentize = 1) {
  // Needed global imprecision estimates
  float lambda = imprecision_estimates["lambda"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  int n_A = MP_A.size();
  int n_B = MP_B.size();
  if(n_A != n_B){
    stop("Lengths for MP_A and MP_B are not equal");
  }
  if(method == "fg"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x);
    float msyy = var(y);
    float msxy = 0;
    for(int i = 0; i < n_B; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy = msxy / (n_B - 1);
    
    float sub_expression_1 = msyy - lambda * msxx;
    float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    float sub_expression_3 = 2 * msxy;
    float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    float b0 = my - b1 * mx;
    NumericVector yhat(n_A);
    for(int i = 0; i < n_A; ++i){
      float tau_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      yhat[i] = b0 + b1 * tau_i;
    }
    NumericVector r(n_A);
    for(int i = 0; i < n_A; ++i){
      r[i] = y[i] - yhat[i];
    }
    if(studentize == 1){
      float mr = mean(r);
      float sr = sd(r);
      for(int i = 0; i < n_A; ++i){
        r[i] = (r[i] - mr) / sr;
      }
    }
    List out = List::create(Named("residuals") = r, Named("fitted") = yhat);
    return out;
  }
  
  else if(method == "clsi"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x) * (n_B - 1) / n_B;
    float msyy = var(y) * (n_A - 1) / n_A;
    float msxy = 0;
    for(int i = 0; i < n_B; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy = msxy / n_B;
    
    float sub_expression_1 = msyy - lambda * msxx;
    float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    float sub_expression_3 = 2 * msxy;
    float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    float b0 = my - b1 * mx;
    NumericVector yhat(n_A);
    for(int i = 0; i < n_A; ++i){
      yhat[i] = b0 + b1 * x[i];
    }
    NumericVector r(n_A);
    for(int i = 0; i < n_A; ++i){
      r[i] = y[i] - yhat[i];
    }
    if(studentize == 1){
      float mr = mean(r);
      float sr = sd(r);
      for(int i = 0; i < n_A; ++i){
        r[i] = (r[i] - mr) / sr;
      }
    }
    List out = List::create(Named("residuals") = r, Named("fitted") = yhat);
    return out;
  }
  
  else if(method == "ols"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x);
    float msxy = 0;
    for(int i = 0; i < n_B; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy = msxy / (n_B - 1);
    
    float b1 = msxy / msxx;
    float b0 = my - b1 * mx;
    NumericVector yhat(n_A);
    for(int i = 0; i < n_A; ++i){
      yhat[i] = b0 + b1 * x[i];
    }
    NumericVector r(n_A);
    for(int i = 0; i < n_A; ++i){
      r[i] = y[i] - yhat[i];
    }
    if(studentize == 1){
      float mr = mean(r);
      float sr = sd(r);
      for(int i = 0; i < n_A; ++i){
        r[i] = (r[i] - mr) / sr;
      }
    }
    List out = List::create(Named("residuals") = r, Named("fitted") = yhat);
    return out;
  }
  
  NumericVector x = MP_B;
  NumericVector y = MP_A;
  float mx = mean(x);
  float my = mean(y);
  float msxx = var(x);
  float msxy = 0;
  for(int i = 0; i < n_B; ++i){
    msxy = msxy + (x[i] - mx) * (y[i] - my);
  }
  msxy = msxy / (n_B - 1);
  
  float b1 = msxy / msxx;
  float b0 = my - b1 * mx;
  NumericVector yhat(n_A);
  for(int i = 0; i < n_A; ++i){
    yhat[i] = b0 + b1 * x[i];
  }
  NumericVector r(n_A);
  for(int i = 0; i < n_A; ++i){
    r[i] = y[i] - yhat[i];
  }
  if(studentize == 1){
    float mr = mean(r);
    float sr = sd(r);
    for(int i = 0; i < n_A; ++i){
      r[i] = (r[i] - mr) / sr;
    }
  }
  List out = List::create(Named("residuals") = r, Named("fitted") = yhat);
  return out;
}


