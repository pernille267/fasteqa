#include <Rcpp.h>
using namespace Rcpp;

//' Obtain the residuals based on data and imprecision estimates
//' 
//' @title Obtain the residuals based on data and imprecision estimates
//' @name residuals_eqa
//' 
//' @param data A \code{list} or \code{data.table} that includes mean values of replicated clinical sample data, grouped by \code{SampleID}. The data for the two compared IVD-MDs under consideration should be found in the columns named \code{MP_A} and \code{MP_B}. The \code{SampleID} should be a character vector, whereas \code{MP_A} and \code{MP_B} must be numeric vectors.
//' @param imprecision_estimates A \code{list} or \code{data.table} carrying imprecision estimates for IVD-MD, ideally matching those provided by the \code{global_precision_estimates()} function. A non-missing positive value for \code{lambda} is the minimum requirement.
//' @param method A \code{character} string specifying the linear regression model applied for obtaining the residuals and fitted values. The default value is \code{fg}. The possible values for \code{method} include:
//' \itemize{
//'   \item{\code{fg: }}{Refers to Deming regression, incorporating elements from the work of W. Fuller and J. Gillard.}
//'   \item{\code{clsi: }}{Refers to Deming regression, based on the CLSI EP14 standard.}
//'   \item{\code{ols: }}{Represents ordinary least squares regression.}
//' }
//' @param studentize An \code{integer}. When set to \code{1} (default), the residuals get studentized, implying the mean of the residuals is subtracted from each residual and then divided by the standard deviation of the residuals. If set to \code{0}, residuals are not studentized. Studentized residuals can be beneficial if there is a need to compare the residuals' distribution to the standard normal distribution or a specific t-distribution.
//'
//' @description This function efficiently extracts the residuals and fitted values based on the provided \code{data}, \code{imprecision_estimates}, and \code{method}.
//'
//' @details While the outcome differences between the three provided methods may be minute, they can still be noteworthy. Suppose 'n' is the count of unique clinical samples. The \code{fg} method estimates latent concentration values and utilizes n - 1 degrees of freedom. In contrast, \code{clsi} does not estimate latent variables and uses 'n' degrees of freedom. The model associated with \code{ols} overlooks the imprecision in \code{MP_B}.
//'
//' @return This function returns a \code{list} comprising two numeric vectors named \code{residuals} and \code{fitted}. These are the residuals and fitted values, respectively, based on the input \code{data}, \code{imprecision_estimates}, and \code{method}.
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
List residuals_eqa(List data, List imprecision_estimates, String method = "fg", int studentize = 1, bool invalid_NA = false) {
  
  if (!(data.containsElementNamed("SampleID") && data.containsElementNamed("MP_A") && data.containsElementNamed("MP_B"))) {
    if(!invalid_NA){
      stop("'data' must contain all the following elements: 'SampleID', 'MP_A', and 'MP_B'. One or more of these elements are missing in the current 'data'.");  
    }
    else{
      return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);  
    }
  }
  
  if (!(imprecision_estimates.containsElementNamed("lambda"))) {
    if(!invalid_NA){
      stop("'imprecision_estimates' must contain 'lambda', but is currently missing in 'imprecision_estimates'.");  
    }
    else{
      return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);  
    }
  }
  
  // Needed global imprecision estimates
  float lambda = imprecision_estimates["lambda"];
  
  // Checks if 'lambda' is an NA value
  bool lambda_na = ISNAN(lambda);
 
  // Throws an error if 'lambda' is an NA value
  if(lambda_na){
    if(!invalid_NA){
      stop("The provided value of 'lambda' is missing. 'lambda' is expected to be a positive float-point value.");
    }
    else{
      return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
    }
  }
  
  // Checks if 'lambda' is a negative value and throws an error if it is
  bool lambda_negative = lambda < 0;
  if(lambda_negative){
    if(!invalid_NA){
      stop("The provided value of 'lambda'[%d] is negative. 'lambda' is expected to be a positive floating-point value.",
           lambda);
    }
    else{
      return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
    }
  }
  
  NumericVector raw_MP_A = data["MP_A"];
  NumericVector raw_MP_B = data["MP_B"];
  
  if(raw_MP_A.size() != raw_MP_B.size()){
    stop("Lengths for 'MP_A'[%i] and 'MP_B'[%i] found in 'data' are inequal",
         raw_MP_A.size(), raw_MP_B.size());
  }
  
  NumericVector MP_A;
  NumericVector MP_B;
  
  // Remove NA-values if present before calculations of residuals and fitted values
  for(int i = 0; i < raw_MP_A.size(); ++i){
    bool is_na_x = ISNAN(raw_MP_B[i]);
    bool is_na_y = ISNAN(raw_MP_A[i]);
    if(is_na_x || is_na_y){
      continue;
    }
    MP_B.push_back(raw_MP_B[i]);
    MP_A.push_back(raw_MP_A[i]);
  }
  
  int n = MP_B.size();
  float n_float = static_cast<float>(n);
  
  if(method == "fg"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x);
    float msyy = var(y);
    float msxy = 0;
    for(int i = 0; i < n; ++i){
      msxy += (x[i] - mx) * (y[i] - my) / (n_float - 1.0);
    }
    
    float sub_expression_1 = msyy - lambda * msxx;
    float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    float sub_expression_3 = 2 * msxy;
    float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    float b0 = my - b1 * mx;
    
    NumericVector fitted(n);
    NumericVector residuals(n);
    for(int i = 0; i < n; ++i){
      float tau_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      fitted[i] = b0 + b1 * tau_i;
      residuals[i] = y[i] - fitted[i];
    }
    if(studentize == 1){
      float mean_residuals = mean(residuals);
      float sd_residuals = sd(residuals);
      for(int i = 0; i < n; ++i){
        residuals[i] = (residuals[i] - mean_residuals) / sd_residuals;
      }
    }
    List out = List::create(Named("residuals") = residuals, Named("fitted") = fitted);
    return out;
  }
  
  else if(method == "clsi"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    
    float mx = mean(x);
    float my = mean(y);
    float msxx = 0;
    float msyy = 0;
    float msxy = 0;
    for(int i = 0; i < n; ++i){
      msxx += pow(x[i] - mx, 2) / n_float;
      msyy += pow(y[i] - my, 2) / n_float;
      msxy += (x[i] - mx) * (y[i] - my) / n_float;
    }
    
    float sub_expression_1 = msyy - lambda * msxx;
    float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    float sub_expression_3 = 2 * msxy;
    float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    float b0 = my - b1 * mx;
    
    NumericVector fitted(n);
    NumericVector residuals(n);
    for(int i = 0; i < n; ++i){
      fitted[i] = b0 + b1 * x[i];
      residuals[i] = y[i] - fitted[i];
    }
    
    if(studentize == 1){
      float mean_residuals = mean(residuals);
      float sd_residuals = sd(residuals);
      for(int i = 0; i < n; ++i){
        residuals[i] = (residuals[i] - mean_residuals) / sd_residuals;
      }
    }
    
    List out = List::create(Named("residuals") = residuals, Named("fitted") = fitted);
    return out;
  }
  
  else if(method == "ols"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x);
    float msxy = 0;
    for(int i = 0; i < n; ++i){
      msxy += (x[i] - mx) * (y[i] - my) / (n_float - 1.0);
    }
    
    float b1 = msxy / msxx;
    float b0 = my - b1 * mx;
    NumericVector fitted(n);
    NumericVector residuals(n);
    for(int i = 0; i < n; ++i){
      fitted[i] = b0 + b1 * x[i];
      residuals[i] = y[i] - fitted[i];
    }
    
    if(studentize == 1){
      float mean_residuals = mean(residuals);
      float sd_residuals = sd(residuals);
      for(int i = 0; i < n; ++i){
        residuals[i] = (residuals[i] - mean_residuals) / sd_residuals;
      }
    }
    List out = List::create(Named("residuals") = residuals, Named("fitted") = fitted);
    return out;
  }
  
  if(!invalid_NA){
    std::string input_method = method.get_cstring();
    stop("The provided input for 'method'[%s] is invalid. Valid inputs for 'method' include 'fg', 'clsi' and 'ols'.",
         input_method.c_str());
  }
  
  List out = List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
  return out;
}


