#include <Rcpp.h>
using namespace Rcpp;

//' @title Calculate Model Residuals for External Quality Assessment Data
//' @name residuals_eqa
//' 
//' @param data A list or data.table containing:
//'        \itemize{
//'           \item \code{SampleID}: A \code{character} vector of sample identifiers
//'           \item \code{MP_A}: A \code{numeric} vector of measurements from method A (response)
//'           \item \code{MP_B}: A \code{numeric} vector of measurements from method B (predictor)
//'        }
//' @param imprecision_estimates A \code{list} or \code{data.table} containing
//'        \itemize{
//'           \item \code{lambda}: A positive numeric value representing the repeatability
//'                                variance ratio between response and predictor.
//'        }
//'        Can be calculated using the \code{global_precision_estimates()} function.
//' @param method A \code{character} string specifying the linear regression model
//'        applied for obtaining the residuals and fitted values.  Options are:
//'        \itemize{
//'           \item \code{'fg'} (default): Deming regression, incorporating results from W. Fuller and J. Gillard.
//'           \item \code{'clsi'}: Deming regression based on the CLSI EP14 standard.
//'           \item \code{'ols'}: Ordinary Least Squares regression.
//'        }
//' @param studentize A non-missing \code{logical} value. If \code{TRUE} (default),
//'        the standardized residuals are computed. Here, this means that the sample mean
//'        of the residuals is subtracted from each residual and this result is divided by
//'        the sample standard deviation of the residuals. Raw residuals are returned
//'        if this is set to \code{FALSE}. Standardized residuals can be
//'        beneficial if one requires to compare the residuals' distribution to the standard
//'        normal distribution or a particular student t-distribution.
//' 
//' @param unit_sd A non-missing \code{logical} value. If \code{TRUE}, the standard
//'        deviation of the residuals are assumed to be equal to \code{1}. 
//'        
//' @param invalid_NA A non-missing \code{logical} value. If \code{TRUE},
//'        a \code{list} is \code{NA}-values is returned instead of throwing
//'        an error.
//'
//' @description
//' This function computes residuals and fitted values for external quality assessment (EQA)
//' clinical sample data using various regression methods. It is designed to work with
//' paired measurement results from two in vitro diagnostic medical devices (IVD-MDs).
//'
//' @details
//' The function handles missing values by removing them before calculations. It supports
//' three regression methods, each with slightly different approaches to calculating residuals
//' and fitted values. The \code{'fg'} method is generally recommended.
//' While the outcome differences between the three avaiable methods may be small,
//' they can still be noteworthy in some cases. 
//'
//' @return 
//' This function returns a \code{list} comprising two \code{numeric} vectors named
//' \code{residuals} and \code{fitted}. These are the residuals and fitted values, respectively,
//' based on the input \code{data}, \code{imprecision_estimates}, and \code{method}.
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
//' 
//' @seealso
//' \code{\link{global_precision_estimates}} for calculating imprecision estimates
//' 
//' @references
//' Fuller, W. A. (1987). Measurement Error Models. John Wiley & Sons.
//' 
//' Gillard J. W. (2010). An Overview of Linear Structural Models in Errors in
//' Variables Regression. REVSTAT-Statistical Journal 8(1), 57-80
//' 

// [[Rcpp::export]]
List residuals_eqa(const List& data,
                   const List& imprecision_estimates,
                   const String& method = "fg",
                   bool studentize = true,
                   bool unit_sd = false,
                   bool invalid_NA = false) {
  
  // Input validation
  if (!data.containsElementNamed("SampleID") || 
      !data.containsElementNamed("MP_A") || 
      !data.containsElementNamed("MP_B")) {
      if (!invalid_NA) {
        stop("'data' must contain 'SampleID', 'MP_A', and 'MP_B'");
      }
      return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
  }
  
  // Checks if lambda exists
  if (!imprecision_estimates.containsElementNamed("lambda")) {
    if (!invalid_NA) {
      stop("'imprecision_estimates' must contain 'lambda'");
    }
    return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
  }
  
  // Extract lambda, now sure that it exists
  double lambda = as<double>(imprecision_estimates["lambda"]);
  
  if (ISNAN(lambda) || lambda < 0) {
    if (!invalid_NA) {
      stop("'lambda' must be a positive double");
    }
    return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
  }
  
  //NumericVector MP_A = data["MP_A"];
  //NumericVector MP_B = data["MP_B"];
  
  //if (MP_A.size() != MP_B.size()) {
  //  if(!invalid_NA){
  //    stop("'MP_A' and 'MP_B' must have the same length");  
  //  }
  //  return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
  //}
  
  // Remove NA values
  //LogicalVector valid = !is_na(MP_A) & !is_na(MP_B);
  //MP_A = MP_A[valid];
  //MP_B = MP_B[valid];
  
  NumericVector raw_MP_A = data["MP_A"];
  NumericVector raw_MP_B = data["MP_B"];
  
  if(raw_MP_A.size() != raw_MP_B.size()){
    if(!invalid_NA){
      stop("Lengths for 'MP_A'[%i] and 'MP_B'[%i] found in 'data' are not equal",
           raw_MP_A.size(), raw_MP_B.size());  
    }
    return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
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
  
  
  int n = MP_A.size();
  double n_double = static_cast<double>(n);
  
  NumericVector fitted(n);
  NumericVector residuals(n);
  
  if (method == "fg" || method == "clsi") {
    double mx = mean(MP_B);
    double my = mean(MP_A);
    double msxx = var(MP_B);
    double msyy = var(MP_A);
    double msxy = sum((MP_B - mx) * (MP_A - my)) / (n_double - 1.0);
    
    //double msxy = 0;
    //for(int i = 0; i < n; ++i){
    //  msxx += pow(x[i] - mx, 2) / n_double;
    //  msyy += pow(y[i] - my, 2) / n_double;
    //  msxy += (x[i] - mx) * (y[i] - my) / n_double;
    //}
    
    double sub_expression_1 = msyy - lambda * msxx;
    double sub_expression_2 = std::sqrt(std::pow(sub_expression_1, 2) + 4 * lambda * std::pow(msxy, 2));
    double b1 = (sub_expression_1 + sub_expression_2) / (2.0 * msxy);
    double b0 = my - b1 * mx;
    
    if (method == "fg") {
      for (int i = 0; i < n; ++i) {
        double tau_i = ((lambda / (lambda + std::pow(b1, 2))) * MP_B[i]) + 
          (b1 / (lambda + std::pow(b1, 2))) * (MP_A[i] - b0);
        fitted[i] = b0 + b1 * tau_i;
        residuals[i] = MP_A[i] - fitted[i];
      }
    } else { // method == "clsi"
      for (int i = 0; i < n; ++i) {
        fitted[i] = b0 + b1 * MP_B[i];
        residuals[i] = MP_A[i] - fitted[i];
      }
    }
  } else if (method == "ols") {
    double mx = mean(MP_B);
    double my = mean(MP_A);
    double msxx = var(MP_B);
    double msxy = sum((MP_B - mx) * (MP_A - my)) / (n_double - 1.0);
    
    //double msxy = 0;
    //for(int i = 0; i < n; ++i){
    //  msxx += pow(x[i] - mx, 2) / n_double;
    //  msyy += pow(y[i] - my, 2) / n_double;
    //  msxy += (x[i] - mx) * (y[i] - my) / n_double;
    //}
    
    double b1 = msxy / msxx;
    double b0 = my - b1 * mx;
    
    for (int i = 0; i < n; ++i) {
      fitted[i] = b0 + b1 * MP_B[i];
      residuals[i] = MP_A[i] - fitted[i];
    }
  } else {
    if (!invalid_NA) {
      stop("Invalid 'method'. Use 'fg', 'clsi', or 'ols'");
    }
    return List::create(Named("residuals") = NA_REAL, Named("fitted") = NA_REAL);
  }
  
  if (studentize) {
    double mean_residuals = mean(residuals);
    double sd_residuals = sd(residuals);
    if(unit_sd){
      sd_residuals = 1.0;
    }
    residuals = (residuals - mean_residuals) / sd_residuals;
  }
  
  return List::create(Named("residuals") = residuals, Named("fitted") = fitted);
}

