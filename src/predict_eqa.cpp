#include <Rcpp.h>

using namespace Rcpp;

//' @title Prediction Interval Estimation for EQA Data via Deming or
//'        OLS Regression
//' @name predict_eqa
//'
//' @description
//' This function estimates prediction intervals for evaluated material data by
//' applying either Deming (\code{method = 'fg'} or \code{method = 'clsi'}) or
//' Ordinary Least Squares (\code{method = 'ols'}) regression methodologies.
//'
//' @param data A \code{list} or \code{data.table}. Must contain:
//'             \itemize{
//'                 \item \code{SampleID: } A \code{character} vector. The
//'                 clinical sample identifiers.
//'                 \item \code{MP_A: } A \code{numeric} vector. The means of
//'                 replicated measurements from IVD-MD \code{MP_A} (response).
//'                 \item \code{MP_B: } A \code{numeric} vector. The means of
//'                 replicated measurements from IVD-MD \code{MP_B} (predictor).
//'             }
//' @param new_data A \code{list} or \code{data.table}. Can contain:
//'                 \itemize{
//'                     \item \code{SampleID: } A \code{character} vector. The
//'                     evaluated material sample identifiers. Optional.
//'                     \item \code{MP_A: } A \code{numeric} vector. The means
//'                     of replicated measurements from IVD-MD \code{MP_A}.
//'                     Optional.
//'                     \item \code{MP_B: } A \code{numeric} vector. The means
//'                     of replicated measurements from IVD-MD \code{MP_B}.
//'                     Mandatory.
//'                 }
//'                 The structure of \code{new_data} determines the output
//'                 type. For constructing prediction band (PB) data, include
//'                 only \code{MP_B}. For constructing commutability
//'                 evaluation data, include all three. Inside checks can only
//'                 be performed if \code{MP_A} is included.
//'                
//' @param imprecision_estimates A \code{list} or \code{data.table}. Must
//'                              include:
//'                              \itemize{
//'                                 \item \code{Var_B: } Pooled variance 
//'                                       estimate of \eqn{\sigma_h^2}.
//'                                 \item \code{lambda: } Estimate of
//'                                       \eqn{\sigma_v^2 / \sigma_h^2}. 
//'                              }
//'                              See details.
//' @param R An \code{integer}. The average number of replicates on which
//'          \code{new_data} is based. \code{count_samplewise_replicates()} can
//'          be employed to generate a suitable input for this parameter
//'          directly.
//' @param R_ratio A \code{double}. The ratio of the number of replicates used
//'                in \code{new_data} and \code{data}. Only relevant if
//'                \code{method = 'fg'} or \code{method = 'ols'}. Defaults to
//'                \code{1L}.
//' @param method A \code{character} string. The desired method for estimating
//'               prediction intervals. Possible prediction estimation
//'               approaches includes:
//'               \itemize{
//'                   \item \code{fg: } Implements standard Deming regression
//'                         with prediction intervals calculated using the
//'                         'Fuller & Gillard' approach. See details.
//'                   \item \code{clsi: } Utilizes standard Deming regression,
//'                         with prediction intervals derived from the CLSI
//'                         method. See details.
//'                   \item \code{ols: } Implements Ordinary Least Squares
//'                         regression.
//'               }
//' @param level A \code{double}. Must be between \code{0} and \code{1}. The
//'              nominal confidence level for the estimated prediction
//'              intervals. Defaults to \code{0.99} (\eqn{99\%}).
//' @param allow_reverse_regression A \code{logical} value. If \code{TRUE} and
//'                                 \code{method = 'ols'}, the response and
//'                                 predictor change roles if
//'                                 \code{lambda < 0.5}. Defaults to
//'                                 \code{FALSE}.
//' @param rounding An \code{integer}. The desired number of decimal places for
//'                 the predictions and prediction intervals. Defaults to
//'                 \code{3L}, offering sufficient precision in relevant
//'                 applications.
//'
//' @details
//' For commutability evaluation purposes, inclusion of all \code{SampleID},
//' \code{MP_A}, and \code{MP_B} in \code{new_data} is required. 
//' If only prediction band data is required, inclusion of \code{MP_B} in
//' \code{new_data} is the sole requirement.
//' 
//' Imprecision Estimates
//' 
//' Estimating prediction intervals may require certain imprecision estimates,
//' which can be estimated from the raw data. Which imprecision estimates that
//' are required depends on which method that is used to construct the
//' prediction intervals. For \code{method = 'fg'} or \code{method = 'ols'},
//' \code{lambda} is required. \code{lambda} is calculated using
//' 
//' \eqn{\lambda = \frac{\hat{\sigma}_v^2}{\hat{\sigma}_h^2}},
//' 
//' where \eqn{\hat{\sigma}_v^2} and \eqn{\hat{\sigma}_h^2} are pooled 
//' variance estimates. If differences in non-selectivity (DINS) or equation
//' error cannot be ruled out, \eqn{\lambda} will generally underestimate the
//' true parameter \eqn{\Lambda}. Keep this in mind!
//' 
//' For \code{method = 'clsi'}, both \code{Var_B} and \code{lambda} are
//' required. \code{lambda} is calculated as before, and \code{Var_B} is just
//' \eqn{\hat{\sigma}_h^2}.
//' 
//' Note that both \code{lambda} and \code{Var_B} always must be given as
//' input. However, if one of them is not necessary to estimate prediction
//' intervals, it may take an arbitrary value. It is nevertheless recommended
//' to always include actual values \code{lambda} and \code{Var_B},
//' independently of whether they are actual used or not.
//' 
//' Note also that \code{global_precision_estimates()} can be used to calculate
//' both \code{lambda} and \code{Var_B}.
//' 
//' @return 
//' A \code{list}. The resulting estimated prediction interval data based on
//' the function inputs. Contain the following elements:
//' \itemize{
//'     \item \code{SampleID: } A \code{character} vector. The ID(s) of the
//'           evaluated material(s).
//'     \item \code{MP_B: } A \code{numeric} vector. The measurement result(s)
//'           of IVD-MD \code{MP_B} (predictor).
//'     \item \code{MP_A: } A \code{numeric} vector. The measurement result(s)
//'           of IVD-MD \code{MP_A} (response).
//'     \item \code{prediction: } A \code{numeric} vector. The predicted
//'           value(s) of \code{MP_A} given the value(s) of \code{MP_B}.
//'     \item \code{lwr: } A \code{numeric} vector. The lower limit(s) of the
//'           estimated prediction interval(s) of \code{MP_A} given the
//'           value(s) of \code{MP_B}.
//'     \item \code{upr: } A \code{numeric} vector. The upper limit(s) of the
//'            estimated prediction interval(s) of \code{MP_A} given the
//'            value(s) of \code{MP_B}.
//'     \item \code{inside: } An \code{integer} vector. The Inside checks.
//'           \code{1} if \code{MP_A} is inside the estimated prediction
//'           interval for \code{MP_A}. Otherwise, \code{0}.
//' }
//' Note: the output \code{list} may be converted to \code{data.table} by using
//' the \code{setDT()} function from the \code{data.table} package in R.
//'
//' @examples
//' # Required packages
//' library(fasteqa)
//' library(data.table)
//' 
//' # Read data and convert to data.table
//' test_data_example <- as.data.table(test_data)
//' 
//' # Log-transform data
//' test_data_example[, MP_A := log(MP_A)]
//' test_data_example[, MP_B := log(MP_B)]
//' 
//' # Use one of the clinical samples as a fictive evaluated material sample
//' test_cs_data <- test_data_example[SampleID != "1"]
//' test_eq_data <- test_data_example[SampleID == "1"]
//' 
//' # Estimate repeatability uncertainty statistics
//' impr_data <- global_precision_estimates(test_cs_data)
//' 
//' # Calculate mean-of-replicates data
//' test_cs_data <- test_cs_data[, fun_of_replicates(.SD)]
//' test_eq_data <- test_eq_data[, fun_of_replicates(.SD)]
//' 
//' # Calculate 95% OLS commutability evaluation data
//' ols_pi <- predict_eqa(data = test_cs_data,
//'                       new_data = test_eq_data,
//'                       imprecision_estimates = impr_data,
//'                       method = "ols",
//'                       level = 0.95,
//'                       allow_reverse_regression = FALSE,
//'                       rounding = 3L)
//' 
//' # Calculate 95% F-G Deming commutability evaluation data
//' fg_pi <- predict_eqa(data = test_cs_data,
//'                      new_data = test_eq_data,
//'                      imprecision_estimates = impr_data,
//'                      method = "fg",
//'                      level = 0.95,
//'                      allow_reverse_regression = FALSE,
//'                      rounding = 3L)
//' 
//' # Calculate 95% CLSI EP14 Deming commutability evaluation data
//' clsi_pi <- predict_eqa(data = test_cs_data,
//'                       new_data = test_eq_data,
//'                       imprecision_estimates = impr_data,
//'                       method = "clsi",
//'                       level = 0.95,
//'                       allow_reverse_regression = FALSE,
//'                       rounding = 3L)
//' 
//' # Convert to data.table objects
//' lapply(X = list(ols_pi,
//'                 fg_pi,
//'                 clsi_pi),
//'        FUN = setDT)
//' 
//' # Gather into one data.table
//' pis <- rbindlist(list("ols" = ols_pi,
//'                       "fg" = fg_pi,
//'                       "clsi" = clsi_pi),
//'                   idcol = "method")
//' 
//' # The result                  
//' print(pis)
//' 


// [[Rcpp::export]]
List predict_eqa(List data, List new_data, List imprecision_estimates,
                 int R = 3, double R_ratio = 1, String method = "fg",
                 double level = 0.99, bool allow_reverse_regression = true, int rounding = 3) {
  
  // Extract required information from imprecision_estimates
  double lambda = imprecision_estimates["lambda"];
  double Var_B = imprecision_estimates["Var_B"];
  
  // Checks if lambda and Var_B are valid. If not, use ols instead.
  if(ISNAN(lambda) || lambda <= 0.0 || ISNAN(Var_B) || Var_B <= 0.0){
    if(method == "fg"){
      method = "ols";  
    }
  }
  
  // Extract traning data
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  
  // Calculate length of traning data
  int n = MP_B.length();
  
  // Converting integer types to double types when needed in division!
  double n_double = static_cast<double>(n);
  double R_double = static_cast<double>(R);
  
  // Default values of existence checks
  int MP_A_exists = 0;
  int MP_B_exists = 0;
  int SampleID_exists = 0;
  
  // Valid names
  CharacterVector candidate_names = CharacterVector::create("MP_A", "MP_B", "SampleID");
  // Actual names
  CharacterVector col_names = new_data.names();
  // Number of columns in new_data
  int n_col_names = col_names.size();
  
  // Iterative checking whether valid columns of new_data is present
  for(int k = 0; k < n_col_names; ++k){
    if(candidate_names[0] == col_names[k]){
      MP_A_exists = 1;
    }
    else if(candidate_names[1] == col_names[k]){
      MP_B_exists = 1;
    }
    else if(candidate_names[2] == col_names[k]){
      SampleID_exists = 1;
    }
  }
  // m is the number of new observations used in prediction of MP_A | MP_B
  int m = 0;
  if(MP_B_exists == 1){
    NumericVector sizer = new_data["MP_B"];
    m = sizer.size();
  }
  else{
    stop("MP_B was not found in new_data. Inclusion of MP_B is mandatory!");
  }
  
  // Empty vectors that are to be filled
  NumericVector new_MP_A(m);
  NumericVector new_MP_B(m);
  CharacterVector new_SampleID(m);
  
  // Based on existence, we define the different columns of new_data
  if(MP_A_exists == 1){
    NumericVector candidate_MP_A = new_data["MP_A"];
    for(int j = 0; j < m; ++j){
      new_MP_A[j] = candidate_MP_A[j];
    }
  }
  if(MP_B_exists == 1){
    NumericVector candidate_MP_B = new_data["MP_B"];
    for(int j = 0; j < m; ++j){
      new_MP_B[j] = candidate_MP_B[j];
    }
  }
  if(SampleID_exists == 1){
    CharacterVector candidate_SampleID = new_data["SampleID"];
    for(int j = 0; j < m; ++j){
      new_SampleID[j] = candidate_SampleID[j];
    }
  }
  
  // Fuller-Gillard (F-G) approach
  if(method == "fg"){
    
    // Standard Statistics
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    double mx = mean(x);
    double my = mean(y);
    double msxx = var(x);
    double msyy = var(y);
    double msxy = 0;
    for(int i = 0; i < n; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy /= (n_double - 1.0);
    
    // Define vectors
    NumericVector latent(n);
    NumericVector nx = new_MP_B;
    NumericVector ny(m);
    NumericVector nl(m);
    
    // Fitting the standard Deming model
    double sub_expression_1 = msyy - lambda * msxx;
    double sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    double sub_expression_3 = 2 * msxy;
    double b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    double b0 = my - b1 * mx;
    
    // Dispersion
    double varb1 = (pow(b1, 2) / (n_double * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
    double hvar = (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda);
    double vvar = lambda * hvar;
    
    // Latent Variable Estimation
    double mu_hat = 0;
    for(int i = 0; i < n; ++i){
      double latent_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      mu_hat += latent_i / n_double;
      latent[i] = latent_i;
    }
    
    // Predicted (x, y)
    for(int j = 0; j < m; ++j){
      ny[j] = b0 + b1 * nx[j];
      nl[j] = mu_hat * (1 - (msxy / (b1 * msxx))) + ny[j] * (msxy / (b1 * msxx));
    }
    
    // Prediction Interval Estimation
    double t_quantile = R::qt((1 - level) / 2, n - 2, 0, 0);
    NumericVector var_pred_error(m);
    NumericVector lwr(m);
    NumericVector upr(m);
    IntegerVector inside(m);
    
    for(int j = 0; j < m; ++j){
      var_pred_error[j] = (1.0 + 1.0 / (n_double - 2.0)) * (varb1 * pow((nl[j] - mu_hat), 2) + (varb1 * hvar * R_ratio) + (1 + (1 / n_double)) * (pow(b1, 2) * hvar + vvar) * R_ratio);
      lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
      upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
    }
    
    // Rounding of Results
    lwr = round(lwr, rounding);
    upr = round(upr, rounding);
    ny = round(ny, rounding);
    
    // Output
    if(MP_A_exists == 1){
      NumericVector oy = new_MP_A;
      CharacterVector fill_SampleID(m);
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] && oy[j] <= upr[j]){
          inside[j] = 1;
          fill_SampleID[j] = NA_STRING;
        }
      }
      if(SampleID_exists == 1){
        List out = List::create(Named("SampleID") = new_SampleID,
                                Named("MP_B") = new_MP_B,
                                Named("MP_A") = new_MP_A,
                                Named("prediction") = ny,
                                Named("lwr") = lwr,
                                Named("upr") = upr,
                                Named("inside") = inside);  
        return out;
      }
      
      List out = List::create(Named("SampleID") = fill_SampleID,
                              Named("MP_B") = new_MP_B,
                              Named("MP_A") = new_MP_A,
                              Named("prediction") = ny,
                              Named("lwr") = lwr,
                              Named("upr") = upr,
                              Named("inside") = inside);
      
      return out;
    }
    
    else if(MP_A_exists == 0){
      CharacterVector fill_SampleID(m);
      NumericVector fill_MP_A(m);
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        fill_MP_A[j] = NA_REAL;
        fill_SampleID[j] = NA_STRING;
      }
      if(SampleID_exists == 1){
        List out = List::create(Named("SampleID") = new_SampleID,
                                Named("MP_B") = new_MP_B,
                                Named("MP_A") = fill_MP_A,
                                Named("prediction") = ny,
                                Named("lwr") = lwr,
                                Named("upr") = upr,
                                Named("inside") = inside);
        
        return out;
      }
      
      List out = List::create(Named("SampleID") = fill_SampleID,
                              Named("MP_B") = new_MP_B,
                              Named("MP_A") = fill_MP_A,
                              Named("prediction") = ny,
                              Named("lwr") = lwr,
                              Named("upr") = upr,
                              Named("inside") = inside);
      
      return out;
    }
    
    else{
      List out = List::create(Named("Output") = NA_STRING);
      return out;
    }
  }
  
  // CLSI, EP14 approach
  if(method == "clsi"){
    
    // Standard Statistics
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    double mx = mean(x);
    double my = mean(y);
    double msxx = var(x) * (n_double - 1) / n_double;
    double msyy = var(y) * (n_double - 1) / n_double;
    double msxy = 0;
    for(int i = 0; i < n; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy /= n_double;
    
    // Define vectors
    NumericVector nx = new_MP_B;
    NumericVector ny(m);
    
    // Fitting the standard Deming model
    double sub_expression_1 = msyy - lambda * msxx;
    double sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    double sub_expression_3 = 2 * msxy;
    double b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    double b0 = my - b1 * mx;
    
    // Dispersion
    double varb1 = (pow(b1, 2) / (n_double * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
    double hvar = Var_B;
    double vvar = lambda * hvar;
    
    // Predicted y
    for(int j = 0; j < m; ++j){
      ny[j] = b0 + b1 * nx[j];
    }
    
    // Prediction Interval Estimation
    double t_quantile = R::qt((1 - level) / 2, n * (R - 1), 0, 0);
    NumericVector var_pred_error(m);
    NumericVector lwr(m);
    NumericVector upr(m);
    IntegerVector inside(m);
    for(int j = 0; j < m; ++j){
      var_pred_error[j] = varb1 * pow((nx[j] - mx), 2) + (1 + (1 / n_double)) * (pow(b1, 2) * hvar + vvar) / R_double;
      lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
      upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
      
    }
    
    // Rounding of Results
    lwr = round(lwr, rounding);
    upr = round(upr, rounding);
    ny = round(ny, rounding);
    
    // Output
    if(MP_A_exists == 1){
      NumericVector oy = new_MP_A;
      CharacterVector fill_SampleID(m);
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] && oy[j] <= upr[j]){
          inside[j] = 1;
          fill_SampleID[j] = NA_STRING;
        }
      }
      if(SampleID_exists == 1){
        List out = List::create(Named("SampleID") = new_SampleID,
                                Named("MP_B") = new_MP_B,
                                Named("MP_A") = new_MP_A,
                                Named("prediction") = ny,
                                Named("lwr") = lwr,
                                Named("upr") = upr,
                                Named("inside") = inside);  
        return out;
      }
      
      List out = List::create(Named("SampleID") = fill_SampleID,
                              Named("MP_B") = new_MP_B,
                              Named("MP_A") = new_MP_A,
                              Named("prediction") = ny,
                              Named("lwr") = lwr,
                              Named("upr") = upr,
                              Named("inside") = inside);
      
      return out;
    }
    
    else if(MP_A_exists == 0){
      CharacterVector fill_SampleID(m);
      NumericVector fill_MP_A(m);
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        fill_MP_A[j] = NA_REAL;
        fill_SampleID[j] = NA_STRING;
      }
      if(SampleID_exists == 1){
        List out = List::create(Named("SampleID") = new_SampleID,
                                Named("MP_B") = new_MP_B,
                                Named("MP_A") = fill_MP_A,
                                Named("prediction") = ny,
                                Named("lwr") = lwr,
                                Named("upr") = upr,
                                Named("inside") = inside);
        
        return out;
      }
      
      List out = List::create(Named("SampleID") = fill_SampleID,
                              Named("MP_B") = new_MP_B,
                              Named("MP_A") = fill_MP_A,
                              Named("prediction") = ny,
                              Named("lwr") = lwr,
                              Named("upr") = upr,
                              Named("inside") = inside);
      
      return out;
    }
    
    else{
      List out = List::create(Named("Output") = NA_STRING);
      return out;
    }
  }
  if(method == "ols"){
    
    // Determine regression direction based on lambda
    const bool reverse_regression = (lambda < 1.0) & allow_reverse_regression;
    NumericVector x = reverse_regression ? MP_A : MP_B;
    NumericVector y = reverse_regression ? MP_B : MP_A;
    
    // Standard Statistics
    double sxx = var(x) * (n_double - 1);
    double mx = mean(x);
    double my = mean(y);
    double sxy = 0;
    double mse = 0;
    for(int i = 0; i < n; ++i){
      sxy = sxy + (x[i] - mx) * (y[i] - my);
    }
    
    // Define vectors
    NumericVector nx = new_MP_B;
    if((MP_A_exists == 1) & reverse_regression){
      nx = new_MP_A;
    }
    NumericVector ny(m);
    
    // Fitting the standard ordinary least squares model
    double b1 = sxy / sxx;
    double b0 = my - b1 * mx;
    
    // Dispersion
    for(int i = 0; i < n; ++i){
      mse = mse + pow(y[i] - b0 - b1 * x[i], 2);
    }
    mse = mse / (n_double - 2) * R_ratio;
    
    // Predicting y
    for(int j = 0; j < m; ++j){
      ny[j] = b0 + b1 * nx[j];
    }
    
    // Prediction Interval Calculation
    double t_quantile = R::qt((1 - level) / 2, n - 2, 0, 0);
    NumericVector var_pred_error(m);
    NumericVector lwr(m);
    NumericVector upr(m);
    IntegerVector inside(m);
    
    for(int j = 0; j < m; ++j){
      var_pred_error[j] = mse * (1 + (1 / n_double) + pow(nx[j] - mx, 2) / sxx);
      lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
      upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
    }
    
    // Rounding of Results
    lwr = round(lwr, rounding);
    upr = round(upr, rounding);
    ny = round(ny, rounding);
    
    // Output
    if(MP_A_exists == 1){
      NumericVector ox = reverse_regression ? new_MP_A : new_MP_B;
      NumericVector oy = reverse_regression ? new_MP_B : new_MP_A;
      CharacterVector fill_SampleID(m);
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] && oy[j] <= upr[j]){
          inside[j] = 1;
          fill_SampleID[j] = NA_STRING;
        }
      }
      if(SampleID_exists == 1){
        List out = List::create(Named("SampleID") = new_SampleID,
                                Named("MP_B") = ox,
                                Named("MP_A") = oy,
                                Named("prediction") = ny,
                                Named("lwr") = lwr,
                                Named("upr") = upr,
                                Named("inside") = inside);  
        return out;
      }
      
      List out = List::create(Named("SampleID") = fill_SampleID,
                              Named("MP_B") = ox,
                              Named("MP_A") = oy,
                              Named("prediction") = ny,
                              Named("lwr") = lwr,
                              Named("upr") = upr,
                              Named("inside") = inside);
      
      return out;
    }
    
    else if(MP_A_exists == 0){
      CharacterVector fill_SampleID(m);
      NumericVector fill_MP_A(m);
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        fill_MP_A[j] = NA_REAL;
        fill_SampleID[j] = NA_STRING;
      }
      if(SampleID_exists == 1){
        List out = List::create(Named("SampleID") = new_SampleID,
                                Named("MP_B") = new_MP_B,
                                Named("MP_A") = fill_MP_A,
                                Named("prediction") = ny,
                                Named("lwr") = lwr,
                                Named("upr") = upr,
                                Named("inside") = inside);
        
        return out;
      }
      
      List out = List::create(Named("SampleID") = fill_SampleID,
                              Named("MP_B") = new_MP_B,
                              Named("MP_A") = fill_MP_A,
                              Named("prediction") = ny,
                              Named("lwr") = lwr,
                              Named("upr") = upr,
                              Named("inside") = inside);
      
      return out;
    }
  }
  List out = List::create(Named("output") = NA_STRING);
  return out;
}
  
  

