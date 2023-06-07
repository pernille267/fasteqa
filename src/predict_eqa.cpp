#include <Rcpp.h>
using namespace Rcpp;

//' Estimate prediction intervals for EQA data with Deming or OLS
//' 
//' @title Estimate prediction intervals for EQA data with Deming or OLS
//' @name predict_eqa
//' 
//' @param data \code{list} or \code{data table} - Mean-of-replicates clinical sample data with elements/columns \code{SampleID}, \code{MP_A} and \code{MP_B}
//' @param new_data \code{list} or \code{data table} - Mean-of-replicates control material data with \code{MP_B}. \code{SampleID} and \code{MP_A} may also be in new_data based on desired output structure:
//' \itemize{
//'   \item{\code{Minimum requirement: }}{Must contain \code{MP_B}}
//'   \item{\code{method = ols: }}{Must at least contain both \code{MP_B} and \code{MP_A}}
//'   \item{\code{Inside checks performed: }}{Must at least contain both \code{MP_B} and \code{MP_A}}
//'   \item{\code{Sample-wise prediction}}{Must at least contain \code{MP_B} and \code{SampleID} (plus \code{MP_A} if method = ols)}
//' }
//' @param imprecision_estimates \code{list} or \code{data table} - Imprecision estimates e.g., that which are outputted by \code{global_precision_estimates()}. Minimum requirement: \code{lambda} and \code{Var_B} must be part of the  
//' @param R integer - Number of replicates, which new_data is based on
//' @param method string - Which method should be used to estimate the prediction intervals. Default is \code{fg}. At the moment, possible values for method is
//' \itemize{
//'   \item{\code{fg: }}{Standard Deming regression, but prediction intervals are calculated using components from J. Gillards work}
//'   \item{\code{clsi: }}{Standard Deming regression, but prediction intervals are calculated based on derivation of Jeff Vaks}
//'   \item{\code{ols: }}{Ordinary least squares regression, but the predictor is chosen so that neglected prediction error variance is reduced}
//' }
//' @param level float - Confidence level of prediction intervals. Should ideally be corrected for simulations testing. Default level is 0.99
//' @param rounding integer - How many decimals should be included in the predictions and prediction intervals. Two or three decimals are often sufficient. Default is 3
//'
//' @description A rich function that calculates prediction intervals for EQA data
//' 
//'
//' @details If possible it is wise to always include \code{SampleID}, \code{MP_A} and \code{MP_B} in new_data. If by some reason only \code{MP_B} is available, one cannot use \code{method = ols}. For example, when construction predicton bands, we will only have \code{MP_B} is available.
//'
//' @return list - prediction interval data based on use-inputs
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   training_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.015, cil = 10, ciu = 70)
//'   test_parameters <- list(n = 5, R = 3, cvx = 0.01, cvy = 0.015, cil = 10, ciu = 70)
//'   training_data <- simulate_eqa_data(training_parameters)
//'   test_data <- simulate_eqa_data(test_parameters)
//'   training_data$$SampleID <- as.character(training_data$SampleID)
//'   training_data$$ReplicateID <- as.character(training_data$ReplicateID)
//'   imprecision <- global_precision_estimates(data)
//'   mean_of_replicates_training_data <- fun_of_replicates(training_data)
//'   mean_of_replicates_test_data <- fun_of_replicates(test_data)
//'   prediction_intervals <- predict_eqa(mean_of_replicates_training data,
//'                                       mean_of_replicates_test_data,
//'                                       imprecision)
//' }

// [[Rcpp::export]]
List predict_eqa(List data, List new_data, List imprecision_estimates, int R = 3, String method = "fg", float level = 0.99, int rounding = 3) {
  
  // Needed global imprecision estimates
  float lambda = imprecision_estimates["lambda"];
  float Var_B = imprecision_estimates["Var_B"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];
  int n = MP_B.length();
  // Converting integer types to float types when needed in division!
  float n_float = static_cast<float>(n);
  float R_float = static_cast<float>(R);
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
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x);
    float msyy = var(y);
    float msxy = 0;
    for(int i = 0; i < n; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy = msxy / (n_float - 1);
    NumericVector latent(n);
    if(MP_B_exists == 0){
      stop("MP_B was not found in new_data. Calculations are terminated");
    }
    NumericVector nx = new_MP_B;
    NumericVector ny(m);
    NumericVector nl(m);
    
    // Fitting the standard Deming model
    float sub_expression_1 = msyy - lambda * msxx;
    float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    float sub_expression_3 = 2 * msxy;
    float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    float b0 = my - b1 * mx;
    
    // Dispersion
    float varb1 = (pow(b1, 2) / (n_float * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
    float hvar = (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda);
    if(hvar < 0.1 * Var_B){
      hvar = 0.1 * Var_B;
    }
    float vvar = lambda * hvar;
    float t_quantile = R::qt((1 - level) / 2, n - 2, 0, 0);
    float mu_hat = 0;
    for(int i = 0; i < n; ++i){
      float latent_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      mu_hat += latent_i / n_float;
      latent[i] = latent_i;
    }
    for(int j = 0; j < m; ++j){
      ny[j] = b0 + b1 * nx[j];
      nl[j] = mu_hat * (1 - (msxy / (b1 * msxx))) + ny[j] * (msxy / (b1 * msxx));
    }
    
    NumericVector var_pred_error(m);
    NumericVector lwr(m);
    NumericVector upr(m);
    IntegerVector inside(m);
    for(int j = 0; j < m; ++j){
      var_pred_error[j] = varb1 * pow((nl[j] - mu_hat), 2) + (varb1 * hvar / R_float) + (1 + (1 / n_float)) * (pow(b1, 2) * hvar + vvar) / R_float;
    }
    for(int j = 0; j < m; ++j){
      lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
      if(lwr[j] < 0){
        lwr[j] = 0;
      }
      upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
    }
    lwr = round(lwr, rounding);
    upr = round(upr, rounding);
    ny = round(ny, rounding);
    if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 1){
      NumericVector oy = new_MP_A;
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
          inside[j] = 1;
        }
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    else if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 0){
      NumericVector oy = new_MP_A;
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
          inside[j] = 1;
        }
        new_SampleID[j] = NA_STRING;
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    
    else if(MP_A_exists == 0 and MP_B_exists == 1 and SampleID_exists == 0){
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        new_SampleID[j] = NA_STRING;
        new_MP_A[j] = NA_REAL;
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    
    else if(MP_A_exists == 0 and MP_B_exists == 1 and SampleID_exists == 1){
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        new_MP_A[j] = NA_REAL;
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    else{
      stop("Unexpected error ecountered when using method = 'fg'!");
    }
  }
  
  // CLSI, EP14 approach
  if(method == "clsi"){
    NumericVector x = MP_B;
    NumericVector y = MP_A;
    float mx = mean(x);
    float my = mean(y);
    float msxx = var(x) * (n_float - 1) / n_float;
    float msyy = var(y) * (n_float - 1) / n_float;
    float msxy = 0;
    for(int i = 0; i < n; ++i){
      msxy = msxy + (x[i] - mx) * (y[i] - my);
    }
    msxy = msxy / n_float;
    if(MP_B_exists == 0){
      stop("MP_B was not found in new_data. Calculations are terminated");
    }
    NumericVector nx = new_MP_B;
    NumericVector ny(m);
    float sub_expression_1 = msyy - lambda * msxx;
    float sub_expression_2 = sqrt(pow(msyy - lambda * msxx, 2) + 4 * lambda * pow(msxy, 2));
    float sub_expression_3 = 2 * msxy;
    float b1 = (sub_expression_1 + sub_expression_2) / sub_expression_3;
    float b0 = my - b1 * mx;
    float varb1 = (pow(b1, 2) / (n_float * pow(msxy, 2))) * ((msxx * msyy) - pow(msxy, 2));
    float hvar = Var_B;
    float vvar = lambda * hvar;
    float t_quantile = R::qt((1 - level) / 2, n * (R - 1), 0, 0);
  
    for(int j = 0; j < m; ++j){
      ny[j] = b0 + b1 * nx[j];
    }
    
    NumericVector var_pred_error(m);
    NumericVector lwr(m);
    NumericVector upr(m);
    IntegerVector inside(m);
    for(int j = 0; j < m; ++j){
      var_pred_error[j] = varb1 * pow((nx[j] - mx), 2) + (1 + (1 / n_float)) * (pow(b1, 2) * hvar + vvar) / R_float;
    }
    for(int j = 0; j < m; ++j){
      lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
      if(lwr[j] < 0){
        lwr[j] = 0;
      }
      upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
    }
    lwr = round(lwr, rounding);
    upr = round(upr, rounding);
    ny = round(ny, rounding);
    if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 1){
      NumericVector oy = new_MP_A;
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
          inside[j] = 1;
        }
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    else if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 0){
      NumericVector oy = new_MP_A;
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
          inside[j] = 1;
        }
        new_SampleID[j] = NA_STRING;
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    
    else if(MP_A_exists == 0 and MP_B_exists == 1 and SampleID_exists == 0){
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        new_SampleID[j] = NA_STRING;
        new_MP_A[j] = NA_REAL;
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    
    else if(MP_A_exists == 0 and MP_B_exists == 1 and SampleID_exists == 1){
      // Checks whether observed value is inside prediction interval
      for(int j = 0; j < m; ++j){
        inside[j] = NA_INTEGER;
        new_MP_A[j] = NA_REAL;
      }
      List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
      return out;
    }
    else{
      stop("Unexpected error ecountered when using method = 'clsi'!");
    }
  }
  if(method == "ols"){
    
    // if lambda > 1, set x = MP_B and y = MP_A
    if(lambda > 1){
      NumericVector x = MP_B;
      NumericVector y = MP_A;
      float sxx = var(x) * (n_float - 1);
      float mx = mean(x);
      float my = mean(y);
      float sxy = 0;
      float mse = 0;
      for(int i = 0; i < n; ++i){
        sxy = sxy + (x[i] - mx) * (y[i] - my);
      }
      float b1 = sxy / sxx;
      float b0 = my - b1 * mx;
      
      for(int i = 0; i < n; ++i){
        mse = mse + pow(y[i] - b0 - b1 * x[i], 2);
      }
      mse = mse / (n_float - 2) / R_float;
      if(MP_B_exists == 0){
        stop("MP_B was not found in new_data. Calculations are terminated");
      }
      NumericVector nx = new_MP_B;
      NumericVector ny(m);
      for(int j = 0; j < m; ++j){
        ny[j] = b0 + b1 * nx[j];
      }
      NumericVector var_pred_error(m);
      NumericVector lwr(m);
      NumericVector upr(m);
      IntegerVector inside(m);
      for(int j = 0; j < m; ++j){
        var_pred_error[j] = mse * (1 + (1 / n_float) + pow(nx[j] - mx, 2) / sxx);
      }
      float t_quantile = R::qt((1 - level) / 2, n - 2, 0, 0);
      for(int j = 0; j < m; ++j){
        lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
        if(lwr[j] < 0){
          lwr[j] = 0;
        }
        upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
      }
      lwr = round(lwr, rounding);
      upr = round(upr, rounding);
      ny = round(ny, rounding);
      
      if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 1){
        NumericVector oy = new_MP_A;
        // Checks whether observed value is inside prediction interval
        for(int j = 0; j < m; ++j){
          if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
            inside[j] = 1;
          }
        }
        List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
        return out;
      }
      else if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 0){
        NumericVector oy = new_MP_A;
        // Checks whether observed value is inside prediction interval
        for(int j = 0; j < m; ++j){
          if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
            inside[j] = 1;
          }
          new_SampleID[j] = NA_STRING;
        }
        List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
        return out;
      }
      else if(MP_A_exists == 0 or MP_B_exists == 0){
        stop("new_data must contain both MP_A and MP_B measurements when method = 'ols'. Calculations are terminated");
      }
      else{
        stop("Unexpected error ecountered when using method = 'ols'!");
      }
    }
    // if lambda < 1, set x = MP_A and y = MP_B
    else{
      NumericVector x = MP_A;
      NumericVector y = MP_B;
      float sxx = var(x) * (n_float - 1);
      float mx = mean(x);
      float my = mean(y);
      float sxy = 0;
      float mse = 0;
      for(int i = 0; i < n; ++i){
        sxy = sxy + (x[i] - mx) * (y[i] - my);
      }
      float b1 = sxy / sxx;
      float b0 = my - b1 * mx;
      
      for(int i = 0; i < n; ++i){
        mse = mse + pow(y[i] - b0 - b1 * x[i], 2);
      }
      mse = mse / (n_float - 2) / R_float;
      if(MP_A_exists == 0){
        stop("MP_A was not found in new_data. Calculations are terminated");
      }
      NumericVector nx = new_MP_A;
      NumericVector ny(m);
      for(int j = 0; j < m; ++j){
        ny[j] = b0 + b1 * nx[j];
      }
      NumericVector var_pred_error(m);
      NumericVector lwr(m);
      NumericVector upr(m);
      IntegerVector inside(m);
      for(int j = 0; j < m; ++j){
        var_pred_error[j] = mse * (1 + (1 / n_float) + pow(nx[j] - mx, 2) / sxx);
      }
      float t_quantile = R::qt((1 - level) / 2, n - 2, 0, 0);
      for(int j = 0; j < m; ++j){
        lwr[j] = ny[j] - t_quantile * sqrt(var_pred_error[j]);
        if(lwr[j] < 0){
          lwr[j] = 0;
        }
        upr[j] = ny[j] + t_quantile * sqrt(var_pred_error[j]);
      }
      lwr = round(lwr, rounding);
      upr = round(upr, rounding);
      ny = round(ny, rounding);
      if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 1){
        NumericVector oy = new_MP_B;
        // Checks whether observed value is inside prediction interval
        for(int j = 0; j < m; ++j){
          if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
            inside[j] = 1;
          }
        }
        List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
        return out;
      }
      else if(MP_A_exists == 1 and MP_B_exists == 1 and SampleID_exists == 0){
        NumericVector oy = new_MP_B;
        // Checks whether observed value is inside prediction interval
        for(int j = 0; j < m; ++j){
          if(oy[j] >= lwr[j] and oy[j] <= upr[j]){
            inside[j] = 1;
          }
          new_SampleID[j] = NA_STRING;
        }
        List out = List::create(Named("SampleID") = new_SampleID, Named("MP_B") = new_MP_B, Named("MP_A") = new_MP_A, Named("prediction") = ny, Named("lwr") = lwr, Named("upr") = upr,  Named("inside") = inside);
        return out;
      }
      else if(MP_A_exists == 0 or MP_B_exists == 0){
        stop("new_data must contain both MP_A and MP_B measurements when method = 'ols'");
      }
      else{
        stop("Unexpected error ecountered when using method = 'ols'!");
      }
    }
  }
  List out = List::create(Named("output") = NA_STRING);
  return out;
}
  
  

