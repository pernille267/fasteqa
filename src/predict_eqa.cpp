#include <Rcpp.h>
using namespace Rcpp;

//' Prediction Interval Estimation for EQA Data via Deming or OLS Regression
//'
//' @title Prediction Interval Estimation for EQA Data via Deming or OLS Regression
//' @name predict_eqa
//'
//' @description This function estimates prediction intervals for External Quality Assessment (EQA) data by applying either Deming (\code{method = 'fg'} or \code{method = 'clsi'}) or Ordinary Least Squares (\code{method = 'ols'}) regression methodologies. It is specifically designed for a single In Vitro Diagnostic Medical Device (IVD-MD) comparison. For processing multiple IVD-MD comparisons simultaneously, consider using the \code{estimate_prediction_data()} function from the commutability package, which executes this function on a per IVD-MD comparison basis.   
//'
//' @param data A \code{list} or \code{data.table} representing mean-of-replicates clinical sample data with the \code{list} elements or \code{data.table} columns: \code{SampleID}, \code{MP_A}, and \code{MP_B}.
//' @param new_data A \code{list} or \code{data.table} representing mean-of-replicates EQA material / reference material data. This should at least include \code{MP_B} but may also contain \code{SampleID} and \code{MP_A} depending on the desired output structure:
//' \itemize{
//'   \item{\code{Minimum requirement}: }{new_data must include \code{MP_B}}
//'   \item{\code{method = 'ols'}: }{new_data must include both \code{MP_B} and \code{MP_A}}
//'   \item{\code{Inside checks performed}: }{new_data must include both \code{MP_B} and \code{MP_A}}
//'   \item{\code{Sample-wise prediction}: }{new_data must include \code{MP_B} and \code{SampleID} (plus \code{MP_A} if method = 'ols')}
//' }
//' @param imprecision_estimates A \code{list} or \code{data.table} that includes necessary imprecision estimates. If \code{method} is \code{'fg'} or \code{'ols'}, this parameter should contain \code{lambda}. In case of \code{method = 'clsi'}, both \code{lambda} and \code{Var_B} must be included. The function \code{global_precision_estimates()} can be employed to generate suitable input for this parameter directly.
//' @param R An \code{integer} indicating the number of replicates on which new_data is based. The convenience function \code{count_samplewise_replicates()} can be employed to generate suitable input for this parameter directly
//' @param R_ratio An \code{float} indicating the ratio of replicates between that number data and that number new_data are based on. Only relevant if \code{method = 'fg'} and if the number of replicates of data and new_data differs.
//' @param method A \code{string} specifying the method for estimating the prediction intervals. Default is \code{'fg'}. Current possible prediction interval estimation methods are:
//' \itemize{
//'   \item{\code{fg: }}{Implements standard Deming regression with prediction intervals calculated using the Fuller and Gillard methods.}
//'   \item{\code{clsi: }}{Utilizes standard Deming regression, with prediction intervals derived from the CLSI method.}
//'   \item{\code{ols: }}{Implements Ordinary Least Squares regression, where the predictor is selected to minimize the variance of the ignored IVD-MD uncertainty. Specifically, \code{MP_A} is used as the predictor when \code{lambda < 1}, otherwise, \code{MP_B} is utilized.}
//' }
//' @param level A \code{float} representing the confidence level for the prediction intervals. It should be between \code{0} and \code{1}. The default setting is \code{0.99}. Please adjust for simultaneous testing if pointwise prediction intervals are used for classifying more than one EQA material / reference material in the same IVD-MD comparison.
//' @param rounding An \code{integer} specifying the desired decimal places for the predictions and prediction intervals. The default setting is three, offering sufficient precision. The maximum limit is six due to the utilization of floating-point numbers.
//'
//' @details For optimal results, we recommended to include \code{SampleID}, \code{MP_A}, and \code{MP_B} in new_data whenever possible. If only \code{MP_B} is available, the use of \code{method = 'ols'} is not possible. For instance, when constructing prediction bands, only \code{MP_B} might be available. Thus, \code{method = 'ols'} may sadly not be used to estimate pointwise prediction intervals using this function.
//'
//' @return A \code{list} comprising estimated prediction interval data based on user inputs. The output \code{list} may be converted to \code{data.table} by employing the \code{setDT()} method from the \code{data.table} package in R.
//'
//' @examples \dontrun{
//' library(fasteqa)
//' # Simulation parameters for clinical sample data
//' training_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.015, cil = 10, ciu = 70)
//' # Simulation parameters for external quality assessment material data (commutable materials)
//' test_parameters <- list(n = 5, R = 3, cvx = 0.01, cvy = 0.015, cil = 10, ciu = 70)
//' # Simulation of clinical sample data using the simulation parameters 'training_parameters'
//' training_data <- simulate_eqa_data(training_parameters)
//' # Simulation of external quality assessment material data using the simulation parameters 'test_parameters'
//' test_data <- simulate_eqa_data(test_parameters)
//' # Convert ID columns to character type because this is the type accepted
//' training_data$SampleID <- as.character(training_data$SampleID)
//' training_data$ReplicateID <- as.character(training_data$ReplicateID)
//' # Estimate imprecision estimates
//' imprecision <- global_precision_estimates(training_data)
//' # Calculate mean-of-replicates data
//' mean_of_replicates_training_data <- fun_of_replicates(training_data)
//' mean_of_replicates_test_data <- fun_of_replicates(test_data)
//' # Estimate prediction intervals using method = 'fg', R = 3L and R_ratio = 1:
//' prediction_intervals <- predict_eqa(mean_of_replicates_training_data,
//'                                     mean_of_replicates_test_data,
//'                                     imprecision)
//' }

// [[Rcpp::export]]
List predict_eqa(List data, List new_data, List imprecision_estimates, int R = 3, float R_ratio = 1, String method = "fg", float level = 0.99, int rounding = 3) {
  
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
    float t_quantile = R::qt((1 - level) / 2, n - 2, 0, 0);
    float mu_hat = 0;
    for(int i = 0; i < n; ++i){
      float latent_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      mu_hat += latent_i / n_float;
      latent[i] = latent_i;
    }
    float hvar = (msyy + (lambda * msxx) - sub_expression_2) / (2 * lambda);
    float vvar = lambda * hvar;
    for(int j = 0; j < m; ++j){
      ny[j] = b0 + b1 * nx[j];
      nl[j] = mu_hat * (1 - (msxy / (b1 * msxx))) + ny[j] * (msxy / (b1 * msxx));
    }
    
    NumericVector var_pred_error(m);
    NumericVector lwr(m);
    NumericVector upr(m);
    IntegerVector inside(m);
    for(int j = 0; j < m; ++j){
      var_pred_error[j] = varb1 * pow((nl[j] - mu_hat), 2) + (varb1 * hvar * R_ratio) + (1 + (1 / n_float)) * (pow(b1, 2) * hvar + vvar) * R_ratio;
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
      mse = mse / (n_float - 2) * R_ratio;
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
      mse = mse / (n_float - 2) * R_ratio;
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
  
  

