#include <Rcpp.h>

using namespace Rcpp;

//' @title Prediction Interval Estimation for EQA Data via Deming or OLS Regression
//' @name predict_eqa
//'
//' @description
//' This function estimates pointwise prediction intervals for
//' External Quality Assessment (EQA) data by applying either
//' Deming (\code{method = 'fg'} or \code{method = 'clsi'}) or
//' Ordinary Least Squares (\code{method = 'ols'}) regression methodologies.
//' It is specifically designed for a single In Vitro Diagnostic Medical Device (IVD-MD) comparison.
//' For processing multiple IVD-MD comparisons simultaneously, consider using the
//' \code{estimate_prediction_data()} function from the \code{commutability} package,
//' which executes this function on a per IVD-MD comparison basis.   
//'
//' @param data A \code{list} or \code{data.table} representing mean-of-replicates clinical sample data
//'        with the \code{list} elements or \code{data.table} columns: \code{SampleID}, \code{MP_A}, and \code{MP_B}.
//' @param new_data A \code{list} or \code{data.table} representing mean-of-replicates EQA material / reference material data.
//'        This should at least include \code{MP_B} but may also contain \code{SampleID} and \code{MP_A} depending on the desired output structure:
//' \itemize{
//'   \item \code{Minimum requirement}: \code{new_data} must include \code{MP_B}.
//'   \item \code{Inside checks performed}: \code{new_data} must include both \code{MP_B} and \code{MP_A}.
//'   \item \code{Sample-wise prediction}: \code{new_data} must include \code{MP_B} and \code{SampleID}.
//' }
//' @param imprecision_estimates A \code{list} or \code{data.table}.
//'        Must include \code{Var_B} and \code{lambda}. See details.
//' @param R An \code{integer} indicating the average number of replicates on which new_data is based.
//'        The convenience function \code{count_samplewise_replicates()} can be employed to generate
//'        suitable input for this parameter directly.
//' @param R_ratio An \code{double}. The ratio of the number of replicates used in \code{new_data} and \code{data}.
//'        Only relevant if \code{method = 'fg'} or \code{method = 'ols'}.
//'        If the same number of replicates is used for measurements in \code{data} and \code{new_data},
//'        use the default value of \code{R_ratio = 1}.
//' @param method A \code{character} string. The method used to estimate the prediction intervals.
//'        The default is \code{'fg'}. Current possible prediction interval estimation approaches are:
//' \itemize{
//'   \item \code{fg: } Implements standard Deming regression with prediction intervals calculated using the 'Fuller & Gillard' approach. See details.
//'   \item \code{clsi: } Utilizes standard Deming regression, with prediction intervals derived from the CLSI method. See details.
//'   \item \code{ols: } Implements Ordinary Least Squares regression.
//' }
//' @param level A \code{double}. The nominal confidence level for the prediction intervals.
//'        Must be a value between \code{0} and \code{1}. The default setting is \code{0.99}.
//'        Please adjust for simultaneous testing if pointwise prediction intervals are used for
//'        classifying more than one EQA material / reference material in the same IVD-MD comparison.
//' @param allow_reverse_regression A \code{logical} value. If set to \code{TRUE}, and \code{method = 'ols'},
//'        the response and predictor change roles if \code{lambda < 1}.
//' @param rounding An \code{integer} specifying the desired number of decimal places for the predictions and prediction intervals.
//'        The default setting is \code{3L}, offering sufficient precision for most applications.
//'        The maximum limit is \code{12L}.
//'
//' @details
//' For commutability evaluation, inclusion of all \code{SampleID}, \code{MP_A}, and \code{MP_B}
//' in \code{new_data} is required. If only prediction bands are required, one only need to include \code{MP_B} part of \code{new_data}.
//' 
//' Imprecision Estimates
//' 
//' Estimating prediction intervals may require certain imprecision estimates, which can be estimated from the raw data.
//' Which imprecision estimates that are required depends on which method that is used to construct the prediction intervals.
//' For \code{method = 'fg'} or \code{method = 'ols'}, \code{lambda} is required. This estimate is calculated by
//' 
//' \eqn{\lambda = \frac{\hat{\sigma}_v^2}{\hat{\sigma}_h^2}}
//' 
//' where \eqn{\hat{\sigma}_v^2} and \eqn{\hat{\sigma}_h^2} are pooled variances.
//' If differences in non-selectivity (DINS) or equation error cannot be ruled out,
//' \eqn{\lambda} will generally underestimate the true \eqn{\Lambda}. Do not forget this.
//' 
//' For \code{method = 'clsi'}, both \code{Var_B} and \code{lambda} are required. \code{lambda} is calculated as before,
//' and \code{Var_B} is just \eqn{\hat{\sigma}_h^2}.
//' 
//' Note that both \code{lambda} and \code{Var_B} always must be given. However, if one of them
//' is not necessary to estimate prediction intervals, it may take an arbitrary value.
//' It is nevertheless recommended to always include actual values \code{lambda} and \code{Var_B},
//' independently of whether they are actual used or not.
//' 
//' Note also that \code{global_precision_estimates()} can be used to calculate both
//' \code{lambda} and \code{Var_B}.
//' 
//'
//' @return 
//' A \code{list} comprising estimated prediction interval data based on user inputs.
//' This list will contain the following elements:
//' \itemize{
//'     \item \code{SampleID: } The ID(s) of the evaluated material(s).
//'     \item \code{MP_B: } The measurement result(s) of the second IVD-MD in the \code{MP_A} - \code{MP_B} comparison (predictor).
//'     \item \code{MP_A: } The measurement result(s) of the first IVD-MD in the \code{MP_A} - \code{MP_B} comparison (response).
//'     \item \code{prediction: } The predicted value(s) of \code{MP_A} given the value(s) of \code{MP_B}.
//'     \item \code{lwr: } The lower limit(s) of the estimated prediction interval(s) of \code{MP_A} given the value(s) of \code{MP_B}.
//'     \item \code{upr: } The upper limit(s) of the estimated prediction interval(s) of \code{MP_A} given the value(s) of \code{MP_B}.
//'     \item \code{inside: } Inside checks. If \code{1}, \code{MP_A} is inside the estimated prediction interval for \code{MP_A}.
//'     On the other hand, if \code{0}, \code{MP_A} is outside the estimated prediction interval for \code{MP_A}.
//' }
//' The output \code{list} may be converted to \code{data.table} by employing the \code{setDT()} method from the \code{data.table} package in R.
//'
//' @examples \dontrun{
//' library(fasteqa)
//' # Simulation parameters for clinical sample data
//' training_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.015, cil = 10, ciu = 70)
//' # Simulation parameters for external quality assessment material data (commutable materials)
//' test_parameters <- list(n = 5, R = 3, cvx = 0.01, cvy = 0.015, cil = 10, ciu = 70)
//' # Simulation of clinical sample data using the simulation parameters 'training_parameters'
//' training_data <- simulate_eqa_data(training_parameters)
//' # Simulation of external quality assessment material data based on 'test_parameters'
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
List predict_eqa(List data, List new_data, List imprecision_estimates,
                 int R = 3, double R_ratio = 1, String method = "fg",
                 double level = 0.99, bool allow_reverse_regression = true, int rounding = 3) {
  
  // Extract required information from imprecision_estimates
  double lambda = imprecision_estimates["lambda"];
  double Var_B = imprecision_estimates["Var_B"];
  
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
  
  

