#include <Rcpp.h>
using namespace Rcpp;

//' Merge all commutability evaluation computations efficiently into one object
//' 
//' @title Merge all computations efficiently into one object
//' @name merge_results
//' @param pb_data \code{list} or \code{data table} - Prediction band data, which must at least contain \code{comparison}, \code{predictor}, \code{prediction}, \code{lwr} and \code{upr}
//' @param ce_data \code{list} or \code{data table} - Commutability evaluation data for evaluated material, which must at least contain \code{comparison}, \code{SampleID}, \code{MP_B}, \code{MP_A}, \code{prediction}, \code{lwr}, \code{upr} and \code{inside}
//' @param zeta_data \code{list} or \code{data table} - Zeta estimates, which must at least contain \code{comparison}, \code{zeta}, \code{lwr}, \code{upr}, \code{zeta_critical}, \code{zeta_conclusion}
//' @param imprecision_data \code{list} or \code{data table} - global imprecision estimates, which must at least contain \code{comparison}, \code{CV_A}, \code{CV_A_lwr}, \code{CV_A_upr}, \code{CV_B}, \code{CV_B_lwr}, \code{CV_B_upr}, \code{lambda}, \code{lambda_lwr} and \code{lambda_upr}. This one may be omitted, if one is uninterested in viewing the imprecision estimates
//' @param rounding \code{integer} - How many decimals should be included in the float type variables? The implicit maximum is 6 decimals, meaning that any integer larger than 6 will produce a warning
//' @param include_imprecision_estimates \code{boolean} - Should imprecision estimates be part of the merged results. Default is \code{FALSE}.
//' @param silence \code{integer} - Should progress reports be printed? For yes, pass \code{0}. For no, pass \code{1}. \code{1} is default.
//' 
//' @description Merge all the essential components (e.g., pb_data, ce_data, zeta_data) of the commutability evaluation analysis into one object, so it will be easier to plot and present commutability evaluation results
//' 
//' @details The merging is done with respect to the comparison column that needs to be common in every data input. The effects of missing values in the comparison variable may yield unexpected results or errors, so ensure that at least this column, is NA-free. 
//' 
//' @return A \code{list} with two main components, merged data results for commutability evaluation of control materials and merged results for prediction bands. The latter may be used to plot commutability evaluation plots
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   data <- simulate_data_eqa(list(n = 25, R = 3, cvx = 0.06, cvy = 0.04))
//'   estimate_zeta(data)
//' }
//'

// [[Rcpp::export]]
List merge_results(List pb_data, List ce_data, List zeta_data, List imprecision_data, int rounding = 3, bool include_imprecision_estimates = true, int silence = 1){
  
  // pb_data as reference
  CharacterVector comparison = pb_data["comparison"];
  CharacterVector copy_comparison = clone(comparison);
  CharacterVector comparison_sorted = copy_comparison.sort();
  CharacterVector unique_comparisons = sort_unique(comparison);
  int n = unique_comparisons.size();
  int N = comparison.size();
  if(silence == 0){
    Rcout << "Progress report 1..." << "\n";
    Rcout << "Number of unique comparisons:  " << n << "\n"; 
    Rcout << "Number of rows:  " << N << "\n";
    Rcout << "__________________________________________________________________" << "\n";
  }
  IntegerVector n_comparison(n);
  IntegerVector pb_new_indices(N);
  // Count absolute frequency of each unique comparison (Part 1 / *)
  int c = 1; // within comparison counter
  int j = 0; // unique comparison counter
  for(int i = 0; i < N; ++i){
    if(i == (N - 1)){
      n_comparison[j] = c;
      break;
    }
    
    if(comparison_sorted[i] == comparison_sorted[i + 1]){
      ++c;
    }
    
    else if(comparison_sorted[i] != comparison_sorted[i + 1]){
      n_comparison[j] = c;
      c = 1;
      ++j;
    }
  }
  
  if(silence == 0){
    Rcout << "Progress report 2..." << "\n";
    Rcout << "There are " << n_comparison << " unique comparisons in input data" << "\n";  
    Rcout << "__________________________________________________________________" << "\n";
  }
  
  // Reset counter
  c = 0;
  
  // Correct placement of numbers after sorting comparison
  for(int i = 0; i < n; ++i){
    CharacterVector to_match(1);
    to_match[0] = unique_comparisons[i];
    LogicalVector matches = in(comparison, to_match);
    int n_matches = sum(matches);
    IntegerVector sub_indices(n_matches);
    int d = 0; // match counter
    for(int j = 0; j < N; ++j){
      if(matches[j] == 1){
        sub_indices[d] = j;
        ++d;  
      }
    }
    
    if(n_matches != n_comparison[i]){
      stop("n_matches is different from n_comparison. Something is wrong...");
    }
    
    for(int j = 0; j < n_comparison[i]; ++j){
      pb_new_indices[c] = sub_indices[j];
      ++c;
    }
    
    if(silence == 0){
      Rcout << "Progress report 3..." << "\n";
      Rcout << "We are matching: " << to_match << "\n";
      Rcout << "Matches for " << to_match << " in original comparison: " << matches << "\n";
      Rcout << "__________________________________________________________________" << "\n";
    }
    
  }
  if(silence == 0){
    Rcout << "Progress report 4..." << "\n";
    Rcout << "pb_new_indices are " << ": "<< pb_new_indices << "\n";
    Rcout << "__________________________________________________________________" << "\n";
  }
  
  NumericVector predictor = pb_data["predictor"];
  NumericVector prediction = pb_data["prediction"];
  NumericVector lwr = pb_data["lwr"];
  NumericVector upr = pb_data["upr"];
  
  NumericVector new_predictor(N);
  NumericVector new_prediction(N);
  NumericVector new_lwr(N);
  NumericVector new_upr(N);
  
  // Update whole data based on sorting of comparison
  for(int i = 0; i < N; ++i){
    new_predictor[i] = predictor[pb_new_indices[i]];
    new_prediction[i] = prediction[pb_new_indices[i]];
    new_lwr[i] = lwr[pb_new_indices[i]];
    new_upr[i] = upr[pb_new_indices[i]];
  }
  
  // ce_data as reference
  CharacterVector comparison_ce = ce_data["comparison"];
  CharacterVector copy_comparison_ce = clone(comparison_ce);
  CharacterVector comparison_ce_sorted = copy_comparison_ce.sort();
  CharacterVector unique_comparisons_ce = sort_unique(comparison_ce);
  int m = unique_comparisons_ce.size();
  int M = comparison_ce.size();
  IntegerVector ce_new_indices(M);
  IntegerVector m_comparison(m);
  
  if(m != n){
    stop("The number of comparisons in ce_data does not match the number of comparisons in pb_data!");
  }
  
  // Reset the two counters
  c = 1; // within counter for comparison
  j = 0; // number of comparisons counter
  
  for(int i = 0; i < M; ++i){
    if(i == (M - 1)){
      m_comparison[j] = c;
      break;
    }
    if(comparison_ce_sorted[i] == comparison_ce_sorted[i + 1]){
      ++c;
    }
    else if(comparison_ce_sorted[i] != comparison_ce_sorted[i + 1]){
      m_comparison[j] = c;
      c = 1;
      ++j;
    }
  }
  
  // reset counter
  c = 0;
  if(silence == 0){
    Rcout << "Progress report 5..." << "\n";
    Rcout << "m_comparison is : " << m_comparison << "\n";
    Rcout << "__________________________________________________________________" << "\n";
  }
  
  
  for(int i = 0; i < m; ++i){
    CharacterVector to_match(1);
    to_match[0] = unique_comparisons_ce[i];
    LogicalVector matches = in(comparison_ce, to_match);
    int m_matches = sum(matches);
    IntegerVector sub_indices(m_matches);
    int d = 0;
    for(int j = 0; j < M; ++j){
      if(matches[j] == 1){
        sub_indices[d] = j;
        ++d;  
      }
    }
    
    if(m_matches != m_comparison[i]){
      stop("m_matches is different from m_comparison. Something is wrong...");
    }
    
    for(int j = 0; j < m_comparison[i]; ++j){
      ce_new_indices[c] = sub_indices[j];
      ++c;
    }
    
    if(silence == 0){
      Rcout << "Progress report 6..." << "\n";
      Rcout << "We are matching: " << to_match << "\n";
      Rcout << "Matches for " << to_match << " in original comparison: " << matches << "\n";
      Rcout << "__________________________________________________________________" << "\n";
    }
  }
  if(silence == 0){
    Rcout << "Progress report 7..." << "\n";
    Rcout << "ce_new_indices are " << ": "<< ce_new_indices << "\n";
    Rcout << "__________________________________________________________________" << "\n";
  }
  
  CharacterVector SampleID = ce_data["SampleID"];
  NumericVector MS_A = ce_data["MP_A"];
  NumericVector MS_B = ce_data["MP_B"];
  NumericVector ce_prediction = ce_data["prediction"];
  NumericVector ce_lwr = ce_data["lwr"];
  NumericVector ce_upr = ce_data["upr"];
  IntegerVector ce_inside = ce_data["inside"];
  
  if(silence == 0){
    Rcout << "Progress report 8..." << "\n";
    Rcout << "SampleIDs of ce_data: " << SampleID << "\n";
    Rcout << "__________________________________________________________________" << "\n";
  }
  
  // New ordering after ordering of 'comparison'
  CharacterVector new_SampleID(M);
  NumericVector new_MS_A(M);
  NumericVector new_MS_B(M);
  NumericVector new_ce_prediction(M);
  NumericVector new_ce_lwr(M);
  NumericVector new_ce_upr(M);
  IntegerVector new_ce_inside(M);
  
  // Update whole data based on sorting of comparison
  for(int i = 0; i < M; ++i){
    new_SampleID[i] = SampleID[ce_new_indices[i]];
    new_MS_A[i] = MS_A[ce_new_indices[i]];
    new_MS_B[i] = MS_B[ce_new_indices[i]];
    new_ce_prediction[i] = ce_prediction[ce_new_indices[i]];
    new_ce_lwr[i] = ce_lwr[ce_new_indices[i]];
    new_ce_upr[i] = ce_upr[ce_new_indices[i]];
    new_ce_inside[i] = ce_inside[ce_new_indices[i]];
  }
  
  //////////////////////////////////////////////
  // Expand zeta_data with respect to pb_data //
  //////////////////////////////////////////////
  
  CharacterVector comparison_zeta = zeta_data["comparison"];
  CharacterVector copy_comparison_zeta = clone(comparison_zeta);
  CharacterVector comparison_zeta_sorted = copy_comparison_zeta.sort();
  IntegerVector pb_new_indices_zeta(n);
  if(n != comparison_zeta.size()){
    stop("Each comparison must correspond to exactly one zeta value. This is not the case, so calculations are terminated");
  }
  
  // Correct placement of numbers after sorting comparison
  for(int i = 0; i < n; ++i){
    CharacterVector to_match(1);
    to_match[0] = comparison_zeta_sorted[i];
    LogicalVector matches = in(comparison_zeta, to_match);
    for(j = 0; j < n; ++j){
       if(matches[j] == 1){
         pb_new_indices_zeta[i] = j;
       } 
    }
  }
  
  // Get data from input list
  NumericVector zeta_pe = zeta_data["zeta"];
  NumericVector zeta_lwr = zeta_data["lwr"];
  NumericVector zeta_upr = zeta_data["upr"];
  NumericVector zeta_crit = zeta_data["zeta_critical"];
  IntegerVector zeta_conc = zeta_data["zeta_conclusion"];
  
  // New order based on comparison sorted
  NumericVector new_order_zeta_pe(n);
  NumericVector new_order_zeta_lwr(n);
  NumericVector new_order_zeta_upr(n);
  NumericVector new_order_zeta_crit(n);
  IntegerVector new_order_zeta_conc(n);
  
  // Update whole data based on sorting of comparison
  for(int i = 0; i < n; ++i){
    new_order_zeta_pe[i] = zeta_pe[pb_new_indices_zeta[i]];
    new_order_zeta_lwr[i] = zeta_lwr[pb_new_indices_zeta[i]];
    new_order_zeta_upr[i] = zeta_upr[pb_new_indices_zeta[i]];
    new_order_zeta_crit[i] = zeta_crit[pb_new_indices_zeta[i]];
    new_order_zeta_conc[i] = zeta_conc[pb_new_indices_zeta[i]];
  }
  
  CharacterVector new_comparison_zeta(N);
  NumericVector new_zeta_pe(N);
  NumericVector new_zeta_lwr(N);
  NumericVector new_zeta_upr(N);
  NumericVector new_zeta_crit(N);
  IntegerVector new_zeta_conc(N);
  
  // reset counter
  c = 0;
  
  // Updated order for imprecision estimates
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n_comparison[i]; ++j){
      new_comparison_zeta[c] = comparison_zeta_sorted[i];
      new_zeta_pe[c] = new_order_zeta_pe[i];
      new_zeta_lwr[c] = new_order_zeta_lwr[i];
      new_zeta_upr[c] = new_order_zeta_upr[i];
      new_zeta_crit[c] = new_order_zeta_crit[i];
      new_zeta_conc[c] = new_order_zeta_conc[i];
      ++c;
    }  
  }
  
  //////////////////////////////////
  // expand zeta_data for ce_data //
  //////////////////////////////////
  
  // Expand zeta_data with respect to ce_data
  IntegerVector ce_new_indices_zeta(m);
  if(m != comparison_zeta.size()){
    stop("Each comparison in zeta_data must correspond to exactly one zeta value. This is not the case in zeta_data, so calculations are terminated");
  }
  
  // Correct placement of numbers after sorting comparison
  for(int i = 0; i < m; ++i){
    ce_new_indices_zeta[i] = pb_new_indices_zeta[i];
  }
  
  // New order based on comparison sorted
  NumericVector ce_new_order_zeta_pe(m);
  NumericVector ce_new_order_zeta_lwr(m);
  NumericVector ce_new_order_zeta_upr(m);
  NumericVector ce_new_order_zeta_crit(m);
  IntegerVector ce_new_order_zeta_conc(m);
  
  // Update whole data based on sorting of comparison
  for(int i = 0; i < m; ++i){
    ce_new_order_zeta_pe[i] = zeta_pe[ce_new_indices_zeta[i]];
    ce_new_order_zeta_lwr[i] = zeta_lwr[ce_new_indices_zeta[i]];
    ce_new_order_zeta_upr[i] = zeta_upr[ce_new_indices_zeta[i]];
    ce_new_order_zeta_crit[i] = zeta_crit[ce_new_indices_zeta[i]];
    ce_new_order_zeta_conc[i] = zeta_conc[ce_new_indices_zeta[i]];
  }
  
  CharacterVector ce_new_comparison_zeta(M);
  NumericVector ce_new_zeta_pe(M);
  NumericVector ce_new_zeta_lwr(M);
  NumericVector ce_new_zeta_upr(M);
  NumericVector ce_new_zeta_crit(M);
  IntegerVector ce_new_zeta_conc(M);
  
  // reset counter
  c = 0;
  
  // Updated order for imprecision estimates
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m_comparison[i]; ++j){
      ce_new_comparison_zeta[c] = comparison_zeta_sorted[i];
      ce_new_zeta_pe[c] = ce_new_order_zeta_pe[i];
      ce_new_zeta_lwr[c] = ce_new_order_zeta_lwr[i];
      ce_new_zeta_upr[c] = ce_new_order_zeta_upr[i];
      ce_new_zeta_crit[c] = ce_new_order_zeta_crit[i];
      ce_new_zeta_conc[c] = ce_new_order_zeta_conc[i];
      ++c;
    }  
  }
  
  
  if(include_imprecision_estimates == true){
    
    // Expand imprecision_data for PB data
    
    CharacterVector comparison_imprecision = imprecision_data["comparison"];
    CharacterVector copy_comparison_imprecision = clone(comparison_imprecision);
    CharacterVector comparison_imprecision_sorted = copy_comparison_imprecision.sort();
    IntegerVector pb_new_indices_imprecision(n);
    if(n != comparison_imprecision.size()){
      stop("Each comparison must correspond to exactly one row of imprecision estimates. This is not the case, so calculations are terminated");
    }
    
    
    // reordering of other variables based on order of comparison
    for(int i = 0; i < n; ++i){
      CharacterVector to_match(1);
      to_match[0] = comparison_imprecision_sorted[i];
      LogicalVector matches = in(comparison_imprecision, to_match);
      for(j = 0; j < n; ++j){
        if(matches[j] == 1){
          pb_new_indices_imprecision[i] = j;
        } 
      }
    }
    
    if(silence == 0){
      Rcout << "Progress report 9..." << "\n";
      Rcout << "New indices based on sorted comparisons are : " << pb_new_indices_imprecision << "\n";
      Rcout << "Original order of comparisons are : " <<  "\n";
      Rcout << "---> " << comparison_imprecision << "\n";
      Rcout << "Sorted order of comparisons are :" << "\n";
      Rcout << "---> "<< comparison_imprecision_sorted << "\n";
      Rcout << "__________________________________________________________________" << "\n";  
    }
    
    // Get data from input list
    NumericVector CV_A = imprecision_data["CV_A"];
    NumericVector CV_A_lwr = imprecision_data["CV_A_lwr"];
    NumericVector CV_A_upr = imprecision_data["CV_A_upr"];
    NumericVector CV_B = imprecision_data["CV_B"];
    NumericVector CV_B_lwr = imprecision_data["CV_B_lwr"];
    NumericVector CV_B_upr = imprecision_data["CV_B_upr"];;
    NumericVector lambda = imprecision_data["lambda"];
    NumericVector lambda_lwr = imprecision_data["lambda_lwr"];
    NumericVector lambda_upr = imprecision_data["lambda_upr"];;
    
    // New order based on comparison sorted
    NumericVector new_order_CV_A(n);
    NumericVector new_order_CV_A_lwr(n);
    NumericVector new_order_CV_A_upr(n);
    NumericVector new_order_CV_B(n);
    NumericVector new_order_CV_B_lwr(n);
    NumericVector new_order_CV_B_upr(n);
    NumericVector new_order_lambda(n);
    NumericVector new_order_lambda_lwr(n);
    NumericVector new_order_lambda_upr(n);
    
    // Updated order for imprecision estimates
    for(int i = 0; i < n; ++i){
      new_order_CV_A[i] = CV_A[pb_new_indices_imprecision[i]];
      new_order_CV_A_lwr[i] = CV_A_lwr[pb_new_indices_imprecision[i]];
      new_order_CV_A_upr[i] = CV_A_upr[pb_new_indices_imprecision[i]];
      new_order_CV_B[i] = CV_B[pb_new_indices_imprecision[i]];
      new_order_CV_B_lwr[i] = CV_B_lwr[pb_new_indices_imprecision[i]];
      new_order_CV_B_upr[i] = CV_B_upr[pb_new_indices_imprecision[i]];
      new_order_lambda[i] = lambda[pb_new_indices_imprecision[i]];
      new_order_lambda_lwr[i] = lambda_lwr[pb_new_indices_imprecision[i]];
      new_order_lambda_upr[i] = lambda_upr[pb_new_indices_imprecision[i]];  
    }
    
    if(silence == 0){
      Rcout << "Progress report 10..." << "\n";
      Rcout << "CV_A unordered: " << CV_A << "\n";
      Rcout << "CV_A ordered: " << new_order_CV_A << "\n";
      Rcout << "CV_B unordered: " << CV_B << "\n";
      Rcout << "CV_B ordered: " << new_order_CV_B << "\n";
      Rcout << "lambda unordered: " << lambda << "\n";
      Rcout << "lambda ordered: " << new_order_lambda << "\n";
      Rcout << "__________________________________________________________________" << "\n";  
    }
    
    ////////////////////////////////////////////////////
    // expansion of imprecision_data based on pb_data //
    ////////////////////////////////////////////////////
    
    CharacterVector new_comparison_imprecision(N);
    NumericVector new_CV_A(N);
    NumericVector new_CV_A_lwr(N);
    NumericVector new_CV_A_upr(N);
    NumericVector new_CV_B(N);
    NumericVector new_CV_B_lwr(N);
    NumericVector new_CV_B_upr(N);
    NumericVector new_lambda(N);
    NumericVector new_lambda_lwr(N);
    NumericVector new_lambda_upr(N);
    
    // reset counter
    c = 0;
    
    // Updated imprecision estimates based on the order found in comparison
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < n_comparison[i]; ++j){
        new_comparison_imprecision[c] = comparison_imprecision_sorted[i];
        new_CV_A[c] = new_order_CV_A[i];
        new_CV_A_lwr[c] = new_order_CV_A_lwr[i];
        new_CV_A_upr[c] = new_order_CV_A_upr[i];
        new_CV_B[c] = new_order_CV_B[i];
        new_CV_B_lwr[c] = new_order_CV_B_lwr[i];
        new_CV_B_upr[c] = new_order_CV_B_upr[i];
        new_lambda[c] = new_order_lambda[i];
        new_lambda_lwr[c] = new_order_lambda_lwr[i];
        new_lambda_upr[c] = new_order_lambda_upr[i];
        ++c;
      }
    }
    
    if(silence == 0){
      Rcout << "Progress report 10..." << "\n";
      Rcout << "lambda values with old order : " << lambda <<"\n";
      Rcout << "lambda values with new order : " << new_lambda << "\n";   
      Rcout << "__________________________________________________________________" << "\n";
    }
    
    /////////////////////////////////////////
    // Expand imprecision_data for ce_data //
    /////////////////////////////////////////
    
    IntegerVector ce_new_indices_imprecision(m);
    if(m != comparison_imprecision_sorted.size()){
      stop("Each comparison in imprecision_data must correspond to exactly one set of imprecision estimates. This is not the case in imprecision_data, so calculations are terminated");
    }
    
    // Correct placement of numbers after sorting comparison
    for(int i = 0; i < m; ++i){
      ce_new_indices_imprecision[i] = pb_new_indices_imprecision[i];
    }
    
    // New order based on comparison sorted
    NumericVector ce_new_order_CV_A(m);
    NumericVector ce_new_order_CV_A_lwr(m);
    NumericVector ce_new_order_CV_A_upr(m);
    NumericVector ce_new_order_CV_B(m);
    NumericVector ce_new_order_CV_B_lwr(m);
    NumericVector ce_new_order_CV_B_upr(m);
    NumericVector ce_new_order_lambda(m);
    NumericVector ce_new_order_lambda_lwr(m);
    NumericVector ce_new_order_lambda_upr(m);
    
    // Updated order for imprecision estimates
    for(int i = 0; i < m; ++i){
      ce_new_order_CV_A[i] = CV_A[ce_new_indices_imprecision[i]];
      ce_new_order_CV_A_lwr[i] = CV_A_lwr[ce_new_indices_imprecision[i]];
      ce_new_order_CV_A_upr[i] = CV_A_upr[ce_new_indices_imprecision[i]];
      ce_new_order_CV_B[i] = CV_B[ce_new_indices_imprecision[i]];
      ce_new_order_CV_B_lwr[i] = CV_B_lwr[ce_new_indices_imprecision[i]];
      ce_new_order_CV_B_upr[i] = CV_B_upr[ce_new_indices_imprecision[i]];
      ce_new_order_lambda[i] = lambda[ce_new_indices_imprecision[i]];
      ce_new_order_lambda_lwr[i] = lambda_lwr[ce_new_indices_imprecision[i]];
      ce_new_order_lambda_upr[i] = lambda_upr[ce_new_indices_imprecision[i]];  
    }
    
    // expansion from m to M
    CharacterVector ce_new_comparison_imprecision(M);
    NumericVector ce_new_CV_A(M);
    NumericVector ce_new_CV_A_lwr(M);
    NumericVector ce_new_CV_A_upr(M);
    NumericVector ce_new_CV_B(M);
    NumericVector ce_new_CV_B_lwr(M);
    NumericVector ce_new_CV_B_upr(M);
    NumericVector ce_new_lambda(M);
    NumericVector ce_new_lambda_lwr(M);
    NumericVector ce_new_lambda_upr(M);
    
    // reset counter
    c = 0;
    
    // Updated imprecision estimates based on the order found in comparison
    for(int i = 0; i < m; ++i){
      for(int j = 0; j < m_comparison[i]; ++j){
        ce_new_comparison_imprecision[c] = comparison_imprecision_sorted[i];
        ce_new_CV_A[c] = ce_new_order_CV_A[i];
        ce_new_CV_A_lwr[c] = ce_new_order_CV_A_lwr[i];
        ce_new_CV_A_upr[c] = ce_new_order_CV_A_upr[i];
        ce_new_CV_B[c] = ce_new_order_CV_B[i];
        ce_new_CV_B_lwr[c] = ce_new_order_CV_B_lwr[i];
        ce_new_CV_B_upr[c] = ce_new_order_CV_B_upr[i];
        ce_new_lambda[c] = ce_new_order_lambda[i];
        ce_new_lambda_lwr[c] = ce_new_order_lambda_lwr[i];
        ce_new_lambda_upr[c] = ce_new_order_lambda_upr[i];
        ++c;
      }
    }
    
    List out_1 = List::create(Named("comparison") = comparison_sorted,
                              Named("zeta") = round(new_zeta_pe, rounding),
                              Named("zeta_ci_lwr") = round(new_zeta_lwr, rounding),
                              Named("zeta_ci_upr") = round(new_zeta_upr, rounding),
                              Named("zeta_upper") = round(new_zeta_crit, rounding),
                              Named("dins_conclusion") = new_zeta_conc,
                              Named("predictor") = round(new_predictor, rounding),
                              Named("prediction") = round(new_predictor, rounding),
                              Named("pi_lwr") = round(new_lwr, rounding),
                              Named("pi_upr") = round(new_upr, rounding),
                              Named("CV_A") = round(new_CV_A, rounding),
                              Named("CV_A_lwr") = round(new_CV_A_lwr, rounding),
                              Named("CV_A_upr") = round(new_CV_A_upr, rounding),
                              Named("CV_B") = round(new_CV_B, rounding),
                              Named("CV_B_lwr") = round(new_CV_B_lwr, rounding),
                              Named("CV_B_upr") = round(new_CV_B_upr, rounding),
                              Named("lambda") = round(new_lambda, rounding),
                              Named("lambda_lwr") = round(new_lambda_lwr, rounding),
                              Named("lambda_upr") = round(new_lambda_upr, rounding));
    
    List out_2 = List::create(Named("comparison") = ce_new_comparison_imprecision,
                              Named("SampleID") = new_SampleID,
                              Named("zeta") = round(ce_new_zeta_pe, rounding),
                              Named("zeta_ci_lwr") = round(ce_new_zeta_lwr, rounding),
                              Named("zeta_ci_upr") = round(ce_new_zeta_upr, rounding),
                              Named("zeta_upper") = round(ce_new_zeta_crit, rounding),
                              Named("dins_conclusion") = ce_new_zeta_conc,
                              Named("MS_B") = round(new_MS_B, rounding),
                              Named("MS_A") = round(new_MS_A, rounding),
                              Named("prediction") = round(new_ce_prediction, rounding),
                              Named("pi_lwr") = round(new_ce_lwr, rounding),
                              Named("pi_upr") = round(new_ce_upr, rounding),
                              Named("pi_inside") = new_ce_inside,
                              Named("CV_A") = round(ce_new_CV_A, rounding),
                              Named("CV_A_lwr") = round(ce_new_CV_A_lwr, rounding),
                              Named("CV_A_upr") = round(ce_new_CV_A_upr, rounding),
                              Named("CV_B") = round(CV_B, rounding),
                              Named("CV_B_lwr") = round(ce_new_CV_B_lwr, rounding),
                              Named("CV_B_upr") = round(ce_new_CV_B_upr, rounding),
                              Named("lambda") = round(ce_new_lambda, rounding));
    
    List out = List::create(Named("merged_pb_data") = out_1, Named("merged_ce_data") = out_2);
    return out;
    
  }
  
  List out_1 = List::create(Named("comparison") = comparison_sorted,
                            Named("zeta") = round(new_zeta_pe, rounding),
                            Named("zeta_ci_lwr") = round(new_zeta_lwr, rounding),
                            Named("zeta_ci_upr") = round(new_zeta_upr, rounding),
                            Named("zeta_upper") = round(new_zeta_crit, rounding),
                            Named("dins_conclusion") = new_zeta_conc,
                            Named("predictor") = round(new_predictor, rounding),
                            Named("prediction") = round(new_predictor, rounding),
                            Named("pi_lwr") = round(new_lwr, rounding),
                            Named("pi_upr") = round(new_upr, rounding));
  
  List out_2 = List::create(Named("comparison") = ce_new_comparison_zeta,
                            Named("SampleID") = new_SampleID,
                            Named("zeta") = round(ce_new_zeta_pe, rounding),
                            Named("zeta_ci_lwr") = round(ce_new_zeta_lwr, rounding),
                            Named("zeta_ci_upr") = round(ce_new_zeta_upr, rounding),
                            Named("zeta_upper") = round(ce_new_zeta_crit, rounding),
                            Named("dins_conclusion") = ce_new_zeta_conc,
                            Named("MS_B") = round(new_MS_B, rounding),
                            Named("MS_A") = round(new_MS_A, rounding),
                            Named("prediction") = round(new_ce_prediction, rounding),
                            Named("pi_lwr") = round(new_ce_lwr, rounding),
                            Named("pi_upr") = round(new_ce_upr, rounding),
                            Named("pi_inside") = new_ce_inside);
  
  List out = List::create(Named("merged_pb_data") = out_1, Named("merged_ce_data") = out_2);
  return out;
}
