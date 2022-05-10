#include <Rcpp.h>
using namespace Rcpp;

//' Estimate pointwise prediction intervals based on Deming regression
//'
//' @title Estimation of prediction intervals or prediction band based on Deming regression
//' @name PredictDeming
//' 
//' @param data A data table or data frame converted to a list so that each column is a unique list element. Must contain the following columns: MP_A and MP_B
//' @param NewData A numeric vector of which x values are used to predict y
//' @param Precision A list of the Precision estimates based on \code{CharacterEstimatePrecision()}
//' @param MaxR The number of replicates for the majority of clinical samples. For example, if 18 of 25 CSs are meausured in triplicate, and the remaining seven is measured in duplicate, we use MaxR = 3
//' @param level the confidence level of the pointwise prediction intervals constructed. If the prediction intervals are used to evaluate commutability of EQAMs, you must remember to adjusts for simultanuous testing, e.g., Bonferroni correction
//' @param type which type of estimation approach should be used to construct the pointwise prediction intervals? \code{type = 0} correspond to J. Gillard and W. Fullers' approach. \code{type = 1} corresponds to CLSI's approach (elaborated upon in EP14). The first-mentioned is the default choice
//' @param silence an integer that controls progress reports. Set to 0 for these progress reports to appear. For additional debugging information, set this to \code{-1}. Default is \code{1}, which means that only the output will be returned, and no messages
//' 
//' @description Estimation of prediction intervals based on the Deming model is relevant in e.g., commutability evaluation of EQAMs. 
//' 
//' @return A list with four elements, which are, nx (inputted new data), fit (predicted y values based on new data), lwr (lower prediction interval bound), upr (upper prediction interval bound). Use \code{setDT()} from data.table to merge the list elements into a data table
//'

// [[Rcpp::export]]
List PredictDeming(List data, NumericVector NewData, List Precision, int MaxR, float level = 0.99, int type = 0, int silence = 1){
  
  float lambda = Precision["lambda"];
  float Var_B = Precision["Var_B"];
  NumericVector x = data["MP_B"];
  NumericVector y = data["MP_A"];
  float R = MaxR;
  int n = x.size();
  float N = x.size();
  int m = y.size();
  int j = NewData.size();
  NumericVector ypred(j);
  if(n != m){
    stop("Error: x and y is not of same length... please fix");
  }
  float mean_x = mean(x);
  float mean_y = mean(y);
  float gmxx = 0;
  float gmyy = 0;
  float gmxy = 0;
  for(int i = 0; i < n; ++i){
    gmxx += (x[i] - mean_x) * (x[i] - mean_x) / (N - 1);
    gmyy += (y[i] - mean_y) * (y[i] - mean_y) / (N - 1);
    gmxy += (x[i] - mean_x) * (y[i] - mean_y) / (N - 1);
  }
  float cmxx = gmxx * (N - 1) / N;
  float cmyy = gmyy * (N - 1) / N;
  float cmxy = gmxy * (N - 1) / N;
  // type = 0: J. Gillard and W. Fuller
  if(type == 0){
    float p1 = gmyy - lambda * gmxx;
    float p2 = sqrt(pow(gmyy - lambda * gmxx, 2) + 4.0 * lambda * pow(gmxy, 2));
    float p3 = 2.0 * gmxy;
    float b1 = (p1 + p2) / p3;
    float b0 = mean_y - b1 * mean_x;
    float AsyVarb1 = (pow(b1, 2) / (N * pow(gmxy, 2))) * ((gmxx * gmyy) - pow(gmxy, 2));
    NumericVector xpred(j);
    float mu_hat = 0;
    float ShAdjusted = (gmyy + (lambda * gmxx) - p2) / (2 * lambda);
    float SvAdjusted = lambda * ShAdjusted;
    if(silence == 0){
      Rcout << "You have choosed to estimate PB using standard approach (J. Gillard and Fuller)" << "\n";
      Rcout << "The base value of sigma_h squared is" << (gmyy + (lambda * gmxx) - p2) / (2 * lambda) << "\n";
      Rcout << "The value of sigma_h squared is : " << ShAdjusted << "\n";
      Rcout << "The value of sigma_v squared is : " << SvAdjusted << "\n";
      Rcout << "The value of var(b1) is : " << AsyVarb1 << "\n";
      Rcout << "The value of b1 is : " << b1 << "\n";
      Rcout << "The value of b0 is : " << b0 << "\n";
      Rcout << "The value of mxx is : " << gmxx << "\n";
      Rcout << "The value of myy is : " << gmyy << "\n";
      Rcout << "The value of mxy is : " << gmxy << "\n";
      Rcout << "The value of R / 3 is : " << MaxR / 3.0 << "\n";
      Rcout << "F-G uses " << n - 2 << " degrees of freedom" << "\n";
    }
    float tquant = R::qt((1 - level)/2, n - 2, 0, 0);
    for(int i = 0; i < j; ++i){
      ypred[i] = b0 + b1 * NewData[i];
      xpred[i] = ((lambda / (lambda + pow(b1, 2))) * NewData[i]) + (b1 / (lambda + pow(b1, 2))) * (ypred[i] - b0);
    }
    for(int i = 0; i < n; ++i){
      float mu_hat_part_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      mu_hat += mu_hat_part_i / N;
    }
    if(silence == 0){
      Rcout << "In the F-G approach we estimate true concentration levels:" << "\n";
      Rcout << "The estimated relationship between true and observed is" << "\n";
      Rcout << "lambda / (lambda + b1^2) x_i + b1 / (lambda + b1^2) * (y_i - b0), which is" << "\n";
      Rcout << "In this particular case (.obs 1) we have that estimated true concentration is :" << xpred[0] << "\n";
      Rcout << "In this particular case (.obs 1) we had the observed measurement :" << NewData[0] << "\n";
      Rcout << "Before dividing by R : " << (1.0 + 1.0 / N) * (pow(b1, 2) * ShAdjusted + SvAdjusted) << "\n";
      Rcout << "After dividing by R : " << (1.0 + 1.0 / N) * (pow(b1, 2) * ShAdjusted + SvAdjusted) / R << "\n";
    }
    NumericVector lwr(j);
    NumericVector upr(j);
    for(int i = 0; i < j; ++i){
      float PredVariance = AsyVarb1 * pow((xpred[i] - mu_hat), 2) + (AsyVarb1 * ShAdjusted / R) + (1.0 + 1.0 / N) * (b1 * b1 * ShAdjusted + SvAdjusted) / R;
      float PredError = sqrt(PredVariance);
      lwr[i] = ypred[i] - tquant * PredError;
      if(lwr[i] < 0){
        lwr[i] = 0;
      }
      upr[i] = ypred[i] + tquant * PredError;
    }
    
    List out = List::create(Named("nx") = NewData, Named("fit") = round(ypred, 5), Named("lwr") = round(lwr, 5), Named("upr") = round(upr, 5));
    return out;
  }
  // type = 1: CLSI (EP14)
  else{
    float p1 = cmyy - lambda * cmxx;
    float p2 = sqrt((cmyy - lambda * cmxx) * (cmyy - lambda * cmxx) + 4 * lambda * cmxy * cmxy);
    float p3 = 2 * cmxy;
    float b1 = (p1 + p2) / p3;
    float b0 = mean_y - b1 * mean_x;
    float AsyVarb1 = ((b1 * b1) / (n * cmxy * cmxy)) * ((cmxx * cmyy) - (cmxy * cmxy));
    float ShAdjusted = Var_B;
    float SvAdjusted = lambda * ShAdjusted;
    if(silence == 0){
      Rcout << "You have choosed to estimate PB using EP14's approach (CLSI)" << "\n";
      Rcout << "The base value of sigma_h squared is" << ShAdjusted << "\n";
      Rcout << "The value of sigma_h squared is : " << ShAdjusted << "\n";
      Rcout << "The value of sigma_v squared is : " << SvAdjusted << "\n";
      Rcout << "The value of var(b1) is : " << AsyVarb1 << "\n";
      Rcout << "The value of b1 is : " << b1 << "\n";
      Rcout << "The value of b0 is : " << b0 << "\n";
      Rcout << "The value of mxx is : " << cmxx << "\n";
      Rcout << "The value of myy is : " << cmyy << "\n";
      Rcout << "The value of mxy is : " << cmxy << "\n";
      Rcout << "The value of R / 3 is : " << float(MaxR) / 3.0 << "\n";
      Rcout << "CLSI uses " << n * (MaxR - 1) << "degrees of freedom" << "\n";
    }
    
    float tquant = R::qt((1 - level)/2, n * (MaxR - 1), 0, 0);
    for(int i = 0; i < j; ++i){
      ypred[i] = b0 + b1 * NewData[i];
    }
    NumericVector lwr(j);
    NumericVector upr(j);
    for(int i = 0; i < j; ++i){
      float PredVariance = AsyVarb1 * ((NewData[i] - mean_x) * (NewData[i] - mean_x)) + (1.0 + 1.0 / N) * (b1 * b1 * ShAdjusted + SvAdjusted) / R;
      float PredError = sqrt(PredVariance);
      lwr[i] = ypred[i] - tquant * PredError;
      upr[i] = ypred[i] + tquant * PredError;
    }
    float clsi_ignores = AsyVarb1 * ShAdjusted;
    if(silence == 0){
      Rcout << "CLSI method ignores this much prediction variance : " << clsi_ignores << "\n";
      Rcout << "Before dividing by R" << (1.0 + 1.0 / N) * (b1 * b1 * ShAdjusted + SvAdjusted) << "\n";
      Rcout << "After dividing by R" << (1.0 + 1.0 / N) * (b1 * b1 * ShAdjusted + SvAdjusted) / R << "\n";
    }
    List out = List::create(Named("nx") = NewData, Named("fit") = round(ypred, 5), Named("lwr") = round(lwr, 5), Named("upr") = round(upr, 5));
    return out;
  }
}