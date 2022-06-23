#include <Rcpp.h>
using namespace Rcpp;

//' Estimate pointwise prediction intervals based on OLS regression
//'
//' @title Estimation of prediction intervals or prediction band based on OLS regression
//' @name PredictOLS
//' 
//' @param data A data table or data frame converted to a list so that each column is a unique list element. Must contain the following columns: MP_A and MP_B
//' @param NewData A numeric vector of which x values are used to predict y
//' @param Precision A list of the Precision estimates based on \code{CharacterEstimatePrecision()}
//' @param MaxR The number of replicates for the majority of clinical samples. For example, if 18 of 25 CSs are meausured in triplicate, and the remaining seven is measured in duplicate, we use MaxR = 3
//' @param level the confidence level of the pointwise prediction intervals constructed. If the prediction intervals are used to evaluate commutability of EQAMs, you must remember to adjusts for simultanuous testing, e.g., Bonferroni correction
//' @param silence an integer that controls progress reports. Set to 0 for these progress reports to appear. For additional debugging information, set this to \code{-1}. Default is \code{1}, which means that only the output will be returned, and no messages
//' 
//' @description Estimation of prediction intervals based on the OLS regression model is relevant in e.g., commutability evaluation of EQAMs
//' 
//' @return A list with four elements, which are, nx (inputted new data), fit (predicted y values based on new data), lwr (lower prediction interval bound), upr (upper prediction interval bound). Use \code{setDT()} from data.table to merge the list elements into a data table


// [[Rcpp::export]]
List PredictOLS(List data, NumericVector NewData, List Precision, int MaxR, float level, int silence = 1) {
  float lambda = Precision["lambda"];
  if(lambda < 1){
    NumericVector x = data["MP_A"];
    NumericVector y = data["MP_B"];
    int n = x.size();
    int m = y.size();
    int j = NewData.size();
    if(n != m){
      stop("Error: x and y is not of same length, this should be impossible...?");
    }
    float sxx = 0;
    float sxy = 0;
    float mean_x = mean(x);
    float mean_y = mean(y);
    for(int i = 0; i < n; ++i){
      sxx += (x[i] - mean_x) * (x[i] - mean_x);
      sxy += (x[i] - mean_x) * (y[i] - mean_y);
    }
    float b1 = sxy / sxx;
    float b0 = mean_y - b1 * mean_x;
    float mse = 0;
    for(int i = 0; i < n; ++i){
      float yhat = b0 + b1 * x[i];
      mse += MaxR * (y[i] - yhat) * (y[i] - yhat) / (n - 2);
    }
    mse = sqrt(mse);
    float tquant = R::qt((1 - level)/2, n - 2, 0, 0);
    NumericVector ypred(j);
    NumericVector lwr(j);
    NumericVector upr(j);
    CharacterVector WhichResponse(j);
    for(int i = 0; i < j; ++i){
      ypred[i] = b0 + b1 * NewData[i];
      float PredError = mse * sqrt(1 + (1 / n) + (NewData[i] - mean_x) * (NewData[i] - mean_x) / sxx);
      lwr[i] = ypred[i] - tquant * PredError;
      upr[i] = ypred[i] + tquant * PredError;
      WhichResponse[i] = "x";
    }
    List out = List::create(Named("nx") = NewData, Named("fit") = round(ypred, 4), Named("lwr") = round(lwr, 4), Named("upr") = round(upr, 4), Named("WhichTarget") = WhichResponse);
    return out;
  }
  else{
    NumericVector x = data["MP_B"];
    NumericVector y = data["MP_A"];
    int n = x.size();
    int m = y.size();
    int j = NewData.size();
    if(n != m){
      stop("Error: x and y is not of same length... ehhmmm..how?");
    }
    float sxx = 0;
    float sxy = 0;
    float mean_x = mean(x);
    float mean_y = mean(y);
    for(int i = 0; i < n; ++i){
      sxx += (x[i] - mean_x) * (x[i] - mean_x);
      sxy += (x[i] - mean_x) * (y[i] - mean_y);
    }
    float b1 = sxy / sxx;
    float b0 = mean_y - b1 * mean_x;
    float mse = 0;
    for(int i = 0; i < n; ++i){
      float yhat = b0 + b1 * x[i];
      mse += (y[i] - yhat) * (y[i] - yhat) / (n - 2) / MaxR;
    }
    mse = sqrt(mse);
    float tquant = R::qt((1 - level)/2, n - 2, 0, 0);
    if(silence == 0){
      Rcout << "The value of rmse is : " << mse << "\n";
      Rcout << "The value of b1 is : " << b1 << "\n";
      Rcout << "The value of b0 is : " << b0 << "\n";
      Rcout << "The value of sxx is : " << sxx << "\n";
      Rcout << "The value of sxy is : " << sxy << "\n";
      Rcout << "OLS ignores variance caused measurement error in x," << "\n";
      Rcout << "so its PB will usually be narrower than Deming!" << "\n"; 
    }
    NumericVector ypred(j);
    NumericVector lwr(j);
    NumericVector upr(j);
    CharacterVector WhichResponse(j);
    for(int i = 0; i < j; ++i){
      ypred[i] = b0 + b1 * NewData[i];
      float PredError = mse * sqrt(1 + (1 / n) + (NewData[i] - mean_x) * (NewData[i] - mean_x) / sxx);
      lwr[i] = ypred[i] - tquant * PredError;
      if(lwr[i] < 0){
        lwr[i] = 0;
      }
      upr[i] = ypred[i] + tquant * PredError;
      WhichResponse[i] = "y";
    }
    List out = List::create(Named("nx") = NewData, Named("fit") = round(ypred, 4), Named("lwr") = round(lwr, 4), Named("upr") = round(upr, 4), Named("WhichTarget") = WhichResponse);
    return out;
  }
}
