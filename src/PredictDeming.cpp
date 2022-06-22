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
//' @param MaxR The number of replicates for the majority of clinical samples. For example, if 18 of 25 CSs are meausured in triplicate, and the remaining seven is measured in duplicate, we use \code{MaxR = 3}
//' @param level the confidence level of the pointwise prediction intervals constructed. If the prediction intervals are used to evaluate commutability of EQAMs, you must remember to adjusts for simultanuous testing, e.g., Bonferroni correction
//' @param type which type of estimation approach should be used to construct the pointwise prediction intervals? \code{type = 0} correspond to J. Gillard and W. Fullers' approach. \code{type = 1} corresponds to CLSI's approach (elaborated upon in EP14). The first-mentioned is the default choice
//' @param rounding how many decimal places should be use in the estimates? Default is 3
//' @param Should estimated true concentration levels be part of the output? Default is 0 which signify no. This will only return something if \code{type = 0}
//' @param silence an integer that controls progress reports. Set to 0 for these progress reports to appear. For additional debugging information, set this to \code{-1}. Default is \code{1}, which means that only the output will be returned, and no messages
//' 
//' @description Estimation of prediction intervals based on the Deming model is relevant in e.g., commutability evaluation of EQAMs. 
//' 
//' @return A list with four elements, which are, nx (inputted new data), fit (predicted y values based on new data), lwr (lower prediction interval bound), upr (upper prediction interval bound). Use \code{setDT()} from data.table to merge the list elements into a data table
//'

// [[Rcpp::export]]
List PredictDeming(List data, NumericVector NewData, List Precision, int MaxR, float level = 0.99, int type = 0, int rounding = 3, int CalculateLatent = 0, int silence = 1){
  
  float lambda = Precision["lambda"];
  float Var_B = Precision["Var_B"];
  NumericVector x = data["MP_B"];
  NumericVector y = data["MP_A"];
  float R = MaxR;
  int n = x.size();
  int m = y.size();
  int j = NewData.size();
  NumericVector ypred(j);
  if(n != m){
    stop("Error: ensure x and y is of same length");
  }
  float mean_x = mean(x);
  float mean_y = mean(y);
  float gmxx = 0;
  float gmyy = 0;
  float gmxy = 0;
  for(int i = 0; i < n; ++i){
    gmxx += pow(x[i] - mean_x, 2) / (n - 1);
    gmyy += pow(y[i] - mean_y, 2) / (n - 1);
    gmxy += (x[i] - mean_x) * (y[i] - mean_y) / (n - 1);
  }
  float cmxx = gmxx * (n - 1) / n;
  float cmyy = gmyy * (n - 1) / n;
  float cmxy = gmxy * (n - 1) / n;
  // type = 0: J. Gillard and W. Fuller
  if(type == 0){
    NumericVector latent (n);
    NumericVector xpred(j);
    float p1 = gmyy - lambda * gmxx;
    float p2 = sqrt(pow(gmyy - lambda * gmxx, 2) + 4 * lambda * pow(gmxy, 2));
    float p3 = 2 * gmxy;
    float b1 = (p1 + p2) / p3;
    float b0 = mean_y - b1 * mean_x;
    float varb1 = (pow(b1, 2) / (n * pow(gmxy, 2))) * ((gmxx * gmyy) - pow(gmxy, 2));
    float mu_hat = 0;
    float hvar = R * (gmyy + (lambda * gmxx) - p2) / (2 * lambda);
    float vvar = lambda * hvar;
    // float varb1 = (n - 1) * (lambda + pow(b1, 2)) * hvar / (n - 2);
    if(silence == 0){
      Rcout << "You have choosed to estimate PB using standard approach (J. Gillard and Fuller)" << "\n";
      Rcout << "The base value of sigma_h squared is : " << hvar << "\n";
      Rcout << "The base value of sigma_v squared is : " << vvar << "\n";
      Rcout << "The value of lambda is : " << lambda << "\n";
      Rcout << "The value of sigma_h squared is : " << hvar / R << "\n";
      Rcout << "The value of sigma_v squared is : " << vvar / R << "\n";
      Rcout << "The value of var(b1) is : " << varb1 << "\n";
      Rcout << "The value of b1 is : " << b1 << "\n";
      Rcout << "The value of b0 is : " << b0 << "\n";
    }
    if(silence == -1){
      Rcout << "The value of mxx is : " << gmxx << "\n";
      Rcout << "The value of myy is : " << gmyy << "\n";
      Rcout << "The value of mxy is : " << gmxy << "\n";
      Rcout << "F-G uses " << n - 2 << " degrees of freedom" << "\n";  
    }
    float tquant = R::qt((1 - level)/2, n - 2, 0, 0);
    for(int i = 0; i < n; ++i){
      float latent_i = ((lambda / (lambda + pow(b1, 2))) * x[i]) + (b1 / (lambda + pow(b1, 2))) * (y[i] - b0);
      mu_hat += latent_i / n;
      latent[i] = latent_i;
    }
    for(int i = 0; i < j; ++i){
      ypred[i] = b0 + b1 * NewData[i];
      xpred[i] = mu_hat * (1 - (gmxy / (b1 * gmxx))) + NewData[i] * (gmxy / (b1 * gmxx));
    }
    if(silence == -1){
      Rcout << "In the F-G approach we estimate true concentration levels:" << "\n";
      Rcout << "Mean of measurements for new samples values are: " << NewData << "\n";
      Rcout << "Predicted true concentration values are: " << xpred << "\n";
    }
    NumericVector lwr(j);
    NumericVector upr(j);
    for(int i = 0; i < j; ++i){
      float PredVariance = varb1 * pow((xpred[i] - mu_hat), 2) + (varb1 * hvar / R) + (1 + (1 / n)) * (pow(b1, 2) * hvar + vvar) / R;
      float PredError = sqrt(PredVariance);
      lwr[i] = ypred[i] - tquant * PredError;
      if(lwr[i] < 0){
        lwr[i] = 0;
      }
      upr[i] = ypred[i] + tquant * PredError;
    }
    List outwol = List::create(Named("nx") = NewData, Named("fit") = round(ypred, rounding), Named("lwr") = round(lwr, rounding), Named("upr") = round(upr, rounding));
    List outwl = List::create(Named("nx") = NewData, Named("fit") = round(ypred, rounding), Named("lwr") = round(lwr, rounding), Named("upr") = round(upr, rounding), Named("latent train") = latent, Named("latent test") = xpred);
    if(CalculateLatent == 0){
      return outwol;
    }
    return outwl;
    
    
  }
  // type = 1: CLSI (EP14)
  else{
    float p1 = cmyy - lambda * cmxx;
    float p2 = sqrt((cmyy - lambda * cmxx) * (cmyy - lambda * cmxx) + 4 * lambda * cmxy * cmxy);
    float p3 = 2 * cmxy;
    float b1 = (p1 + p2) / p3;
    float b0 = mean_y - b1 * mean_x;
    float varb1 = (pow(b1, 2) / (n * pow(cmxy, 2))) * (cmxx * cmyy - pow(cmxy, 2));
    float hvar = Var_B;
    float vvar = lambda * hvar;
    if(silence == 0){
      Rcout << "You have choosed to estimate PB using EP14's approach (CLSI)" << "\n";
      Rcout << "The base value of sigma_h squared is : " << hvar << "\n";
      Rcout << "The base value of sigma_v squared is : " << vvar << "\n";
      Rcout << "The value of lambda is : " << lambda << "\n";
      Rcout << "The value of sigma_h squared is : " << hvar / R << "\n";
      Rcout << "The value of sigma_v squared is : " << vvar / R << "\n";
      Rcout << "The value of var(b1) is : " << varb1 << "\n";
      Rcout << "The value of b1 is : " << b1 << "\n";
      Rcout << "The value of b0 is : " << b0 << "\n";
    }
    if(silence == -1){
      Rcout << "The value of mxx is : " << cmxx << "\n";
      Rcout << "The value of myy is : " << cmyy << "\n";
      Rcout << "The value of mxy is : " << cmxy << "\n";
      Rcout << "CLSI uses " << n * (MaxR - 1) << "degrees of freedom" << "\n";
    }
    float tquant = R::qt((1 - level)/2, n * (MaxR - 1), 0, 0);
    for(int i = 0; i < j; ++i){
      ypred[i] = b0 + b1 * NewData[i];
    }
    NumericVector lwr(j);
    NumericVector upr(j);
    for(int i = 0; i < j; ++i){
      float PredVariance = varb1 * pow(NewData[i] - mean_x, 2) + (1 + (1 / n)) * (pow(b1, 2) * hvar + vvar) / R;
      float PredError = sqrt(PredVariance);
      lwr[i] = ypred[i] - tquant * PredError;
      upr[i] = ypred[i] + tquant * PredError;
    }
    float clsi_ignores = varb1 * hvar;
    if(silence == -1){
      Rcout << "CLSI method ignores this much prediction variance : " << clsi_ignores << "\n";
    }
    List out = List::create(Named("nx") = NewData, Named("fit") = round(ypred, rounding), Named("lwr") = round(lwr, rounding), Named("upr") = round(upr, rounding));
    return out;
  }
}
