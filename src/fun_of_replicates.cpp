#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// (*) Helper function to calculate the specified statistic
double calculate_stat(const std::vector<double>& values, const std::string& fun) {
  if (values.empty()) {
    return NA_REAL;
  }

  double result = 0.0;
  size_t n = values.size();

  if (fun == "mean") {
    result = std::accumulate(values.begin(), values.end(), 0.0) / n;
  } else if (fun == "var") {
    if (std::adjacent_find(values.begin(), values.end(),
                           std::not_equal_to<>()) == values.end()) {
      result = NA_REAL;
    }
    else{
      double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
      double sqSum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
      result = (sqSum - n * mean * mean) / (n - 1);
    }

  } else if (fun == "median") {
    std::vector<double> sortedValues = values;
    std::sort(sortedValues.begin(), sortedValues.end());
    result = (n % 2 == 0) ? (sortedValues[n/2 - 1] + sortedValues[n/2]) / 2 : sortedValues[n/2];
  } else if (fun == "sd") {
    if (std::adjacent_find(values.begin(), values.end(),
                           std::not_equal_to<>()) == values.end()) {
      result = NA_REAL;
    }
    else{
      double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
      double sqSum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
      result = std::sqrt((sqSum - n * mean * mean) / (n - 1));
    }
  } else if (fun == "cv") {
    if (std::adjacent_find(values.begin(), values.end(),
                           std::not_equal_to<>()) == values.end()) {
      result = NA_REAL;
    }
    else{
      double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
      double sqSum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
      double stDev = std::sqrt((sqSum - n * mean * mean) / (n - 1));
      result = stDev / mean;
    }
  } else if (fun == "min") {
    result = *std::min_element(values.begin(), values.end());
  } else if (fun == "max") {
    result = *std::max_element(values.begin(), values.end());
  } else {
    Rcpp::stop("Invalid statistic type specified");
  }

  return result;
}

//' @title Estimate a Statistic for Each Sample
//' @name fun_of_replicates
//'
//' @param data A \code{data.table} or \code{list} object. Must contain
//'        \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
//' @param fun A \code{character} string. Which statistic is to be calculated for
//'        each SampleID. Possible choices include \code{mean}, \code{var} (variance),
//'        \code{sd} (standard deviaton), \code{cv} (coefficient of variation),
//'        \code{median}, \code{min} (minimum) and \code{max} (maximum).
//'
//' @description
//' For each sample, uses the replicated measurements to estimate
//' a given statistic.
//'
//' @details
//' This function handles \code{NA}-values automatically. If all replicated measurements
//' for a given sample are \code{NA}-values, the corresponding estimated statistic will also
//' be \code{NA}. 
//'
//' @return
//' A \code{list} with elements \code{SampleID}, \code{MP_A} and \code{MP_B}.
//' \code{MP_A} and \code{MP_B} contains the estimated sample-wise statistics.
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   test_data <- simulate_eqa_data(list(n = 25, R = 5, cvx = 0.01,
//'                                       cvy = 0.01, cil = 2, ciu = 10))
//'   print(as.data.frame(fun_of_replicates(test_data, "mean")))
//'   print(as.data.frame(fun_of_replicates(test_data, "var")))
//'   print(as.data.frame(fun_of_replicates(test_data, "cv")))
//' }
// [[Rcpp::export]]
List fun_of_replicates(List data, std::string fun = "mean") {
  CharacterVector sampleID = data["SampleID"];
  NumericVector MP_A = data["MP_A"];
  NumericVector MP_B = data["MP_B"];

  std::unordered_map<std::string, std::vector<double>> groupedDataA;
  std::unordered_map<std::string, std::vector<double>> groupedDataB;
  std::vector<std::string> uniqueIDsOrdered;

  int n = sampleID.size();

  for (int i = 0; i < n; ++i) {
    std::string id = as<std::string>(sampleID[i]);
    if (groupedDataA.find(id) == groupedDataA.end()) {
      uniqueIDsOrdered.push_back(id);
    }
    if (!NumericVector::is_na(MP_A[i])) {
      groupedDataA[id].push_back(MP_A[i]);
    }
    if (!NumericVector::is_na(MP_B[i])) {
      groupedDataB[id].push_back(MP_B[i]);
    }
  }

  CharacterVector uniqueIDs(uniqueIDsOrdered.size());
  NumericVector statA(uniqueIDsOrdered.size());
  NumericVector statB(uniqueIDsOrdered.size());

  for (size_t i = 0; i < uniqueIDsOrdered.size(); ++i) {
    const std::string& id = uniqueIDsOrdered[i];
    uniqueIDs[i] = id;

    statA[i] = calculate_stat(groupedDataA[id], fun);
    statB[i] = calculate_stat(groupedDataB[id], fun);
  }

  return List::create(
    Named("SampleID") = uniqueIDs,
    Named("MP_A") = statA,
    Named("MP_B") = statB
  );
}


