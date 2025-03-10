#include <Rcpp.h>
using namespace Rcpp;

// (*) Helper function for obtaining unique elements maintaining their original order
CharacterVector unique_preserve_order_csr(CharacterVector x) {
  std::unordered_set<std::string> seen;
  CharacterVector result;
  for(int i = 0; i < x.length(); i++) {
    std::string current = as<std::string>(x[i]);
    if(seen.find(current) == seen.end()) {
      seen.insert(current);
      result.push_back(x[i]);
    }
  }
  return result;
}

//' @title Count the Number of Replicated for Each Sample
//' @name count_samplewise_replicates
//'
//' @param data A \code{list} or \code{data.table} representing CS data or EQAM data
//'        with the \code{list} elements or \code{data.table} columns:
//'        \code{SampleID}, \code{MP_A}, and \code{MP_B}.
//' @param summary A \code{character} specifying the summary statistic of the sample-wise
//'        numbers of replicates. Default is \code{'mode'}. Possible summary statistics include:
//'        \itemize{
//'           \item \code{none}: Returns an \code{integer} vector containing every
//'           sample-wise number of replicatees for \code{data}.
//'           \item \code{mode}: Returns the sample mode.
//'           \item \code{median}: Returns the sample median.
//'           \item \code{mean}: Returns the sample mean.
//'           \item \code{ceiling}: Returns the rounded up (to nearest integer) sample mean.
//'           \item \code{floor}: Returns the rounded down (to nearest integer) sample mean.
//'           \item \code{round}: Returns the rounded (to nearest integer) sample mean.
//'        }
//' @param invalid_NA A \code{logical} that determines the behavior of the function in response
//'        to invalid \code{data} or \code{summary} input.
//'        If \code{TRUE}, the function will return an \code{NA}-value when encountering
//'        invalid input or computation errors, rather than throwing an error.
//'        While generally not recommended due to potential masking of issues,
//'        this may be useful in certain scenarios where uninterrupted execution is desired.
//' @param silence An \code{integer} that dictates the verbosity level of the console output.
//'        If \code{silence} is set to a value less than 1, various verbose output,
//'        including debugging reports, will be displayed depending on the specific
//'        value of \code{silence}.
//' 
//' @description
//' This function counts the number of replicated measurements done on each sample
//' within an IVD-MD comparison for either clinical sample (CS) data or
//' external quality asessment material (EQAM) data. Alternatively, a summary
//' statistic can be returned to summarize the numbers of replicates across all
//' samples in the \code{data}.
//'
//' @details
//' The function \code{predict_eqa()} hinges on the count of replicates
//' performed for its \code{data} parameter. When the \code{method} is either
//' \code{fg} or \code{ols}, it is essential to tally the replicates for both
//' \code{data} and \code{new_data} inputs. This handy function offers a
//' streamlined solution. It allows for the direct usage of \code{predict_eqa()}
//' without requiring manual counting of the sample-wise replicates,
//' thereby enhancing efficiency and ease of use.
//' 
//' Note: If the count of replicates for a particular sample is 0, it will be
//' removed from the output. 
//'
//' @return
//' Returns a \code{list} that contains a single element, \code{R_i}.
//' The type of \code{R_i} depends on the \code{summary} parameter.
//' It is typically an single integer value.
//' However, if \code{summary} is set to \code{'none'}, \code{R_i} is an
//' integer vector because it list the number of replicates for each sample.
//' When \code{summary} is set to \code{'mean'},
//' \code{R_i} is returned as a \code{double}.
//'
//' @examples \dontrun{
//' library(fasteqa)
//' # Simulation parameters for clinical sample data
//' cs_parameters <- list(n = 25, R = 3,
//'                       cvx = 0.01, cvy = 0.015,
//'                       cil = 10, ciu = 70)
//' # Use the simulation parameters to simulate toy clinical sample data                      
//' cs_data <- simulate_eqa_data(cs_parameters)
//' # Calculate mode of sample-wise number of replicates
//' mode_R <- count_samplewise_replicates(cs_data, summary = 'mode')$R_i
//' # Count sample-wise number of replicates
//' samplewise_R <- count_samplewise_replicates(cs_data, summary = 'none')$R_i                                    
//' }


// [[Rcpp::export]]
List count_samplewise_replicates(List data, String summary = "mode", bool invalid_NA = true, int silence = 1) {
  
  if (!(data.containsElementNamed("SampleID") && data.containsElementNamed("ReplicateID"))) {
    if(!invalid_NA){
      stop("data does not contain at least one of 'SampleID' or 'ReplicateID'");  
    }
    else{
      return List::create(Named("R_i") = NA_INTEGER);  
    }
  }
  
  if (!(data.containsElementNamed("MP_B") && data.containsElementNamed("MP_A"))) {
    if(!invalid_NA){
      stop("data does not contain at least one of 'MP_B' and 'MP_A'");
    }
    else{
      return List::create(Named("R_i") = NA_INTEGER);
    }
  }
  
  // Extract relevant information from 'data'
  NumericVector x_raw = data["MP_B"];
  NumericVector y_raw = data["MP_A"];
  CharacterVector SampleID_raw = data["SampleID"];
  CharacterVector ReplicateID_raw = data["ReplicateID"];
  CharacterVector SampleID;
  CharacterVector ReplicateID;
  
  // Remove SampleID + ReplicateID corresponding with MP_A / MP_B NA-values
  for(int i = 0; i < x_raw.size(); ++i){
    bool is_na_x = ISNAN(x_raw[i]);
    bool is_na_y = ISNAN(y_raw[i]);
    if(is_na_x || is_na_y){
      continue;
    }
    SampleID.push_back(SampleID_raw[i]);
    ReplicateID.push_back(ReplicateID_raw[i]);
  }
  // Number of measurements in 'data' after NA values are removed
  int N = SampleID.size();
  // Number of unique samples in 'data' after NA values are removed
  CharacterVector unique_SampleID = unique_preserve_order_csr(SampleID);
  int n = unique_SampleID.size();
  IntegerVector R_i(n);
  
  if(n > N){
    stop("n[%d] (number of unique samples) is larger than N[%d] (number of unique measurements)!", n, N);
  }
  else if(n == N){
    for(int i = 0; i < n; ++i){
      R_i[i] = 1;
    }
    if(summary == "none"){
      List out = List::create(Named("R_i") = R_i);
      return out;
    }
    else if(summary == "mode" || summary == "median"){
      List out = List::create(Named("R_i") = 1);
      return out;
    }
    else if(summary == "ceiling" || summary == "floor" || summary == "round"){
      List out = List::create(Named("R_i") = 1);
      return out;
    }
    else if(summary == "mean"){
      List out = List::create(Named("R_i") = 1);
      return out;
    }
    else{
      if(!invalid_NA){
        std::string input_summary = summary.get_cstring();
        stop("summary = '%s' is not valid. Valid entries are: %s %s 'none', 'mode', 'median', 'floor', 'ceiling', 'round' and 'mean'.", input_summary.c_str(), "\n", " *");      
      }
      else{
        return List::create(Named("R_i") = NA_INTEGER);
      }      
    }
    
    
  }
  
  for (int i = 0; i < n; ++i) {
    std::string sample = as<std::string>(unique_SampleID[i]);
    int count = 0;
    for (int j = 0; j < N; ++j) {
      if (sample == as<std::string>(SampleID[j])) {
        count++;
      }
    }
    R_i[i] = count;
  }
  
  int R_mode = R_i[0];
  
  if(summary == "mode"){
    std::map<int, int> counts;
    for(int i = 0; i < n; ++i) {
      counts[R_i[i]]++;
    }
    
    int max_count = 0;
    for(auto &p : counts) {
      if(p.second > max_count) {
        max_count = p.second;
        R_mode = p.first;
      }
    }
  }
  
  if(summary == "none"){
    List out = List::create(Named("R_i") = R_i);
    return out;
  }
  
  else if(summary == "mode" || summary == "median"){
    if(summary == "mode"){
      List out = List::create(Named("R_i") = R_mode);
      return out;
    }
    else{
      List out = List::create(Named("R_i") = median(R_i));
      return out;
    }
    return List::create(Named("R_i") = NA_INTEGER);
  }
  else if(summary == "ceiling" || summary == "floor" || summary == "round" || summary == "mean"){
    double R_mean = mean(R_i);
    if(summary == "ceiling"){
      List out = List::create(Named("R_i") = std::ceil(R_mean));
      return out;
    }
    else if(summary == "floor"){
      List out = List::create(Named("R_i") = std::floor(R_mean));
      return out;
    }
    else if(summary == "round"){
      List out = List::create(Named("R_i") = std::round(R_mean));
      return out;
    }
    else{
      List out = List::create(Named("R_i") = R_mean);
      return out;
    }
    return List::create(Named("R_i") = NA_INTEGER);
  }
  else{
    if(!invalid_NA){
      std::string input_summary = summary.get_cstring();
      stop("summary = '%s' is not valid. Valid entries are: %s %s 'none', 'mode', 'median', 'floor', 'ceiling', 'round' and 'mean'.", input_summary.c_str(), "\n", " *");      
    }
    else{
      return List::create(Named("R_i") = NA_INTEGER);
    }     
  }
  
  return List::create(Named("R_i") = NA_INTEGER);
}
