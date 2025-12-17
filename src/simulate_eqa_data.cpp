#include <Rcpp.h>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm> // for std::max

using namespace Rcpp;

//' @title Simulation of External Quality Assessment (EQA) Data (Deprecated)
//' @name simulate_eqa_data
//' 
//' @param parameters A \code{list} of parameters used to simulate the EQA data.
//'        You must at least specify one parameter for this function to run.
//'        Except that one mandatory parameter, you may optionally choose the
//'        remaining of the parameters. These are the optimal parameters that you may include into the list:
//' \itemize{
//'   \item \code{n: } The number of clinical samples
//'   \item \code{R: } The number of replicates on each sample
//'   \item \code{cvx: } The analytical CV of x measurements
//'   \item \code{cvy: } The analytical CV of y measurements
//'   \item \code{cil: } The lower range of the concentration interval
//'   \item \code{ciu: } The upper range of the concentration interval
//'   \item \code{dist: } The distribution to simulate latent variables from.
//'                       Possbile choices include 'unif' 
//'                       (uniform distribution, default),
//'                       'norm' (normal distribution) ,
//'                       'lst' (location-scale student t-distribution),
//'                       'lnorm' (log-normal distribution)
//'   \item \code{df_tau: } The degrees of freedom for the 'lst' distribution
//'                         if the distribution of latent variables are
//'                         location-scale student t-distributed ('lst').
//'                         Defaults to 5 if not otherwise is specified.
//'   \item \code{eta: } The heteroscedasticity factor
//'   \item \code{eta0: } The proportion of "base standard deviation"
//'   \item \code{qpos: } Position of systematic differences in
//'                       non-selectivity. 0 signify lower concentration
//'                       interval and 1 signify upper concentration interval.
//'   \item \code{qran: } Interquantile range where systematic differences in
//'                       non-selectivity should have its effect
//'   \item \code{prop: } average proportion of clinical samples affected by
//'                       random differences in non-selectivity
//'   \item \code{mmax: } The maximum relocation magnitude multiplier.
//'   \item \code{b0: } For systematic linear bias between IVD-MDs.
//'                     The intercept. This defaults to 0.
//'   \item \code{b1: } For systematic linear DINS between IVD-MDs.
//'                     The slope. This defaults to 1.
//'   \item \code{c0: } For systematic linear non-selectivity in IVD-MD one.
//'                     The intercept. This defaults to 0.
//'   \item \code{c1: } For systematic linear non-selectivity in IVD-MD one.
//'                     The slope. This defaults to 1.
//'   \item \code{error_dist: } The distribution to simulate measurement error
//'                             components from. Possible choices include
//'                             'norm' (normal distribution, default) and 'lt'
//'                             (location student t-distribution)
//'   \item \code{dfx: } The degrees of freedom for the measurement error
//'                      components in IVD-MD one if \code{error_dist = 'lt'}.
//'                      Defaults to 5 if not otherwise is specified.
//'   \item \code{dfy: } The degrees of freedom for the measurement error
//'                      components in IVD-MD two if \code{error_dist = 'lt'}.
//'                      Defaults to 5 if not otherwise is specified.
//'   \item \code{md_method: } Method for simulation missing data. Possible
//'                            choices include 'none'
//'                            (no missing data is simulated, default),
//'                            'mar' (Missing at Random),
//'                            'mnar' (Missing Not at Random) and 
//'                            'marmnar' (Missing systematically and randomly)
//'   \item \code{mar_prob: } The probability (value between 0 and 1) of having
//'                           a missing measurement. Only relevant if
//'                           \code{md_method = 'mar'} or 
//'                           \code{md_method = 'marmnar'}. If not specified,
//'                           but \code{md_method  = 'mar'} or
//'                           \code{md_method = 'marmnar'}, it defaults to 0.05.
//'   \item \code{mnar_threshold: } The lower bound threshold (a real value)
//'                                 for when a measurement should be missing.
//'                                 Only relevant if \code{md_method = 'mnar'}
//'                                 or \code{md_method = 'marmnar'}. If not
//'                                 specified, but \code{md_method = 'mnar'} or
//'                                 \code{md_method = 'marmnar'}, it defaults
//'                                 to \code{cil}. Alternatively,
//'                                 if not specified,
//'                                 but \code{md_method = 'mnar0'} or
//'                                 \code{md_method = 'marmnar0'}, it defaults
//'                                 to 0.
//'   
//' }
//' @param silence \code{Integer} should temporary calculation results be
//'        printed? This may be useful for debugging or strange curiosity.
//'        \code{silence = 1} signify that printing will be suppressed,
//'        which is the default. \code{silence = 0} allows such printing. 
//' 
//' @description
//' (Deprecated) Simulates a dataset with at most n x R rows measurements.
//' 
//' @details
//' (Deprecated) Simulates method comparison data for \eqn{\leq} \code{n}
//' samples (e.g., clinical samples, pools, external quality assessment
//' samples, reference material samples), where each sample is measured
//' \eqn{\leq} \code{R} times (replicated measurements). In theory, we simulate
//' \eqn{(x_{ir}, y_{ir})} from
//' 
//' \eqn{x_{ir} = f(\xi_i) + h_{ir}, \newline
//'      y_{ir} = g(f(\xi_i)) + v_{ir}, \newline
//'      \tau_i = g(f(\xi_i)) = f(\xi_i) + e_{i}}.
//' 
//' If \eqn{f} is linear in \eqn{\xi_i}, we have that
//' \eqn{\tau_i = f(\xi_i) = c_0 + c_1 \cdot \xi_i + e_{i}}.
//' We generally do not care much about the shape of \eqn{f(\xi_i)}.
//' We are more interested in the values \eqn{\tau_i} and \eqn{g}. Assuming
//' that \eqn{g} is linear, we have that
//' \eqn{g(\tau_i) = \beta_0 + \beta_1 \cdot \tau_i + e_{i}}.
//' Thus, we simulate \eqn{(x_{ir}, y_{ir})} from
//' 
//' \eqn{x_{ir} = \tau_i + h_{ir}, \newline
//'      y_{ir} = g(\tau_i) + v_{ir}, \newline
//'      g(\tau_i) = \beta_0 + \beta_1 \cdot \tau_i + e_{i}}.
//' 
//' The form of \eqn{f} is specified through parameters \code{c0} and
//' \code{c1}, whereas \eqn{g} is specified through numerous parameters such as
//' \code{b0}, \code{b1}, \code{qpos}, \code{qran}, \code{mmax}, \code{prop}.
//' The latter four parameters influence the random variable \eqn{e_{i}}.
//' 
//' \eqn{\tau_i} is modelled through \code{cil}, \code{ciu}, \code{dist} and
//' \code{df_tau}. \eqn{h_{ir}} and \eqn{v_{ir}} are measurement error
//' components modelled through \code{cvx} and \code{cvy}. These coefficients
//' of variation, can be functions of \code{error_dist}, \code{dfx},
//' \code{dfy}, \code{eta} and \code{eta0}. 
//' In order to convert the outputted \code{list} to a data frame format,
//' use either \code{as.data.frame()}, \code{as.data.table()},
//' \code{as.tibble()}. The most efficient way 
//' 
//' This function is deprecated. Use \code{simulate_eqa_data_new()} instead.
//' It is kept to allow for reproducibility for some research.
//'
//' @return
//' A \code{list} of length 4:
//' \itemize{
//'   \item \code{SampleID:} An \code{integer} vector. Sample identifiers.
//'   \item \code{ReplicateID:} An \code{integer} vector. Replicate identifiers.
//'   \item \code{MP_A:} An \code{numeric} vector. Simulated measurements from
//'                      one of the IVD-MDs in comparison.
//'   \item \code{MP_B:} An \code{numeric} vector. Simulated measurements from
//'                      the other IVD-MD in comparison.
//' } 
//'
//' @examples \dontrun{
//'   
//'   # Load data.table package from library
//'   library(data.table)
//'   
//'   # Simulate 25 clinical samples measured in triplicate affected by
//'   # random differences in non-selectivity
//'   parameters_css <- list(n = 25, R = 3, prop = 0.1, mmax = 5, cil = 25, ciu = 75)
//'   simulated_css <- simulate_eqa_data(parameters = parameters_css)
//'   
//'   # Simulate 3 external quality assessment material samples
//'   # measured in duplicate not affected by differences in non-selectivity
//'   parameters_eqams <- list(n = 3, R = 2, b0 = 0.1, b1 = 1.1)
//'   simulated_eqams <- simulate_eqa_data(parameters = parameters_eqams)
//'   
//'   # We can assume that tau_i ~ lst(df_tau = 10, mu_tau = 50, var_tau = 78.583)
//'   parameters_css <- c(parameters_css, dist = "lst", df_tau = 10)
//'   simulated_css_lst <- simulate_eqa_data(parameters = parameters_css)
//'   
//'   # We can convert the list objects to data.table objects using setDT()
//'   setDT(simulated_css)
//'   setDT(simulated_eqams)
//'   setDT(simulated_css_lst)
//'   
//'   # Print results
//'   print(simulated_css)
//'   print(simulated_eqams)
//'   print(simulated_css_lst)
//'   
//' }
//'

// [[Rcpp::export]]
List simulate_eqa_data(List parameters, int silence = 1){
  
  // Assume that no parameters are given as default
  int n_exists = 0;
  int R_exists = 0;
  int cvx_exists = 0;
  int cvy_exists = 0;
  int cil_exists = 0;
  int ciu_exists = 0;
  int dist_exists = 0;
  int df_tau_exists = 0;
  int eta_exists = 0;
  int eta0_exists = 0;
  int qpos_exists = 0;
  int qran_exists = 0;
  int prop_exists = 0;
  int mmax_exists = 0;
  int b0_exists = 0;
  int b1_exists = 0;
  int c0_exists = 0;
  int c1_exists = 0;
  int error_dist_exists = 0;
  int dfx_exists = 0;
  int dfy_exists = 0;
  int md_method_exists = 0;
  int mar_prob_exists = 0;
  int mnar_threshold_exists = 0;
  
  // Checks which of the parameters found in the 'parameters' argument
  CharacterVector given_parameters = parameters.names();
  int number_given_parameters = given_parameters.size();
  
  for(int i = 0; i < number_given_parameters; ++i){
    
    // Print to console which parameters are given
    if(silence == 0){
      Rcout << i + 1 <<". given parameter of " << number_given_parameters << " in total: " << given_parameters[i] << "\n";  
    }
    
    CharacterVector candidate_param(24);
    candidate_param[0] = "n";
    candidate_param[1] = "R";
    candidate_param[2] = "cvx";
    candidate_param[3] = "cvy";
    candidate_param[4] = "cil";
    candidate_param[5] = "ciu";
    candidate_param[6] = "dist";
    candidate_param[7] = "df_tau";
    candidate_param[8] = "eta";
    candidate_param[9] = "eta0";
    candidate_param[10] = "qpos";
    candidate_param[11] = "qran";
    candidate_param[12] = "prop";
    candidate_param[13] = "mmax";
    candidate_param[14] = "b0";
    candidate_param[15] = "b1";
    candidate_param[16] = "c0";
    candidate_param[17] = "c1";
    candidate_param[18] = "error_dist";
    candidate_param[19] = "dfx";
    candidate_param[20] = "dfy";
    candidate_param[21] = "md_method";
    candidate_param[22] = "mar_prob";
    candidate_param[23] = "mnar_threshold";
    
    
    if(candidate_param[0] == given_parameters[i]){
      n_exists = 1;
    }
    else if(candidate_param[1] == given_parameters[i]){
      R_exists = 1;
    }
    else if(candidate_param[2] == given_parameters[i]){
      cvx_exists = 1;
    }
    else if(candidate_param[3] == given_parameters[i]){
      cvy_exists = 1;
    }
    else if(candidate_param[4] == given_parameters[i]){
      cil_exists = 1;
    }
    else if(candidate_param[5] == given_parameters[i]){
      ciu_exists = 1;
    }
    else if(candidate_param[6] == given_parameters[i]){
      dist_exists = 1;
    }
    else if(candidate_param[7] == given_parameters[i]){
      df_tau_exists = 1;
    }
    else if(candidate_param[8] == given_parameters[i]){
      eta_exists = 1;
    }
    else if(candidate_param[9] == given_parameters[i]){
      eta0_exists = 1;
    }
    else if(candidate_param[10] == given_parameters[i]){
      qpos_exists = 1;
    }
    else if(candidate_param[11] == given_parameters[i]){
      qran_exists = 1;
    }
    else if(candidate_param[12] == given_parameters[i]){
      prop_exists = 1;
    }
    else if(candidate_param[13] == given_parameters[i]){
      mmax_exists = 1;
    }
    else if(candidate_param[14] == given_parameters[i]){
      b0_exists = 1;
    }
    else if(candidate_param[15] == given_parameters[i]){
      b1_exists = 1;
    }
    else if(candidate_param[16] == given_parameters[i]){
      c0_exists = 1;
    }
    else if(candidate_param[17] == given_parameters[i]){
      c1_exists = 1;
    }
    else if(candidate_param[18] == given_parameters[i]){
      error_dist_exists = 1;
    }
    else if(candidate_param[19] == given_parameters[i]){
      dfx_exists = 1;
    }
    else if(candidate_param[20] == given_parameters[i]){
      dfy_exists = 1;
    }
    else if(candidate_param[21] == given_parameters[i]){
      md_method_exists = 1;
    }
    else if(candidate_param[22] == given_parameters[i]){
      mar_prob_exists = 1;
    }
    else if(candidate_param[23] == given_parameters[i]){
      mnar_threshold_exists = 1;
    }
  }
  
  if(silence == 0){
    Rcout << "Does n exists? 1 for yes and 0 for no : " << n_exists << "\n";
    Rcout << "Does R exists? 1 for yes and 0 for no : " << R_exists << "\n";
    Rcout << "Does cvx exists? 1 for yes and 0 for no : " << cvx_exists << "\n";
    Rcout << "Does cvy exists? 1 for yes and 0 for no : " << cvy_exists << "\n";
    Rcout << "Does cil exists? 1 for yes and 0 for no : " << cil_exists << "\n";
    Rcout << "Does ciu exists? 1 for yes and 0 for no : " << ciu_exists << "\n";
    Rcout << "Does dist exists? 1 for yes and 0 for no : " << dist_exists << "\n";
    Rcout << "Does df_tau exists? 1 for yes and 0 for no : " << df_tau_exists << "\n";
    Rcout << "Does eta exists? 1 for yes and 0 for no : " << eta_exists << "\n";
    Rcout << "Does eta0 exists? 1 for yes and 0 for no : " << eta0_exists << "\n";
    Rcout << "Does qpos exists? 1 for yes and 0 for no : " << qpos_exists << "\n";
    Rcout << "Does qran exists? 1 for yes and 0 for no : " << qran_exists << "\n";
    Rcout << "Does prop exists? 1 for yes and 0 for no : " << prop_exists << "\n";
    Rcout << "Does mmax exists? 1 for yes and 0 for no : " << mmax_exists << "\n";
    Rcout << "Does b0 exists? 1 for yes and 0 for no : " << b0_exists << "\n";
    Rcout << "Does b1 exists? 1 for yes and 0 for no : " << b1_exists << "\n";
    Rcout << "Does c0 exists? 1 for yes and 0 for no : " << c0_exists << "\n";
    Rcout << "Does c1 exists? 1 for yes and 0 for no : " << c1_exists << "\n";
    Rcout << "Does error_dist exists? 1 for yes and 0 for no : " << error_dist_exists << "\n";
    Rcout << "Does dfx exists? 1 for yes and 0 for no : " << dfx_exists << "\n";
    Rcout << "Does dfy exists? 1 for yes and 0 for no : " << dfy_exists << "\n";
    Rcout << "Does md_method exists? 1 for yes and 0 for no : " << md_method_exists << "\n";
    Rcout << "Does mar_prob exists? 1 for yes and 0 for no : " << mar_prob_exists << "\n";
    Rcout << "Does mnar_threshold exists? 1 for yes and 0 for no : " << mnar_threshold_exists << "\n";
  }
  
  // Base parameters
  int n = 0;
  int R = 0;
  float cvx = 0;
  float cvy = 0;
  float cil = 0;
  float ciu = 0;
  
  // Heteroscedasticity parameters
  float eta = 0;
  float eta0 = 0;
  
  // Random and systematic differences in non-selectivity parameters
  int qpos = 0;
  float qran = 0;
  float prop = 0;
  float mmax = 0;
  
  // Consistent linear differences in non-selectivity parameters
  float b0 = 0;
  float b1 = 1;
  float c0 = 0;
  float c1 = 1;
  
  // Modelling assumptions
  std::string dist = "unif";
  std::string error_dist = "norm";
  
  // Parameters for modelling assumptions
  float sigma_tau = 0;
  float mu_tau = 0;
  float df_tau = 5.0;
  float dfx = 5.0;
  float dfy = 5.0;
  
  // Parameters for modelling missing data
  std::string md_method = "none";
  float mar_prob = 0;
  float mnar_threshold = 0;
  
  
  if(n_exists == 1){
    int reg_n = parameters["n"];
    n = n + reg_n;  
  }
  else{
    int gen_n = R::rpois(25);
    IntegerVector ns (2);
    ns[0] = gen_n;
    ns[1] = 30;
    gen_n = min(ns);
    ns[0] = gen_n;
    ns[1] = 20;
    n = max(ns);
  }
  if(R_exists == 1){
    int reg_R = parameters["R"];
    R = reg_R;     
  }
  else{
    float u = R::runif(0,1);
    if(u < 0.10){
      R += 2;
    }
    else if(u < 0.95){
      R += 3; 
    }
    else{
      R += 4;
    }
  }
  if(cvx_exists == 1){
    float reg_cvx = parameters["cvx"];
    cvx += reg_cvx;
  }
  else{
    cvx += R::rbeta(2, 5) / 10.0;
  }
  if(cvy_exists == 1){
    float reg_cvy = parameters["cvy"];
    cvy += reg_cvy;
  }
  else{
    cvy += R::rbeta(2, 5) / 10.0;
  }
  if(cil_exists == 1){
    float reg_cil = parameters["cil"];
    cil += reg_cil;
  }
  else{
    cil += R::rf(1.057057, 8.15) * 44;
    if(cil < 0.01){
      cil += R::rf(1.057057, 8.15) * 44;
      if(cil < 0.01){
        cil += R::rf(1.057057, 8.15) * 44;
        if(cil < 0.01){
          cil += R::rf(1.057057, 8.15) * 44;
        }
      }
    }
  }
  if(ciu_exists == 1){
    float reg_ciu = parameters["ciu"];
    ciu += reg_ciu;
  }
  else{
    if(cil > 0){
      float multiplier = R::rbeta(0.78, 11) * 44;
      ciu += cil + cil * multiplier;
    }
    else{
      stop("cil is negative or zero. This should not be possible! Check if something is wrong");
    }
  }
  
  if(dist_exists == 1){
    std::string reg_dist = parameters["dist"];
    dist = reg_dist;
    if(dist == "norm"){
      mu_tau += 0.5 * (cil + ciu);
      sigma_tau += 0.2149292 * (ciu - cil);
    }
    else if(dist == "lst"){
      if(df_tau_exists == 1){
        float reg_df_tau = parameters["df_tau"];
        df_tau = reg_df_tau;
      }
      mu_tau += 0.5 * (cil + ciu);
      sigma_tau += 0.49 * (ciu - cil) * (1.0 / R::qt(0.99, df_tau, 0, 0));
    }
    else if(dist == "lnorm"){
      mu_tau += 0.5 * (log(cil + 0.99 * (ciu - cil)) + log(cil + 0.01 * (ciu - cil)));
      sigma_tau += 0.2149292 * log((cil + 0.99 * (ciu - cil))/(cil + 0.01 * (ciu - cil)));
    }
    else{
      mu_tau += 0.5 * (cil + ciu);
      sigma_tau += (1.0 / sqrt(12.0)) * (ciu - cil);
    }
  }
  
  if(error_dist_exists == 1){
    std::string reg_error_dist = parameters["error_dist"];
    error_dist = reg_error_dist;
    if(error_dist != "norm" and (error_dist == "lt" or error_dist == "t" or error_dist == "lst")){
      error_dist = "lt";
      if(dfx_exists == 1){
        float reg_dfx = parameters["dfx"];
        dfx = reg_dfx;
      }
      if(dfy_exists == 1){
        float reg_dfy = parameters["dfy"];
        dfy = reg_dfy;
      }
    }
  }
  
  if(eta_exists == 1){
    float reg_eta = parameters["eta"];
    eta += reg_eta;
  }
  else{
    eta = 1;
  }
  if(eta0_exists == 1){
    float reg_eta0 = parameters["eta0"];
    eta0 = reg_eta0;
  }
  else{
    eta0 = 1;
  }
  
  if(b0_exists == 1){
    float reg_b0 = parameters["b0"];
    b0 = reg_b0;
  }
  
  if(b1_exists == 1){
    float reg_b1 = parameters["b1"];
    b1 = reg_b1;
  }
  
  if(c0_exists == 1){
    float reg_c0 = parameters["c0"];
    c0 = reg_c0;
  }
  
  if(c1_exists == 1){
    float reg_c1 = parameters["c1"];
    c1 = reg_c1;
  }
  
  if(qran_exists == 1 and prop_exists == 1){
    Rcout << "Both qran and prop are found in parameters" << "\n";
    Rcout << "Make sure only one of them is used" << "\n";
    Rcout << "prop is removed in favor of qran" << "\n";
    prop_exists -= 1;
  }
  if(qpos_exists == 1){
    float reg_qpos = parameters["qpos"];
    qpos += reg_qpos;
  }
  else{
    qpos -= 1;
  }
  if(qran_exists == 1){
    float reg_qran = parameters["qran"];
    qran += reg_qran;
  }
  else{
    qran += 0;
  }
  if(prop_exists == 1){
    float reg_prop = parameters["prop"];
    prop += reg_prop;
  }
  else{
    prop = 0;
  }
  if(mmax_exists == 1){
    float reg_mmax = parameters["mmax"];
    mmax += reg_mmax;
  }
  
  if(md_method_exists == 1){
    std::string reg_md_method = parameters["md_method"];
    md_method = reg_md_method;
    
    if(md_method == "mar"){
      if(mar_prob_exists == 1){
        float reg_mar_prob = parameters["mar_prob"];
        mar_prob += reg_mar_prob;
      }
      else{
        mar_prob += 0.05;
      }
    }
    else if(md_method == "mnar"){
      if(mnar_threshold_exists == 1){
        float reg_mnar_threshold = parameters["mnar_threshold"];
        mnar_threshold += reg_mnar_threshold;
      }
      else{
        mnar_threshold += cil;
      }
    }
    else if(md_method == "mnar0"){
      if(mnar_threshold_exists == 1){
        float reg_mnar_threshold = parameters["mnar_threshold"];
        mnar_threshold += reg_mnar_threshold;
      }
    }
    else if(md_method == "marmnar"){
      if(mnar_threshold_exists == 1){
        float reg_mnar_threshold = parameters["mnar_threshold"];
        mnar_threshold += reg_mnar_threshold;
      }
      else{
        mnar_threshold += cil;
      }
      if(mar_prob_exists == 1){
        float reg_mar_prob = parameters["mar_prob"];
        mar_prob += reg_mar_prob;
      }
      else{
        mar_prob += 0.05;
      }
    }
    
    else if(md_method == "marmnar0"){
      if(mnar_threshold_exists == 1){
        float reg_mnar_threshold = parameters["mnar_threshold"];
        mnar_threshold += reg_mnar_threshold;
      }
      if(mar_prob_exists == 1){
        float reg_mar_prob = parameters["mar_prob"];
        mar_prob += reg_mar_prob;
      }
      else{
        mar_prob += 0.05;
      }
    }
    else{
      mnar_threshold -= 999999.99;
    }
  }
  else{
    mnar_threshold -= 999999.99;
  }

  NumericVector unsorted_tau(n);
  NumericVector tau(n);
  int nR = n * R;
  IntegerVector SampleID(nR);
  IntegerVector ReplicateID(nR);
  NumericVector MP_A(nR);
  NumericVector MP_B(nR);
  
  float average = 0.5 * (cil + ciu);
  if(dist == "lnorm"){
    average = exp(mu_tau + pow(sigma_tau, 2.0) / 2.0);
  }
  
  int above = R::rbinom(1, 0.5);
  if(above == 1){
    above = 1;
  }
  else{
    above = -1;
  }
  
  if(silence == 0){
    Rcout << "The value of n is : " << n << "\n";
    Rcout << "The value of R is : " << R << "\n";
    Rcout << "The value of cvx is : " << cvx << "\n";
    Rcout << "The value of cvy is : " << cvy << "\n";
    Rcout << "The value of cil is : " << cil << "\n";
    Rcout << "The value of ciu is : " << ciu << "\n";
    Rcout << "The value of eta is : " << eta << "\n";
    Rcout << "The value of eta0 is : " << eta0 << "\n";
    Rcout << "The value of qpos is : " << qpos << "\n";
    Rcout << "The value of qran is : " << qran << "\n";
    Rcout << "The value of prop is : " << prop << "\n";
    Rcout << "The value of mmax is : " << mmax << "\n";    
  }
  
  float base_x = average * cvx;
  float base_y = average * cvy;
  
  float beg_sdx = base_x * eta0;
  float end_sdx = base_x * eta * eta0;
  float seg_sdx = 0;
  
  if(eta0 > 0 and eta > 0){
    seg_sdx = (end_sdx - beg_sdx) / n;
  }

  float beg_sdy = base_y * eta0;
  float end_sdy = base_y * eta * eta0;
  float seg_sdy = 0;
  
  if(eta0 > 0 and eta > 0){
    seg_sdy = (end_sdy - beg_sdy) / n;
  }
  
  if(dist == "norm"){
    for(int i = 0; i < n; ++i){
      unsorted_tau[i] = R::rnorm(mu_tau, sigma_tau);  
    }  
  }
  else if(dist == "lst"){
    for(int i = 0; i < n; ++i){
      unsorted_tau[i] = mu_tau + sigma_tau * R::rt(df_tau);
    }  
  }
  else if(dist == "lnorm"){
    for(int i = 0; i < n; ++i){
      unsorted_tau[i] = R::rlnorm(mu_tau, sigma_tau);  
    }  
  }
  else{
    for(int i = 0; i < n; ++i){
      unsorted_tau[i] = R::runif(cil, ciu);  
    }
  }
    
  
  if(eta0_exists == 1 and eta_exists == 1){
    tau = unsorted_tau.sort();  
  }
  
  for(int i = 0; i < n; ++i){
    
    float hetero_extra_x = i * seg_sdx;
    float hetero_extra_y = i * seg_sdy;
    float sdx = beg_sdx + hetero_extra_x;
    float sdy = beg_sdy + hetero_extra_y;
    
    if(silence == 0){
      Rcout << "sdx for i = " << i << " is " << sdx << "\n";
      Rcout << "sdy for i = " << i << " is " << sdy << "\n";  
    }
    
    
    if(sdx <= 0 or sdy <= 0){
      stop("Your choices of eta and eta0 resulted in negative standard deviations. Calculations are terminated");
    }
    
    if(eta0_exists == 0 or eta_exists == 0){
      tau[i] = unsorted_tau[i];  
    }

    float limit = 0;
    float relocating_magnitude = 0;
    int relocate_sample_i = 0;
    
    // if prop exists, and CS is dins-affected, randomly select relocation magnitude from from Beta(2, 2) * mmax 
    if(prop_exists == 1){
      relocate_sample_i += R::rbinom(1, prop);
      float sign_relocated_sample_i = R::runif(0, 1);
      if(sign_relocated_sample_i < 0.5){
        sign_relocated_sample_i = -1;
      }
      else{
        sign_relocated_sample_i = 1;
      }
      relocating_magnitude += R::rbeta(2, 2) * sign_relocated_sample_i * mmax;
    }
    else if(qran_exists == 1){
      if(qpos == 0){
        limit += cil + qran * (ciu - cil); 
        if(tau[i] <= limit){
          relocate_sample_i = relocate_sample_i + 1;
          float rel_diff = (limit - tau[i]) / (limit - cil);
          relocating_magnitude = rel_diff * mmax * above;
        }
      }
      else if(qpos == 1){
        limit += ciu - qran * (ciu - cil);
        if(tau[i] >= limit){
          relocate_sample_i = relocate_sample_i + 1;
          float rel_diff = (tau[i] - limit) / (ciu - limit);
          relocating_magnitude = rel_diff * mmax * above;
        }
      }
    }
    int lower = i * R;
    int upper = R * (1 + i) - 1;
    IntegerVector idx = Rcpp::Range(lower, upper);
    if(error_dist == "norm"){
      for(int r = 0; r < R; ++r){
        SampleID[idx[r]] = i + 1;
        ReplicateID[idx[r]] = r + 1;
        MP_B[idx[r]] = c0 + tau[i] * c1;
        MP_A[idx[r]] = b0 + MP_B[idx[r]] * b1 + R::rnorm(0, sdy) + relocate_sample_i * relocating_magnitude * sqrt(pow(sdx, 2) + pow(sdy, 2));
        MP_B[idx[r]] = MP_B[idx[r]] + R::rnorm(0, sdx);
        if(md_method == "mar" or md_method == "mnar" or md_method == "mnar0" or md_method == "marmnar" or md_method == "marmnar0"){
          bool mar_r_A = R::runif(0, 1) < mar_prob;
          bool mar_r_B = R::runif(0, 1) < mar_prob;
          bool mnar_r_A = MP_A[idx[r]] < mnar_threshold;
          bool mnar_r_B = MP_B[idx[r]] < mnar_threshold;
          if(mar_r_A or mnar_r_A){
            MP_A[idx[r]] = NA_REAL;
          }
          if(mar_r_B or mnar_r_B){
            MP_B[idx[r]] = NA_REAL;
          }
        }
        else{
          if(MP_A[idx[r]] < 0){
            MP_A[idx[r]] = MP_A[idx[r]] * (-1);
          }
          if(MP_B[idx[r]] < 0){
            MP_B[idx[r]] = MP_B[idx[r]] * (-1);
          }  
        }
      }  
    }
    else if(error_dist == "lt"){
      for(int r = 0; r < R; ++r){
        SampleID[idx[r]] = i + 1;
        ReplicateID[idx[r]] = r + 1;
        MP_B[idx[r]] = c0 + tau[i] * c1;
        MP_A[idx[r]] = b0 + MP_B[idx[r]] * b1 + sdy * R::rt(dfy) + relocate_sample_i * relocating_magnitude * sqrt(pow(sdx, 2) + pow(sdy, 2));
        MP_B[idx[r]] = MP_B[idx[r]] + sdx * R::rt(dfx);
        if(md_method == "mar" or md_method == "mnar" or md_method == "mnar0" or md_method == "marmnar" or md_method == "marmnar0"){
          bool mar_r_A = R::runif(0, 1) < mar_prob;
          bool mar_r_B = R::runif(0, 1) < mar_prob;
          bool mnar_r_A = MP_A[idx[r]] < mnar_threshold;
          bool mnar_r_B = MP_B[idx[r]] < mnar_threshold;
          if(mar_r_A or mnar_r_A){
            MP_A[idx[r]] = NA_REAL;
          }
          if(mar_r_B or mnar_r_B){
            MP_B[idx[r]] = NA_REAL;
          }
        }
        else{
          if(MP_A[idx[r]] < 0){
            MP_A[idx[r]] = MP_A[idx[r]] * (-1);
          }
          if(MP_B[idx[r]] < 0){
            MP_B[idx[r]] = MP_B[idx[r]] * (-1);
          }  
        }
      }  
    }
  }
  
  List out = List::create(Named("SampleID") = SampleID, Named("ReplicateID") = ReplicateID, Named("MP_A") = round(MP_A, 6), Named("MP_B") = round(MP_B, 6));
  return out;
}

// Allow for custom functions
NumericVector custom_function(NumericVector x, Rcpp::Function func){
  NumericVector output = func(x);
  return output;
} 

//' @title Simulation of External Quality Assessment (EQA) Data
//' @name sim_eqa_data
//'
//' @param parameters A \code{list} of parameters and their value. These values
//'        are used to simulate the EQA data. Optional parameters that you may
//'        include into the list are:
//' \itemize{
//'   \item \code{n: } The number of samples.
//'   \item \code{R: } The number of replicates on each sample.
//'   \item \code{cvx: } The repeatability coefficient of variation for IVD-MD 
//'                      \code{MP_B}.
//'   \item \code{cvy: } The repeatability coefficient of variation for IVD-MD
//'                      \code{MP_A}.
//'   \item \code{cil: } The lower bound of the concentration interval.
//'   \item \code{ciu: } The upper bound of the concentration interval.
//'   \item \code{dist: } The distribution to simulate latent variables from.
//'                       Possbile choices include
//'                       \code{unif} (uniform distribution, default),
//'                       \code{norm} (normal distribution),
//'                       \code{lst} (location-scale student t-distribution),
//'                       \code{lnorm} (log-normal distribution).
//'   \item \code{df_tau: } The degrees of freedom for the \code{lst}
//'                         distribution if the distribution of latent variables
//'                         are location-scale student t-distributed ('lst').
//'                         Defaults to \eqn{5} if not otherwise is specified.
//'   \item \code{eta: } The heteroscedasticity factor \eqn{\eta}.
//'   \item \code{eta0: } The proportion of "base standard deviation". See
//'                       details.
//'   \item \code{qpos: } Position of systematic differences in non-selectivity.
//'                       \code{0} signify lower concentration interval and
//'                       \code{1} signify upper.
//'   \item \code{qran: } Interquantile range where systematic differences in
//'                       non-selectivity should have its effect
//'   \item \code{prop: } average proportion of clinical samples affected by
//'                       random differences in non-selectivity.
//'   \item \code{mmax: } The maximum relocation magnitude multiplier,
//'                       \eqn{m_{\max}}. This number is a multiplier of the
//'                       base standard deviation. This will only have an effect
//'                       if \code{prop}, or \code{qpos} and \code{qran} are
//'                       specified as well.
//'   \item \code{b0: } For systematic linear bias between IVD-MDs.
//'                     The intercept. Defaults to 0.
//'   \item \code{b1: } For systematic linear bias between IVD-MDs. The 
//'                     slope. Defaults to 1.
//'   \item \code{c0: } For systematic linear non-selectivity in IVD-MD 1.
//'                     The intercept. Defaults to 0.
//'   \item \code{c1: } For systematic linear non-selectivity in IVD-MD 1.
//'                     The slope. Defaults to 1.
//'   \item \code{error_dist: } The distribution to simulate measurement error
//'                             components from. Possible choices include 'norm'
//'                             (normal distribution, default) and
//'                             'lt' (location student t-distribution).
//'   \item \code{dfx: } The degrees of freedom for the measurement error
//'                      components in IVD-MD 1 if \code{error_dist = 'lt'}.
//'                      Defaults to \eqn{5} if not otherwise is specified.
//'   \item \code{dfy: } The degrees of freedom for the measurement error
//'                      components in IVD-MD 2 if \code{error_dist = 'lt'}.
//'                      Defaults to \eqn{5} if not otherwise is specified.
//'   \item \code{md_method: } Method for simulating missing data. Possible
//'                            choices include \code{'none'}
//'                            (no missing data is simulated, default),
//'                            \code{'mar'} (missing at random),
//'                            \code{'mnar'} (missing not at random) and
//'                            \code{'marmnar'}
//'                            (missing both at random and not at random)
//'   \item \code{mar_prob: } The probability (value between 0 and 1) of having
//'                           a missing measurement. Only relevant if
//'                           \code{md_method} is \code{'mar'} or
//'                           \code{'marmnar'}. If not specified, but 
//'                           \code{md_method = 'mar'} or
//'                           \code{md_method = 'marmnar'}, it defaults to 0.05.
//'   \item \code{mnar_threshold: } The lower bound threshold (\code{double})
//'                                 for when a measurement should be missing.
//'                                 Only relevant if \code{md_method} is
//'                                 \code{'mnar'} or \code{'marmnar'}. If not
//'                                 specified, but \code{md_method = 'mnar'}
//'                                 or \code{md_method = 'marmnar'}, it defaults
//'                                 to \code{cil}. Alternatively, if not
//'                                 specified, but \code{md_method = 'mnar0'}
//'                                 or \code{md_method = 'marmnar0'},
//'                                 it defaults to 0.
//'   \item \code{g: } A custom function to simulate \eqn{\mathrm{E}[y_{ir}]}
//'                    from. It is not recommended to use this for linear
//'                    \code{g}, because it is much slower. This is generally
//'                    used to simulate data from non-linear functions other
//'                    than the three built-in alternatives
//'                    (see \code{type} argument).
//'   \item \code{obs_tau: } A \code{numeric} vector of \eqn{\tau_i} values.
//'                          (true but latent values for IVD-MD one).
//'                          If this parameter is specified, \code{n} will take
//'                          the length of \code{obs_tau}. Moreover, all
//'                          parameters used to generate \eqn{\tau_i} are
//'                          ignored.
//' }
//' @param type A \code{integer}. Custom built-in functions \code{g}. Set to
//'        \code{0}, to use a linear \code{g}. Otherwise, particular custom
//'        non-linear functions can be used. The four different built-in
//'        \code{g}:
//'        \itemize{
//'           \item \code{type = 0: } \eqn{g(\tau_i) = \beta_0 + \beta_1 \cdot
//'                                        \tau_i}
//'           \item \code{type = 1: } \eqn{g(\tau_i) = \tau_i + 0.90 \cdot
//'                                        \mathrm{sin}(0.40 \cdot
//'                                        \tau_i ^ {1.06})}
//'           \item \code{type = 2: } \eqn{g(\tau_i) = \tau_i + 0.05 \cdot
//'                                        \exp(0.16 \cdot \tau_i ^ {1.35})}
//'           \item \code{type = 3: } \eqn{g(\tau_i) = \tau_i - \exp[
//'                                        -0.125 \cdot (\tau_i - 1.50)^2]}
//'        }
//'        If \code{g} exists in \code{parameters}, this argument is ignored.
//' @param AR A \code{logical} value. If \code{TRUE}, data is simulated including
//'        replicated measurements. Otherwise, mean of replicated measurements
//'        are returned (MOR).
//' @param include_parameters A \code{logical} value. If \code{TRUE}, the used
//'        parameters in the data simulation is saved and placed in a seperate
//'        list as part of the output.
//'
//' @description
//' Simulates a dataset with at most n x R rows measurements.
//' 
//' @details
//' Simulates method comparison data for \code{n} samples
//' (e.g., clinical samples, pools, external quality assessment samples, reference material samples),
//' where each sample is measured at most \code{R} times (replicated measurements).
//' In theory, we simulate measurement pairs \eqn{(x_{ir}, y_{ir})} from:
//' 
//' \eqn{x_{ir} = f(\xi_i) + h_{ir} \newline
//'      y_{ir} = g(f(\xi_i)) + v_{ir} \newline
//'      g(f(\xi_i)) = f(\xi_i) + e_{i}}.
//' 
//' If \eqn{f} is linear in \eqn{\xi_i}, we have that
//' \eqn{\tau_i = f(\xi_i) = c_0 + c_1 \, \xi_i}.
//' We generally do not care much about the shape of \eqn{f(\xi_i)}. We are
//' more interested in the values \eqn{\tau_i} and the shape of \eqn{g}.
//' Assuming that \eqn{g} is linear, we have that
//' 
//' \eqn{g(\tau_i) = \beta_0 + \beta_1 \cdot \tau_i + e_{i}}.
//' 
//' Thus, we simulate \eqn{(x_{ir}, y_{ir})} from:
//' 
//' \eqn{x_{ir} = \tau_i + h_{ir} \newline
//'      y_{ir} = g(\tau_i) + v_{ir}},
//' 
//' where \eqn{g} is some function.
//' 
//' If we wish to specify \eqn{f} (not very informative) we can specify it
//' through \code{c0} and \code{c1}. \eqn{g} is specified through numerous
//' parameters such as \code{b0}, \code{b1}, \code{qpos}, \code{qran},
//' \code{mmax}, \code{prop} and \code{type}. One can also specify \eqn{g}
//' using a custom R function. See explanation of \code{parameters}.
//' 
//' \eqn{\tau_i} is a random variable modelled through \code{cil}, \code{ciu},
//' \code{dist} and \code{df_tau}.
//' 
//' \eqn{h_{ir}} and \eqn{v_{ir}} are the measurement error components and
//' are modelled through \code{cvx}, \code{cvy}. The coefficients of variation,
//' \code{cvx}, \code{cvy} can be functions of \code{error_dist}, \code{dfx},
//' \code{dfy}, \code{eta} and \code{eta0}.
//' 
//' In order to convert the outputted \code{list} to a data frame structure,
//' use either \code{as.data.frame()}, \code{as.data.table()} or
//' \code{as.tibble()}.  The most efficient way to convert the output is
//' however \code{setDT()} from the \code{data.table} package.
//' 
//' Note, there is no guarantee that every simulated \eqn{(x_{ir}, y_{ir})}
//' will take positive values. To avoid negative values, you will need to
//' select values of \code{cil}, \code{ciu}, \code{cvx}, \code{cvy}, and
//' \code{mmax} that result in a negligble probability of generating negative
//' values. In general, avoid using large values of \code{cvx} and \code{cvy}
//' if \code{cil} is close to zero and \code{ciu} - \code{ciu} is large.
//' Negative values can alternatively be replaced with \code{NA} values if
//' \code{md_method = 'mnar0'}.
//'
//' @return
//' A \code{list} of length three or four:
//' \itemize{
//'   \item \code{SampleID:} An \code{integer} vector. Sample identifiers.
//'   \item \code{ReplicateID:} An \code{integer} vector. Replicate identifiers.
//'                             Only appears if \code{AR = TRUE}.
//'   \item \code{MP_A:} An \code{numeric} vector. Simulated measurements from
//'                      one of the IVD-MDs in comparison.
//'   \item \code{MP_B:} An \code{numeric} vector. Simulated measurements from
//'                      the other IVD-MD in comparison.
//' }
//'
//' @examples
//' # Load package
//' library(fasteqa)
//' 
//' # Simulate 25 clinical samples measure in triplicate affected by random
//' # differences in non-selectivity
//' 
//' parameters_cs <- list(n = 25,
//'                       R = 3,
//'                       cvx = 0.01,
//'                       cvy = 0.01,
//'                       cil = 2,
//'                       ciu = 10,
//'                       prop = 0.2,
//'                       mmax = 5)
//' simulated_cs <- sim_eqa_data(parameters_cs,
//'                              type = 0,
//'                              AR = TRUE)
//' print(as.data.frame(simulated_cs))
//' 
//' 
//' # Simulate 30 clinical samples and 3 External Quality Assessment material
//' # samples (all measured in duplicate) that are commutable with the
//' # clinical samples.
//' 
//' parameters_cs <- list(n = 30,
//'                       R = 2,
//'                       cvx = 0.01,
//'                       cvy = 0.01,
//'                       cil = 5,
//'                       ciu = 90,
//'                       qran = 0.30,
//'                       qpos = 1,
//'                       mmax = 7.5)
//'
//' simulated_cs <- sim_eqa_data(parameters_cs,
//'                              type = 0,
//'                              AR = TRUE,
//'                              include_parameters = TRUE)
//' true_values_eqa_samples <- c(10, 40, 75)                             
//' parameters_eq <- simulated_cs$parameters
//' parameters_eq <- c(list(obs_tau = true_values_eqa_samples),
//'                    parameters_eq)
//' simulated_cs <- simulated_cs$simulated_data
//' simulated_eq <- sim_eqa_data(parameters_eq,
//'                              type = 0,
//'                              AR = TRUE)
//' print(as.data.frame(simulated_cs))
//' print(as.data.frame(simulated_eq))
//'                              
//'
 
 
 // [[Rcpp::export]]
 List sim_eqa_data(Nullable<List> parameters = R_NilValue,
                   int type = 0,
                   bool AR = false,
                   bool include_parameters = false){
   
   
   List params;
   if (parameters.isNotNull()) {
     params = as<List>(parameters);
   } else {
     // Default parameters if nothing is specified in parameters
     params = List::create(Named("n") = 25,
                           Named("R") = 3,
                           Named("cvx") = 0.01,
                           Named("cvy") = 0.01,
                           Named("cil") = 2.0,
                           Named("ciu") = 10.0);
   }
   
   // Assume that no parameters are given as default
   int n_exists = 0;
   int R_exists = 0;
   int cvx_exists = 0;
   int cvy_exists = 0;
   int cil_exists = 0;
   int ciu_exists = 0;
   int dist_exists = 0;
   int df_tau_exists = 0;
   int eta_exists = 0;
   int eta0_exists = 0;
   int qpos_exists = 0;
   int qran_exists = 0;
   int prop_exists = 0;
   int mmax_exists = 0;
   int b0_exists = 0;
   int b1_exists = 0;
   int c0_exists = 0;
   int c1_exists = 0;
   int error_dist_exists = 0;
   int dfx_exists = 0;
   int dfy_exists = 0;
   int md_method_exists = 0;
   int mar_prob_exists = 0;
   int mnar_threshold_exists = 0;
   int qdir_exists = 0;
   int obs_tau_exists = 0;
   int qlim_exists = 0;
   int qnor_exists = 0;
   int sdx_exists = 0;
   int sdy_exists = 0;
   int g_exists = 0;
   
   // Checks which of the parameters found in the 'parameters' argument
   CharacterVector given_parameters = params.names();
   int number_given_parameters = given_parameters.size();
   
   for(int i = 0; i < number_given_parameters; ++i){
     
     CharacterVector candidate_param(31);
     candidate_param[0] = "n";
     candidate_param[1] = "R";
     candidate_param[2] = "cvx";
     candidate_param[3] = "cvy";
     candidate_param[4] = "cil";
     candidate_param[5] = "ciu";
     candidate_param[6] = "dist";
     candidate_param[7] = "df_tau";
     candidate_param[8] = "eta";
     candidate_param[9] = "eta0";
     candidate_param[10] = "qpos";
     candidate_param[11] = "qran";
     candidate_param[12] = "prop";
     candidate_param[13] = "mmax";
     candidate_param[14] = "b0";
     candidate_param[15] = "b1";
     candidate_param[16] = "c0";
     candidate_param[17] = "c1";
     candidate_param[18] = "error_dist";
     candidate_param[19] = "dfx";
     candidate_param[20] = "dfy";
     candidate_param[21] = "md_method";
     candidate_param[22] = "mar_prob";
     candidate_param[23] = "mnar_threshold";
     candidate_param[24] = "qdir";
     candidate_param[25] = "obs_tau";
     candidate_param[26] = "qlim";
     candidate_param[27] = "qnor";
     candidate_param[28] = "sdx";
     candidate_param[29] = "sdy";
     candidate_param[30] = "g";
     
     
     
     if(candidate_param[0] == given_parameters[i]){
       n_exists = 1;
     }
     else if(candidate_param[1] == given_parameters[i]){
       R_exists = 1;
     }
     else if(candidate_param[2] == given_parameters[i]){
       cvx_exists = 1;
     }
     else if(candidate_param[3] == given_parameters[i]){
       cvy_exists = 1;
     }
     else if(candidate_param[4] == given_parameters[i]){
       cil_exists = 1;
     }
     else if(candidate_param[5] == given_parameters[i]){
       ciu_exists = 1;
     }
     else if(candidate_param[6] == given_parameters[i]){
       dist_exists = 1;
     }
     else if(candidate_param[7] == given_parameters[i]){
       df_tau_exists = 1;
     }
     else if(candidate_param[8] == given_parameters[i]){
       eta_exists = 1;
     }
     else if(candidate_param[9] == given_parameters[i]){
       eta0_exists = 1;
     }
     else if(candidate_param[10] == given_parameters[i]){
       qpos_exists = 1;
     }
     else if(candidate_param[11] == given_parameters[i]){
       qran_exists = 1;
     }
     else if(candidate_param[12] == given_parameters[i]){
       prop_exists = 1;
     }
     else if(candidate_param[13] == given_parameters[i]){
       mmax_exists = 1;
     }
     else if(candidate_param[14] == given_parameters[i]){
       b0_exists = 1;
     }
     else if(candidate_param[15] == given_parameters[i]){
       b1_exists = 1;
     }
     else if(candidate_param[16] == given_parameters[i]){
       c0_exists = 1;
     }
     else if(candidate_param[17] == given_parameters[i]){
       c1_exists = 1;
     }
     else if(candidate_param[18] == given_parameters[i]){
       error_dist_exists = 1;
     }
     else if(candidate_param[19] == given_parameters[i]){
       dfx_exists = 1;
     }
     else if(candidate_param[20] == given_parameters[i]){
       dfy_exists = 1;
     }
     else if(candidate_param[21] == given_parameters[i]){
       md_method_exists = 1;
     }
     else if(candidate_param[22] == given_parameters[i]){
       mar_prob_exists = 1;
     }
     else if(candidate_param[23] == given_parameters[i]){
       mnar_threshold_exists = 1;
     }
     else if(candidate_param[24] == given_parameters[i]){
       qdir_exists = 1;
     }
     else if(candidate_param[25] == given_parameters[i]){
       obs_tau_exists = 1;
     }
     else if(candidate_param[26] == given_parameters[i]){
       qlim_exists = 1;
     }
     else if(candidate_param[27] == given_parameters[i]){
       qnor_exists = 1;
     }
     else if(candidate_param[28] == given_parameters[i]){
       sdx_exists = 1;
     }
     else if(candidate_param[29] == given_parameters[i]){
       sdy_exists = 1;
     }
     else if(candidate_param[30] == given_parameters[i]){
       g_exists = 1;
     }
   }
   
   // Base parameters
   int n = 25;
   int R = 3;
   double cvx = 0;
   double cvy = 0;
   double cil = 0;
   double ciu = 0;
   
   // Heteroscedasticity parameters
   double eta = 0;
   double eta0 = 0;
   
   // Random and systematic differences in non-selectivity parameters
   int qpos = 0;
   int qdir = 0;
   double qran = 0;
   double qlim = 0;
   double qnor = 1;
   double prop = 0;
   double mmax = 0;
   
   // Systematic linear differences in non-selectivity parameters
   double b0 = 0;
   double b1 = 1;
   
   // Systematic linear non-selectivity
   double c0 = 0;
   double c1 = 1;
   
   // Modelling assumptions
   std::string dist = "unif";
   std::string error_dist = "norm";
   
   // Parameters for modelling assumptions
   double sigma_tau = 0;
   double mu_tau = 0;
   double df_tau = 5.0;
   double dfx = 5.0;
   double dfy = 5.0;
   
   // Parameters for modelling missing data
   std::string md_method = "none";
   double mar_prob = 0;
   double mnar_threshold = 0;
   
   // Generation of n
   if(n_exists == 1){
     int reg_n = params["n"];
     n = reg_n;
   }
   else{
     int gen_n = R::rpois(5) + 19;
     IntegerVector ns (2);
     ns[0] = gen_n;
     ns[1] = 30;
     gen_n = min(ns);
     ns[0] = gen_n;
     ns[1] = 20;
     n = max(ns);
   }
   // Generation of R
   if(R_exists == 1){
     int reg_R = params["R"];
     R = reg_R;
   }
   else{
     double u = R::runif(0,1);
     if(u < 0.10){
       R = 2;
     }
     else if(u < 0.95){
       R = 3;
     }
     else{
       R = 4;
     }
   }
   
   if(cil_exists == 1){
     double reg_cil = params["cil"];
     cil = reg_cil;
   }
   else{
     if(type != 0){
       cil = 2;
     }
     else{
       cil = R::rf(1.057057, 8.15) * 44;
       if(cil < 0.01){
         cil = R::rf(1.057057, 8.15) * 44;
         if(cil < 0.01){
           cil = R::rf(1.057057, 8.15) * 44;
           if(cil < 0.01){
             cil = R::rf(1.057057, 8.15) * 44;
           }
         }
       }
     }
   }
   if(ciu_exists == 1){
     double reg_ciu = params["ciu"];
     if(reg_ciu <= cil){
       Rcout << "ciu is" << reg_ciu << "which is smaller than" << cil << "\n";
       stop("ciu must be larger than cil!");
     }
     ciu = reg_ciu;
   }
   else{
     if(cil > 0){
       double multiplier = R::rbeta(0.78, 11) * 44;
       if(type != 0){
         ciu += 10.0;
       }
       else{
         ciu += cil + cil * multiplier;
       }
       
     }
     else{
       stop("cil is negative or zero. This should not be possible! Check if something is wrong");
     }
   }
   
   // Generation of cvx
   if(cvx_exists == 1){
     double reg_cvx = params["cvx"];
     cvx += reg_cvx;
   }
   else{
     if(sdx_exists == 1){
       double reg_sdx = params["sdx"];
       cvx += reg_sdx / (0.5 * (cil + ciu));
     }
     else{
       cvx += R::rbeta(2, 5) / 10.0;
     }
   }
   // Generation of cvy
   if(cvy_exists == 1){
     double reg_cvy = params["cvy"];
     cvy += reg_cvy;
   }
   else{
     if(sdy_exists == 1){
       double reg_sdy = params["sdy"];
       cvy += reg_sdy / (0.5 * (cil + ciu));
     }
     else{
       cvy += R::rbeta(2, 5) / 10.0;
     }
   }
   
   // Calculate mu_tau and sigma_tau
   if(dist_exists == 1){
     std::string reg_dist = params["dist"];
     dist = reg_dist;
     if(dist == "norm"){
       mu_tau += 0.5 * (cil + ciu);
       sigma_tau += 0.5 * (ciu - cil) * (1.0 / R::qnorm(0.99, 0.0, 1.0, 1, 0));
     }
     else if(dist == "lst"){
       if(df_tau_exists == 1){
         double reg_df_tau = params["df_tau"];
         df_tau = reg_df_tau;
       }
       mu_tau += 0.5 * (cil + ciu);
       sigma_tau += 0.5 * (ciu - cil) * (1.0 / R::qt(0.99, df_tau, 0, 0));
     }
     else if(dist == "lnorm"){
       mu_tau += 0.5 * log(ciu * cil);
       sigma_tau += 0.5 * log(ciu / cil) / R::qnorm(0.99, 0, 1, 1, 0);
     }
     else{
       mu_tau += 0.5 * (cil + ciu);
       sigma_tau += (1.0 / sqrt(12.0)) * (ciu - cil);
     }
   }
   
   // Check dfx and dfy
   if(error_dist_exists == 1){
     std::string reg_error_dist = params["error_dist"];
     error_dist = reg_error_dist;
     if(error_dist != "norm" and (error_dist == "lt" or error_dist == "t" or error_dist == "lst")){
       error_dist = "lt";
       if(dfx_exists == 1){
         double reg_dfx = params["dfx"];
         dfx = reg_dfx;
       }
       if(dfy_exists == 1){
         double reg_dfy = params["dfy"];
         dfy = reg_dfy;
       }
     }
   }
   
   if(eta_exists == 1){
     double reg_eta = params["eta"];
     eta += reg_eta;
   }
   else{
     eta = 1;
   }
   if(eta0_exists == 1){
     double reg_eta0 = params["eta0"];
     eta0 = reg_eta0;
   }
   else{
     eta0 = 1;
   }
   
   if(b0_exists == 1){
     double reg_b0 = params["b0"];
     b0 = reg_b0;
   }
   
   if(b1_exists == 1){
     double reg_b1 = params["b1"];
     b1 = reg_b1;
   }
   
   if(c0_exists == 1){
     double reg_c0 = params["c0"];
     c0 = reg_c0;
   }
   
   if(c1_exists == 1){
     double reg_c1 = params["c1"];
     c1 = reg_c1;
   }
   
   if(qran_exists == 1 and prop_exists == 1){
     prop_exists -= 1;
   }
   if(qpos_exists == 1){
     double reg_qpos = params["qpos"];
     qpos += reg_qpos;
   }
   else{
     qpos -= 1;
   }
   if(qran_exists == 1){
     double reg_qran = params["qran"];
     qran += reg_qran;
   }
   else{
     qran += 0;
   }
   
   if(qdir_exists == 1){
     int reg_qdir = params["qdir"];
     qdir = reg_qdir;
     if(qdir >= 1){
       qdir = 1;
     }
     else if(qdir <= -1){
       qdir = -1;
     }
     else{
       qdir = 1;
     }
   }
   else{
     int above = R::rbinom(1, 0.5);
     if(above == 1){
       qdir += 1;
     }
     else{
       qdir -= 1;
     }
   }
   
   if(prop_exists == 1){
     double reg_prop = params["prop"];
     prop += reg_prop;
   }
   else{
     prop = 0;
   }
   if(mmax_exists == 1){
     double reg_mmax = params["mmax"];
     mmax += reg_mmax;
   }
   
   if(md_method_exists == 1){
     std::string reg_md_method = params["md_method"];
     md_method = reg_md_method;
     
     if(md_method == "mar"){
       if(mar_prob_exists == 1){
         double reg_mar_prob = params["mar_prob"];
         mar_prob += reg_mar_prob;
       }
       else{
         mar_prob += 0.05;
       }
     }
     else if(md_method == "mnar"){
       if(mnar_threshold_exists == 1){
         double reg_mnar_threshold = params["mnar_threshold"];
         mnar_threshold += reg_mnar_threshold;
       }
       else{
         mnar_threshold += cil;
       }
     }
     else if(md_method == "mnar0"){
       if(mnar_threshold_exists == 1){
         double reg_mnar_threshold = params["mnar_threshold"];
         mnar_threshold += reg_mnar_threshold;
       }
     }
     else if(md_method == "marmnar"){
       if(mnar_threshold_exists == 1){
         double reg_mnar_threshold = params["mnar_threshold"];
         mnar_threshold += reg_mnar_threshold;
       }
       else{
         mnar_threshold += cil;
       }
       if(mar_prob_exists == 1){
         double reg_mar_prob = params["mar_prob"];
         mar_prob += reg_mar_prob;
       }
       else{
         mar_prob += 0.05;
       }
     }
     
     else if(md_method == "marmnar0"){
       if(mnar_threshold_exists == 1){
         double reg_mnar_threshold = params["mnar_threshold"];
         mnar_threshold += reg_mnar_threshold;
       }
       if(mar_prob_exists == 1){
         double reg_mar_prob = params["mar_prob"];
         mar_prob += reg_mar_prob;
       }
       else{
         mar_prob += 0.05;
       }
     }
     else{
       mnar_threshold -= 999999.99;
     }
   }
   else{
     mnar_threshold -= 999999.99;
   }
   
   if(obs_tau_exists == 1){
     NumericVector reg_obs_tau = params["obs_tau"];
     n = reg_obs_tau.size();
     dist = "none";
   }
   
   NumericVector unsorted_tau(n);
   NumericVector tau(n);
   
   int nR = n * R;
   IntegerVector SampleID(nR);
   IntegerVector SampleID2(n);
   IntegerVector ReplicateID(nR);
   NumericVector MP_A(nR);
   NumericVector MP_A2(n);
   NumericVector MP_B(nR);
   NumericVector MP_B2(n);
   
   double average = 0.5 * (cil + ciu);
   if(dist == "lnorm"){
     average = exp(mu_tau + pow(sigma_tau, 2.0) / 2.0);
   }
   
   double base_x = average * cvx;
   double base_y = average * cvy;
   
   double beg_sdx = base_x * eta0;
   double end_sdx = base_x * eta * eta0;
   double seg_sdx = 0;
   
   if(eta0 > 0 and eta > 0){
     seg_sdx = (end_sdx - beg_sdx) / n;
   }
   
   double beg_sdy = base_y * eta0;
   double end_sdy = base_y * eta * eta0;
   double seg_sdy = 0;
   
   if(eta0 > 0 and eta > 0){
     seg_sdy = (end_sdy - beg_sdy) / n;
   }
   
   if(dist == "none"){
     NumericVector reg_obs_tau = params["obs_tau"];
     for(int i = 0; i < n; ++i){
       unsorted_tau[i] = reg_obs_tau[i];
     }
   }
   else if(dist == "norm"){
     for(int i = 0; i < n; ++i){
       unsorted_tau[i] = R::rnorm(mu_tau, sigma_tau);
     }
   }
   else if(dist == "lst"){
     for(int i = 0; i < n; ++i){
       unsorted_tau[i] = mu_tau + sigma_tau * R::rt(df_tau);
     }
   }
   else if(dist == "lnorm"){
     for(int i = 0; i < n; ++i){
       unsorted_tau[i] = R::rlnorm(mu_tau, sigma_tau);
     }
   }
   else{
     for(int i = 0; i < n; ++i){
       unsorted_tau[i] = R::runif(cil, ciu);
     }
   }
   
   
   if(eta0_exists == 1 and eta_exists == 1){
     tau = unsorted_tau.sort();
   }
   
   for(int i = 0; i < n; ++i){
     SampleID2[i] = i + 1;
     int na_count_A = 0;
     int na_count_B = 0;
     double hetero_extra_x = i * seg_sdx;
     double hetero_extra_y = i * seg_sdy;
     double sdx = beg_sdx + hetero_extra_x;
     double sdy = beg_sdy + hetero_extra_y;
     
     if(sdx <= 0 or sdy <= 0){
       stop("Your choices of eta and eta0 resulted in negative standard deviations. Calculations are terminated");
     }
     
     if(eta0_exists == 0 or eta_exists == 0){
       tau[i] = unsorted_tau[i];
     }
     
     double limit = 0;
     double relocating_magnitude = 0;
     int relocate_sample_i = 0;
     
     // if prop exists, and CS is dins-affected, randomly select relocation magnitude from from Beta(2, 2) * mmax
     if(prop_exists == 1){
       relocate_sample_i += R::rbinom(1, prop);
       double sign_relocated_sample_i = R::runif(0, 1);
       if(sign_relocated_sample_i < 0.5){
         sign_relocated_sample_i = -1;
       }
       else{
         sign_relocated_sample_i = 1;
       }
       relocating_magnitude += R::rbeta(2, 2) * sign_relocated_sample_i * mmax;
     }
     else if(qran_exists == 1){
       if(qpos == 0){
         if(qlim_exists == 1){
           double reg_qlim = params["qlim"];
           limit += reg_qlim;
         }
         else{
           limit += cil + qran * (ciu - cil);
         }
         if(qnor_exists == 1){
           double reg_qnor = params["qnor"];
           qnor = reg_qnor;
         }
         else{
           qnor = limit - cil;
         }
         qlim = limit;
         if(tau[i] <= limit){
           relocate_sample_i = relocate_sample_i + 1;
           double rel_diff = R::pbeta((limit - tau[i]) / qnor, 2, 2, 1, 0) / 2.0;
           relocating_magnitude = 2 * rel_diff * mmax * qdir;
         }
       }
       else if(qpos == 1){
         if(qlim_exists == 1){
           double reg_qlim = params["qlim"];
           limit += reg_qlim;
         }
         else{
           limit += ciu - qran * (ciu - cil);
         }
         if(qnor_exists == 1){
           double reg_qnor = params["qnor"];
           qnor = reg_qnor;
         }
         else{
           qnor = ciu - limit;
         }
         qlim = limit;
         if(tau[i] >= limit){
           relocate_sample_i = relocate_sample_i + 1;
           double rel_diff = R::pbeta((tau[i] - limit) / qnor, 2, 2, 1, 0) / 2.0;
           relocating_magnitude = 2 * rel_diff * mmax * qdir;
         }
       }
     }
     int lower = i * R;
     int upper = R * (1 + i) - 1;
     IntegerVector idx = Rcpp::Range(lower, upper);
     if(error_dist == "norm"){
       if(g_exists == 1){
         MP_B[idx] = custom_function(MP_B[idx], params["g"]);
       }
       for(int r = 0; r < R; ++r){
         SampleID[idx[r]] = i + 1;
         ReplicateID[idx[r]] = r + 1;
         MP_B[idx[r]] = c0 + tau[i] * c1;
         
         if(g_exists == 1){
           MP_A[idx[r]] = MP_B[idx[r]] + R::rnorm(0, sdy);
         }
         else if(type == 1){
           MP_A[idx[r]] = MP_B[idx[r]] + 0.9 * sin(0.4 * pow(MP_B[idx[r]], 1.06)) + R::rnorm(0, sdy);
         }
         else if(type == 2){
           MP_A[idx[r]] = MP_B[idx[r]] + 0.05 * exp(0.16 * pow(MP_B[idx[r]], 1.35)) + R::rnorm(0, sdy);
         }
         else if(type == 3){
           MP_A[idx[r]] = MP_B[idx[r]] - exp(-0.5 * pow((MP_B[idx[r]] - 1.5)/2.0, 2)) + R::rnorm(0, sdy);
         }
         else{
           MP_A[idx[r]] = b0 + MP_B[idx[r]] * b1 + R::rnorm(0, sdy) + relocate_sample_i * relocating_magnitude * sqrt(pow(sdx, 2) + pow(sdy, 2));
         }
         MP_B[idx[r]] = MP_B[idx[r]] + R::rnorm(0, sdx);
         if(!AR){
           MP_A2[i] += MP_A[idx[r]];
           MP_B2[i] += MP_B[idx[r]];
         }
         if(md_method == "mar" or md_method == "mnar" or md_method == "mnar0" or md_method == "marmnar" or md_method == "marmnar0"){
           bool mar_r_A = R::runif(0, 1) < mar_prob;
           bool mar_r_B = R::runif(0, 1) < mar_prob;
           bool mnar_r_A = MP_A[idx[r]] < mnar_threshold;
           bool mnar_r_B = MP_B[idx[r]] < mnar_threshold;
           if(mar_r_A or mnar_r_A){
             MP_A[idx[r]] = NA_REAL;
             if(!AR){
               na_count_A += 1;
               MP_A2[i] -= MP_A[idx[r]];
             }
             
           }
           if(mar_r_B or mnar_r_B){
             MP_B[idx[r]] = NA_REAL;
             if(!AR){
               na_count_B += 1;
               MP_B2[i] -= MP_B[idx[r]];
             }
           }
         }
         else{
           if(MP_A[idx[r]] < 0){
             MP_A[idx[r]] = MP_A[idx[r]] * (-1);
           }
           if(MP_B[idx[r]] < 0){
             MP_B[idx[r]] = MP_B[idx[r]] * (-1);
           }
         }
       }
     }
     else if(error_dist == "lt"){
       if(g_exists == 1){
         MP_B[idx] = custom_function(MP_B[idx], params["g"]);
       }
       for(int r = 0; r < R; ++r){
         SampleID[idx[r]] = i + 1;
         ReplicateID[idx[r]] = r + 1;
         MP_B[idx[r]] = c0 + tau[i] * c1;
         if(g_exists == 1){
           MP_A[idx[r]] = MP_B[idx[r]] + R::rnorm(0, sdy);
         }
         else if(type == 1){
           MP_A[idx[r]] = MP_B[idx[r]] + 0.9 * sin(0.4 * pow(MP_B[idx[r]], 1.06)) + sdy * R::rt(dfy);
         }
         else if(type == 2){
           MP_A[idx[r]] = MP_B[idx[r]] + 0.05 * exp(0.16 * pow(MP_B[idx[r]], 1.35)) + sdy * R::rt(dfy);
         }
         else if(type == 3){
           MP_A[idx[r]] = MP_B[idx[r]] - exp(-0.5 * pow((MP_B[idx[r]] - 1.5)/2.0, 2)) + sdy * R::rt(dfy);
         }
         else{
           MP_A[idx[r]] = b0 + MP_B[idx[r]] * b1 + sdy * R::rt(dfy) + relocate_sample_i * relocating_magnitude * sqrt(pow(sdx, 2) + pow(sdy, 2));
         }
         MP_B[idx[r]] = MP_B[idx[r]] + sdx * R::rt(dfx);
         if(!AR){
           MP_A2[i] += MP_A[idx[r]];
           MP_B2[i] += MP_B[idx[r]];
         }
         if(md_method == "mar" or md_method == "mnar" or md_method == "mnar0" or md_method == "marmnar" or md_method == "marmnar0"){
           bool mar_r_A = R::runif(0, 1) < mar_prob;
           bool mar_r_B = R::runif(0, 1) < mar_prob;
           bool mnar_r_A = MP_A[idx[r]] < mnar_threshold;
           bool mnar_r_B = MP_B[idx[r]] < mnar_threshold;
           if(mar_r_A or mnar_r_A){
             MP_A[idx[r]] = NA_REAL;
           }
           if(mar_r_B or mnar_r_B){
             MP_B[idx[r]] = NA_REAL;
           }
         }
         else{
           if(MP_A[idx[r]] < 0){
             MP_A[idx[r]] = MP_A[idx[r]] * (-1);
           }
           if(MP_B[idx[r]] < 0){
             MP_B[idx[r]] = MP_B[idx[r]] * (-1);
           }
         }
       }
     }
     if(!AR){
       if(R - na_count_A > 0){
         MP_A2[i] = MP_A2[i] / (R - na_count_A);
       }
       else{
         MP_A2[i] = NA_REAL;
       }
       if(R - na_count_A > 0){
         MP_B2[i] = MP_B2[i] / (R - na_count_B);
       }
       else{
         MP_B2[i] = NA_REAL;
       }
     }
   }
   
   List used_parameters = List::create(Named("n") = n,
                                       Named("R") = R,
                                       Named("cvx") = cvx,
                                       Named("cvy") = cvy,
                                       Named("sdx") = base_x,
                                       Named("sdy") = base_y,
                                       Named("cil") = cil,
                                       Named("ciu") = ciu,
                                       Named("dfx") = dfx,
                                       Named("dfy") = dfy,
                                       Named("qpos") = qpos,
                                       Named("qran") = qran,
                                       Named("qdir") = qdir,
                                       Named("qlim") = qlim,
                                       Named("qnor") = qnor,
                                       Named("mmax") = mmax,
                                       Named("md_method") = md_method,
                                       Named("mar_prob") = mar_prob,
                                       Named("mnar_threshold") = mnar_threshold);
   
   if(!AR){
     if(include_parameters){
       List sim_data = List::create(Named("SampleID") = SampleID2,
                                    Named("MP_A") = round(MP_A2, 6),
                                    Named("MP_B") = round(MP_B2, 6));
       List out = List::create(Named("simulated_data") = sim_data,
                               Named("parameters") = used_parameters);
       return out;
     }
     List out = List::create(Named("SampleID") = SampleID2,
                             Named("MP_A") = round(MP_A2, 6),
                             Named("MP_B") = round(MP_B2, 6));
     return out;
   }
   if(include_parameters){
     List sim_data = List::create(Named("SampleID") = SampleID,
                                  Named("ReplicateID") = ReplicateID,
                                  Named("MP_A") = round(MP_A, 6),
                                  Named("MP_B") = round(MP_B, 6));
     List out = List::create(Named("simulated_data") = sim_data,
                             Named("parameters") = used_parameters);
     return out;
   }
   List out = List::create(Named("SampleID") = SampleID,
                           Named("ReplicateID") = ReplicateID,
                           Named("MP_A") = round(MP_A, 6),
                           Named("MP_B") = round(MP_B, 6));
   return out;
 }

// Helper function to call a custom R function safely (expects scalar in/out)
// [[Rcpp::plugins(cpp11)]]
inline double call_custom_g(double input_val, const Function& func) {
  NumericVector input_vec = NumericVector::create(input_val);
  NumericVector output_vec = func(input_vec);
  if (output_vec.size() != 1) {
    Rcpp::stop("Custom function 'g' must return a single numeric value for a single numeric input.");
  }
  if (NumericVector::is_na(output_vec[0])) {
    Rcpp::warning("Custom function 'g' returned NA for input %f", input_val);
    return NA_REAL;
  }
  return output_vec[0];
}


//' @title Simulation of External Quality Assessment (EQA) Data
//' @name sim_eqa_data2
//'
//' @param parameters A \code{list} of parameters and their value. These values
//'        are used to simulate the EQA data. Optional parameters that you may
//'        include into the list are:
//' \itemize{
//'   \item \code{n: } The number of samples. Defaults to 25.
//'   \item \code{R: } The number of replicates on each sample. Must be >= 1. Defaults to 3.
//'   \item \code{cvx: } The base repeatability coefficient of variation for IVD-MD
//'                      \code{MP_B}. Defaults to 0.01.
//'   \item \code{cvy: } The base repeatability coefficient of variation for IVD-MD
//'                      \code{MP_A}. Defaults to 0.01.
//'   \item \code{cil: } The lower bound of the concentration interval. Defaults to 2.0.
//'   \item \code{ciu: } The upper bound of the concentration interval. Must be > \code{cil}. Defaults to 10.0.
//'   \item \code{dist: } The distribution to simulate latent variables \eqn{\tau_i} from.
//'                       Possible choices include
//'                       \code{"unif"} (uniform distribution, default),
//'                       \code{"norm"} (normal distribution),
//'                       \code{"lst"} (location-scale student t-distribution),
//'                       \code{"lnorm"} (log-normal distribution).
//'   \item \code{df_tau: } The degrees of freedom for the \code{lst}
//'                         distribution if \code{dist = "lst"}.
//'                         Defaults to 5.0 if not otherwise specified. Must be > 0.
//'   \item \code{eta: } The heteroscedasticity factor \eqn{\eta}. Controls how strongly SD depends on concentration. See Details. Defaults to 0.
//'   \item \code{eta0: } The base heteroscedasticity factor \eqn{\eta_0}. See Details. Defaults to 1.0.
//'   \item \code{qpos: } Position of systematic differences in non-selectivity.
//'                       \code{0} signify lower concentration interval and
//'                       \code{1} signify upper. Only used if \code{qran} is specified and \code{prop} is not.
//'   \item \code{qran: } Interquantile range (proportion of total range, 0 to 1) where systematic differences in
//'                       non-selectivity should have its effect. E.g., 0.3 means the lower or upper 30% of the range. Only used if \code{prop} is not specified.
//'   \item \code{prop: } Average proportion (0 to 1) of clinical samples affected by
//'                       random differences in non-selectivity. If specified, overrides \code{qpos}/\code{qran}. Defaults to 0.
//'   \item \code{mmax: } The maximum relocation magnitude multiplier,
//'                       \eqn{m_{\max}}. This number is a multiplier of the
//'                       *base* standard deviation \code{sqrt(base_sd_x^2 + base_sd_y^2)}. Only effective
//'                       if \code{prop > 0} or \code{qran > 0}. Defaults to 0.
//'   \item \code{b0: } For systematic linear bias between IVD-MDs (intercept). Defaults to 0.
//'   \item \code{b1: } For systematic linear bias between IVD-MDs (slope). Defaults to 1.
//'   \item \code{c0: } Intercept for linear transformation from latent \eqn{\xi_i} to \eqn{\tau_i}. Defaults to 0. (Often ignored).
//'   \item \code{c1: } Slope for linear transformation from latent \eqn{\xi_i} to \eqn{\tau_i}. Defaults to 1. (Often ignored).
//'   \item \code{error_dist: } The distribution to simulate measurement error
//'                             components (\eqn{h_{ir}, v_{ir}}) from. Possible choices include \code{"norm"}
//'                             (normal distribution, default) and
//'                             \code{"lt"} (location-scale student t-distribution).
//'   \item \code{dfx: } The degrees of freedom for the measurement error
//'                      components in IVD-MD B if \code{error_dist = "lt"}.
//'                      Defaults to 5.0 if not otherwise specified. Must be > 0.
//'   \item \code{dfy: } The degrees of freedom for the measurement error
//'                      components in IVD-MD A if \code{error_dist = "lt"}.
//'                      Defaults to 5.0 if not otherwise specified. Must be > 0.
//'   \item \code{md_method: } Method for simulating missing data. Possible
//'                            choices include \code{"none"} (no missing data, default),
//'                            \code{"mar"} (missing at random),
//'                            \code{"mnar"} (missing not at random based on value vs threshold),
//'                            \code{"mnar0"} (missing not at random if value < threshold, treats negative values as missing if threshold <= 0),
//'                            \code{"marmnar"} (combination of MAR and MNAR),
//'                            \code{"marmnar0"} (combination of MAR and MNAR0).
//'   \item \code{mar_prob: } The probability (value between 0 and 1) of having
//'                           a measurement missing at random. Only relevant if
//'                           \code{md_method} includes \code{"mar"}. Defaults to 0.05 if needed.
//'   \item \code{mnar_threshold: } The threshold (\code{double}) below which
//'                                 a measurement is considered missing for MNAR.
//'                                 Only relevant if \code{md_method} includes \code{"mnar"} or \code{"mnar0"}.
//'                                 Defaults to \code{cil} if \code{md_method = "mnar"} or \code{"marmnar"}.
//'                                 Defaults to 0.0 if \code{md_method = "mnar0"} or \code{"marmnar0"}.
//'   \item \code{g: } A custom R function defining the relationship \eqn{g(\tau_i) = E[y_{ir} | \tau_i]}.
//'                    It should take a single numeric value (\eqn{\tau_i}) and return a single numeric value.
//'                    Overrides \code{type}, \code{b0}, \code{b1}, and non-selectivity parameters (\code{prop}, \code{qpos}, \code{qran}, \code{mmax}).
//'                    Using this can be slower than built-in types.
//'   \item \code{obs_tau: } A \code{numeric} vector of pre-defined \eqn{\tau_i} values.
//'                          If specified, \code{n} becomes \code{length(obs_tau)}, and \code{dist}, \code{cil}, \code{ciu}, \code{df_tau} are ignored.
//' }
//' @param type A \code{integer} specifying built-in functions for \eqn{g(\tau_i)}. Ignored if \code{g} is provided in \code{parameters}.
//'        \itemize{
//'           \item \code{type = 0: } Linear relationship \eqn{g(\tau_i) = \beta_0 + \beta_1 \tau_i + \mathrm{nonselectivity}}. (Default)
//'           \item \code{type = 1: } Non-linear \eqn{g(\tau_i) = \tau_i + 0.90 \cdot \sin(0.40 \cdot \tau_i^{1.06})}
//'           \item \code{type = 2: } Non-linear \eqn{g(\tau_i) = \tau_i + 0.05 \cdot \exp(0.16 \cdot \tau_i^{1.35})}
//'           \item \code{type = 3: } Non-linear \eqn{g(\tau_i) = \tau_i - \exp[-0.125 \cdot (\tau_i - 1.50)^2]}
//'        }
//' @param AR A \code{logical} value. If \code{TRUE}, data includes all replicates. If \code{FALSE}, the mean of replicates (MOR) for each sample is returned (handling NAs appropriately).
//' @param include_parameters A \code{logical} value. If \code{TRUE}, the parameters used in the simulation are returned as a list element named \code{parameters}. Defaults to \code{FALSE}.
//'
//' @description
//' Simulates method comparison data based on specified parameters, allowing for various complexities like non-linear relationships, heteroscedasticity, non-selectivity, and missing data.
//'
//' @details
//' The simulation generates measurement pairs \eqn{(x_{ir}, y_{ir})} for sample \eqn{i} and replicate \eqn{r}, based on the model:
//' \deqn{x_{ir} = \tau_i + h_{ir}}
//' \deqn{y_{ir} = g(\tau_i) + v_{ir}}
//' where \eqn{\tau_i} represents the true (latent) value for sample \eqn{i} for method B (often assumed \eqn{\tau_i = c_0 + c_1 \xi_i}), \eqn{g(\tau_i)} is the expected value for method A given \eqn{\tau_i}, and \eqn{h_{ir}} and \eqn{v_{ir}} are measurement errors.
//'
//' \strong{True Values (\eqn{\tau_i}):} Generated based on \code{dist}, \code{cil}, \code{ciu}, unless \code{obs_tau} is provided.
//'
//' \strong{Relationship (\eqn{g}):} Defined by \code{type} or a custom R function \code{g}. For \code{type = 0}, it incorporates linear bias (\code{b0}, \code{b1}) and potential non-selectivity effects (\code{prop}, \code{qpos}, \code{qran}, \code{mmax}). Non-selectivity effects add an offset for specific samples.
//'
//' \strong{Measurement Errors (\eqn{h_{ir}, v_{ir}}):} Generated based on \code{error_dist} (\code{"norm"} or \code{"lt"}) with standard deviations derived from \code{cvx}, \code{cvy} and potentially modified by heteroscedasticity parameters \code{eta} and \code{eta0}. The model used is \eqn{SD(\tau_i) = \mathrm{base\_sd} \times (\eta_0 + (\tau_i / \mathrm{average\_tau})^\eta)}, where \eqn{\mathrm{base\_sd}} is \code{average_tau * cv}, and \code{average_tau} is \code{0.5*(cil+ciu)}. If \code{eta = 0} and \code{eta0 = 1}, this results in constant CV.
//'
//' \strong{Non-Selectivity:} If \code{prop > 0}, a random proportion of samples have their \eqn{g(\tau_i)} value shifted by a random amount up to \code{mmax} times the base combined SD. If \code{qran > 0} (and \code{prop = 0}), samples in the specified quantile range (\code{qpos}, \code{qran}) have their \eqn{g(\tau_i)} shifted systematically based on their position within that range, up to \code{mmax}.
//'
//' \strong{Missing Data:} Can be introduced using \code{md_method}. \code{"mar"} introduces NAs randomly. \code{"mnar"} introduces NAs if the simulated value falls below \code{mnar_threshold}. \code{"mnar0"} is similar to \code{"mnar"} but defaults the threshold to 0 and handles potentially negative simulated values as missing below this threshold.
//'
//' \strong{Negative Values:} If \code{md_method} is not \code{"mnar0"}, any simulated negative measurement values \eqn{x_{ir}} or \eqn{y_{ir}} are clamped at 0.0. If \code{md_method = "mnar0"}, values below the \code{mnar_threshold} (defaulting to 0) are set to \code{NA_REAL}. Choose parameters (\code{cil}, \code{cvx}, \code{cvy}, \code{mmax}) carefully to avoid excessive negative values if they are physically implausible.
//'
//' \strong{Output Conversion:} The returned list can be efficiently converted to other formats using \code{as.data.frame()}, \code{as.data.table::data.table()}, \code{tibble::as_tibble()}, or \code{data.table::setDT()}.
//'
//' @return A \code{list}.
//' If \code{include_parameters = FALSE}:
//' Contains elements \code{SampleID}, \code{MP_A}, \code{MP_B}. If \code{AR = TRUE}, also includes \code{ReplicateID}.
//' If \code{include_parameters = TRUE}:
//' A list containing two elements:
//' \itemize{
//'   \item \code{simulated_data:} The list described above.
//'   \item \code{parameters:} A list of the parameter values actually used in the simulation.
//' }
//'
//' @examples
//' # Load package (assuming it's built and loaded)
//' # library(yourPackageName) # Replace with your package name
//'
//' # Default simulation (n=25, R=3, linear, const CV)
//' default_sim <- sim_eqa_data()
//' print(head(as.data.frame(default_sim)))
//'
//' # Simulation with heteroscedasticity and t-distributed errors
//' params_hetero_t <- list(n = 50, R = 2, cvx = 0.02, cvy = 0.03,
//'                         cil = 10, ciu = 100, eta = 0.5, eta0 = 0.5,
//'                         error_dist = "lt", dfx = 4, dfy = 6)
//' sim_hetero_t <- sim_eqa_data(params_hetero_t, AR = TRUE)
//' # plot(sim_hetero_t$MP_B, sim_hetero_t$MP_A) # Requires conversion first
//'
//' # Simulation with random non-selectivity (type 0) and MAR missing data
//' params_nonselect_mar <- list(n = 40, R = 4, cvx = 0.01, cvy = 0.01,
//'                              cil = 5, ciu = 50, prop = 0.15, mmax = 4,
//'                              md_method = "mar", mar_prob = 0.1)
//' sim_nonselect_mar <- sim_eqa_data(params_nonselect_mar, type = 0, AR = TRUE,
//'                                   include_parameters = TRUE)
//' print(head(as.data.frame(sim_nonselect_mar$simulated_data)))
//' print(sim_nonselect_mar$parameters)
//'
//' # Simulation with systematic non-selectivity (upper range) and MNAR missing data
//' params_syst_nonselect_mnar <- list(n = 30, R = 2, cvx = 0.015, cvy = 0.015,
//'                                    cil = 2, ciu = 20, qpos = 1, qran = 0.25, mmax = 5,
//'                                    md_method = "mnar") # threshold defaults to cil=2
//' sim_syst_nonselect_mnar <- sim_eqa_data(params_syst_nonselect_mnar, type = 0, AR = FALSE)
//' print(head(as.data.frame(sim_syst_nonselect_mnar)))
//'
//' # Simulation using a custom non-linear function g and observed tau
//' my_g <- function(tau) { tau + 0.1 * tau^2 }
//' true_taus <- c(1, 3, 5, 7, 9, 11, 13, 15)
//' params_custom_g <- list(R = 5, cvx = 0.005, cvy = 0.008, obs_tau = true_taus, g = my_g)
//' sim_custom_g <- sim_eqa_data(params_custom_g, AR = TRUE)
//' print(head(as.data.frame(sim_custom_g)))
//'
// [[Rcpp::export]]
List sim_eqa_data2(Nullable<List> parameters = R_NilValue,
                  int type = 0,
                  bool AR = false, // Changed default to TRUE as it seems more common
                  bool include_parameters = false) {
 
 // --- Parameter Initialization and Validation ---
 List params;
 if (parameters.isNotNull()) {
   params = as<List>(parameters);
 }
 
 // Use containsElementNamed for safer checking
 int n = params.containsElementNamed("n") ? as<int>(params["n"]) : 25;
 int R = params.containsElementNamed("R") ? as<int>(params["R"]) : 3;
 double cvx = params.containsElementNamed("cvx") ? as<double>(params["cvx"]) : 0.01;
 double cvy = params.containsElementNamed("cvy") ? as<double>(params["cvy"]) : 0.01;
 double cve = params.containsElementNamed("cve") ? as<double>(params["cve"]) : 0.00;
 double zeta0 = params.containsElementNamed("zeta0") ? as<double>(params["zeta0"]) : 1.00;
 double cil = params.containsElementNamed("cil") ? as<double>(params["cil"]) : 2.0;
 double ciu = params.containsElementNamed("ciu") ? as<double>(params["ciu"]) : 10.0;
 std::string dist = params.containsElementNamed("dist") ? as<std::string>(params["dist"]) : "unif";
 double df_tau = params.containsElementNamed("df_tau") ? as<double>(params["df_tau"]) : 5.0;
 double eta = params.containsElementNamed("eta") ? as<double>(params["eta"]) : 0.0;
 double eta0 = params.containsElementNamed("eta0") ? as<double>(params["eta0"]) : 1.0;
 double prop = params.containsElementNamed("prop") ? as<double>(params["prop"]) : 0.0;
 int qpos = params.containsElementNamed("qpos") ? as<int>(params["qpos"]) : 0; // Default 0 (lower), used if qran > 0 & prop = 0
 double qran = params.containsElementNamed("qran") ? as<double>(params["qran"]) : 0.0;
 double mmax = params.containsElementNamed("mmax") ? as<double>(params["mmax"]) : 0.0;
 double b0 = params.containsElementNamed("b0") ? as<double>(params["b0"]) : 0.0;
 double b1 = params.containsElementNamed("b1") ? as<double>(params["b1"]) : 1.0;
 double c0 = params.containsElementNamed("c0") ? as<double>(params["c0"]) : 0.0;
 double c1 = params.containsElementNamed("c1") ? as<double>(params["c1"]) : 1.0;
 std::string error_dist = params.containsElementNamed("error_dist") ? as<std::string>(params["error_dist"]) : "norm";
 double dfx = params.containsElementNamed("dfx") ? as<double>(params["dfx"]) : 5.0;
 double dfy = params.containsElementNamed("dfy") ? as<double>(params["dfy"]) : 5.0;
 std::string md_method = params.containsElementNamed("md_method") ? as<std::string>(params["md_method"]) : "none";
 double mar_prob = params.containsElementNamed("mar_prob") ? as<double>(params["mar_prob"]) : 0.05; // Default if needed later
 double mnar_threshold = params.containsElementNamed("mnar_threshold") ? as<double>(params["mnar_threshold"]) : R_PosInf; // Default placeholder
 bool obs_tau_exists = params.containsElementNamed("obs_tau");
 bool g_exists = params.containsElementNamed("g");

 // --- Input Validation ---
 if (R < 1) Rcpp::stop("R (number of replicates) must be >= 1.");
 if (cvx < 0) Rcpp::stop("cvx must be non-negative.");
 if (cvy < 0) Rcpp::stop("cvy must be non-negative.");
 if (ciu <= cil) Rcpp::stop("ciu must be greater than cil.");
 if ((dist == "lst" || dist == "t") && df_tau <= 0) Rcpp::stop("df_tau must be positive for lst distribution.");
 if (error_dist == "lt" && dfx <= 0) Rcpp::stop("dfx must be positive for lt error distribution.");
 if (error_dist == "lt" && dfy <= 0) Rcpp::stop("dfy must be positive for lt error distribution.");
 if (prop < 0 || prop > 1) Rcpp::stop("prop must be between 0 and 1.");
 if (qran < 0 || qran > 1) Rcpp::stop("qran must be between 0 and 1.");
 if (md_method.find("mar") != std::string::npos && (mar_prob < 0 || mar_prob > 1)) Rcpp::stop("mar_prob must be between 0 and 1.");
 
 
 // Adjust defaults/logic based on inputs
 NumericVector obs_tau_values;
 if (obs_tau_exists) {
   obs_tau_values = as<NumericVector>(params["obs_tau"]);
   n = obs_tau_values.size();
   dist = "observed"; // Mark dist as observed
   if (n == 0) Rcpp::stop("obs_tau is empty.");
 }
 
 bool use_mar = (md_method.find("mar") != std::string::npos);
 bool use_mnar = (md_method.find("mnar") != std::string::npos); // Catches mnar, mnar0, marmnar, marmnar0
 bool use_mnar0_logic = (md_method.find("mnar0") != std::string::npos); // Catches mnar0, marmnar0
 
 if (use_mnar && !params.containsElementNamed("mnar_threshold")) {
   mnar_threshold = use_mnar0_logic ? 0.0 : cil;
 } else if (!use_mnar) {
   // If MNAR is not used, set threshold effectively to negative infinity
   // so the check `value < mnar_threshold` is always false.
   mnar_threshold = R_NegInf;
 }
 // If MAR is needed but probability wasn't given, use default
 if (use_mar && !params.containsElementNamed("mar_prob")) {
   // mar_prob already defaulted to 0.05 above
 } else if (!use_mar) {
   mar_prob = 0.0; // Ensure it's 0 if MAR not selected
 }
 
 
 // Determine if random non-selectivity (prop) or systematic (qran) is used
 bool use_prop_nonselectivity = (prop > 0 && mmax != 0);
 bool use_qran_nonselectivity = (!use_prop_nonselectivity && qran > 0 && mmax != 0);
 
 // --- Simulate True Values (tau) ---
 NumericVector tau(n);
 double mu_tau = 0.0, sigma_tau = 0.0; // Parameters for generating tau if needed
 
 if (dist == "observed") {
   tau = clone(obs_tau_values); // Use provided values
 } else if (dist == "unif") {
   tau = runif(n, cil, ciu);
 } else {
   // Calculate parameters needed for norm, lst, lnorm
   if (dist == "norm") {
     mu_tau = 0.5 * (cil + ciu);
     sigma_tau = 0.5 * (ciu - cil) / R::qnorm(0.99, 0.0, 1.0, true, false);
     tau = rnorm(n, mu_tau, sigma_tau);
   } else if (dist == "lst" || dist == "t") { // Treat t as lst
     dist = "lst"; // Standardize name
     mu_tau = 0.5 * (cil + ciu);
     sigma_tau = 0.5 * (ciu - cil) / R::qt(0.99, df_tau, true, false);
     NumericVector rt_vals = rt(n, df_tau);
     tau = mu_tau + sigma_tau * rt_vals;
   } else if (dist == "lnorm") {
     // Ensure cil > 0 for lnorm
     if (cil <= 0) Rcpp::stop("cil must be > 0 for dist = 'lnorm'.");
     mu_tau = 0.5 * (log(ciu) + log(cil)); // Mean of log(tau)
     sigma_tau = 0.5 * (log(ciu) - log(cil)) / R::qnorm(0.99, 0.0, 1.0, true, false); // SD of log(tau)
     tau = rlnorm(n, mu_tau, sigma_tau);
   } else {
     Rcpp::stop("Unknown distribution specified for 'dist'.");
   }
   // Optional: Clamp generated values to bounds? (Could distort distribution)
   // tau = pmax(cil, pmin(ciu, tau)); // Consider if necessary
 }
 
 // --- Prepare for Simulation ---
 int nR = n * R;
 IntegerVector SampleID(nR);
 IntegerVector ReplicateID(nR);
 NumericVector MP_A(nR);
 NumericVector MP_B(nR);
 
 IntegerVector SampleID_mor(n); // For AR = FALSE
 NumericVector MP_A_mor(n);     // For AR = FALSE
 NumericVector MP_B_mor(n);     // For AR = FALSE
 
 double average_tau = 0.5 * (cil + ciu); // Reference for heteroscedasticity
 if (dist == "lnorm" && average_tau > 0) { // Use geometric mean center if log-normal approx
   average_tau = exp(mu_tau); // More representative center on original scale
 } else if (dist == "observed"){
   if (n > 0) average_tau = Rcpp::mean(tau); // Use empirical mean if observed
 }
 if (average_tau <= 0) average_tau = 1.0; // Avoid division by zero/negatives
 
 double base_x = average_tau * cvx; // Base SD for method B
 double base_y = average_tau * cvy; // Base SD for method A
 double base_e = average_tau * cve;
 double base_sd_combined = std::sqrt(base_x * base_x + base_y * base_y);
 if (base_sd_combined == 0 && (use_prop_nonselectivity || use_qran_nonselectivity)) {
   Rcpp::warning("mmax specified but base SDs are zero. Non-selectivity effect will be zero.");
   base_sd_combined = 1.0; // Avoid division by zero later
 }
 
 
 // Parameters for systematic non-selectivity (qran)
 double qlim_lower = cil;
 double qlim_upper = ciu;
 double qrange_width = ciu - cil;
 double qnorm_factor = 1.0; // Normalization factor for pbeta calculation
 int qdir = 1; // Direction of shift (+1 or -1)
 
 if (use_qran_nonselectivity) {
   qrange_width = (ciu - cil) * qran;
   if (qpos == 0) { // Lower end
     qlim_upper = cil + qrange_width;
     qnorm_factor = qrange_width;
   } else { // Upper end (qpos == 1 or other)
     qlim_lower = ciu - qrange_width;
     qnorm_factor = qrange_width;
   }
   if (qnorm_factor <= 0) qnorm_factor = 1.0; // Avoid division by zero
   
   // Determine direction randomly if not specified via 'qdir' param (undocumented)
   if (!params.containsElementNamed("qdir")) {
     qdir = (R::runif(0, 1) < 0.5) ? -1 : 1;
   } else {
     qdir = as<int>(params["qdir"]);
     qdir = (qdir >= 0) ? 1 : -1; // Ensure +1 or -1
   }
 }
 
 
 // --- Main Simulation Loop ---
 for (int i = 0; i < n; ++i) {
   double tau_i = tau[i];
   double true_b = c0 + tau_i * c1; // True value underlying method B for sample i
   
   // Calculate sample-specific SDs (Heteroscedasticity)
   // Model: sd = base_sd * (eta0 + (tau_i / average_tau)^eta)
   // Note: If average_tau is <= 0, power calculation is skipped.
   double sdx = base_x;
   double sdy = base_y;
   double sde = base_e;
   if (zeta0 > 1.0) {
     double vare = (pow(sdy, 2.0) + pow(b1, 2.0) * pow(sdx, 2.0)) * (zeta0 - 1);
     sde = sqrt(vare);
   }
   else if (zeta0 < 1.0) {
     sde = 0.0;
   }
   
   if (average_tau > 0 && (eta != 0.0 || eta0 != 1.0)) {
     double rel_conc = tau_i / average_tau;
     if (rel_conc < 0) rel_conc = 0; // Avoid issues with negative tau in power
     double hetero_factor = eta0 + ((eta != 0.0) ? pow(rel_conc, eta) : 0.0);
     sdx = base_x * hetero_factor;
     sdy = base_y * hetero_factor;
   }
   // Ensure SDs are not negative (can happen with extreme eta/eta0)
   if (sdx < 0) sdx = 0;
   if (sdy < 0) sdy = 0;
   
   
   // Calculate non-selectivity offset for this sample (only if type=0 and not custom g)
   double nonselectivity_offset = 0.0;
   bool is_affected = false; // Flag if sample i is affected by non-selectivity
   
   if (!g_exists && type == 0) {
     if (use_prop_nonselectivity) {
       if (R::runif(0, 1) < prop) {
         is_affected = true;
         double sign_reloc = (R::runif(0, 1) < 0.5) ? -1.0 : 1.0;
         // Use Beta(2,2) for magnitude distribution centered at 0.5 * max
         nonselectivity_offset = R::rbeta(2.0, 2.0) * sign_reloc * mmax * base_sd_combined;
       }
     } else if (use_qran_nonselectivity) {
       double relative_pos = 0.0; // Position within the affected quantile range (0 to 1)
       if (qpos == 0 && tau_i < qlim_upper) { // Affected in lower range
         is_affected = true;
         // Calculate relative position (0=cil, 1=qlim_upper)
         relative_pos = (tau_i - cil) / qnorm_factor;
         
       } else if (qpos != 0 && tau_i > qlim_lower) { // Affected in upper range
         is_affected = true;
         // Calculate relative position (0=qlim_lower, 1=ciu)
         relative_pos = (tau_i - qlim_lower) / qnorm_factor;
       }
       
       if (is_affected) {
         relative_pos = std::max(0.0, std::min(1.0, relative_pos)); // Clamp to [0, 1]
         // Use pbeta-like shape for magnitude based on position
         // Centered shape using 2*pbeta - 1 might be complex, use linear for now?
         // Or simpler: scale magnitude linearly with position?
         // Let's try: magnitude increases towards the *outer* edge of the range
         double magnitude_factor;
         if (qpos == 0) { // Magnitude increases as tau decreases towards cil
           magnitude_factor = R::pbeta(1.0 - relative_pos, 2.0, 2.0, true, false); // Use pbeta for shape
         } else { // Magnitude increases as tau increases towards ciu
           magnitude_factor = R::pbeta(relative_pos, 2.0, 2.0, true, false);
         }
         // Scale by mmax and base SD, apply direction
         nonselectivity_offset = magnitude_factor * mmax * base_sd_combined * static_cast<double>(qdir);
       }
     }
   } // end if (!g_exists && type == 0)
   
   // --- Replicate Loop ---
   double sum_a = 0.0, sum_b = 0.0;
   int count_a = 0, count_b = 0;
   double random_dins = R::rnorm(0.0, sde);
   
   for (int r = 0; r < R; ++r) {
     int current_idx = i * R + r;
     SampleID[current_idx] = i + 1;
     ReplicateID[current_idx] = r + 1;
     
     // Calculate expected value for A based on g(tau) or type
     double expected_a;
     if (g_exists) {
       Function actual_g_func = as<Function>(params["g"]); 
       expected_a = call_custom_g(true_b, actual_g_func) + random_dins; // Uses true_b as input to g
     } else if (type == 1) {
       expected_a = true_b + 0.9 * sin(0.4 * pow(fmax(0.0, true_b), 1.06)) + random_dins; // Use fmax to avoid pow issues if true_b<0
     } else if (type == 2) {
       expected_a = true_b + 0.05 * exp(0.16 * pow(fmax(0.0, true_b), 1.35)) + random_dins;
     } else if (type == 3) {
       // Original formula: exp(-0.5 * pow((true_b - 1.5)/2.0, 2))
       // R's dnorm(x, mu, sd) = 1/(sqrt(2*pi)*sd) * exp(- (x-mu)^2 / (2*sd^2) )
       // Let mu=1.5, sd=sqrt(2). Then exp(- (x-1.5)^2 / (2*sd^2)) = exp(- (x-1.5)^2 / 4)
       // The formula in the code was exp[-0.125 * (tau_i - 1.50)^2] = exp[-(tau_i-1.5)^2 / 8] -> sd=2
       // Let's use the code's formula:
       expected_a = true_b - exp(-0.125 * (true_b - 1.5) * (true_b - 1.5)) + random_dins;
     } else { // type == 0 (linear + potential non-selectivity)
       expected_a = b0 + true_b * b1 + nonselectivity_offset + random_dins;
     }
     
     // Generate measurement errors
     double error_x = 0.0, error_y = 0.0;
     if (error_dist == "lt") {
       if (sdx > 0) error_x = sdx * R::rt(dfx);
       if (sdy > 0) error_y = sdy * R::rt(dfy);
     } else { // Default to norm
       if (sdx > 0) error_x = R::rnorm(0.0, sdx);
       if (sdy > 0) error_y = R::rnorm(0.0, sdy);
     }
     
     // Calculate simulated measurements
     double current_b = true_b + error_x;
     double current_a = expected_a + error_y;
     
     // Handle Missing Data (MAR / MNAR)
     bool missing_a = false;
     bool missing_b = false;
     
     if (use_mar) {
       if (R::runif(0, 1) < mar_prob) missing_a = true;
       if (R::runif(0, 1) < mar_prob) missing_b = true;
     }
     if (use_mnar) {
       // Check against threshold BEFORE clamping/NA conversion for mnar0
       if (!missing_a && current_a < mnar_threshold) missing_a = true;
       if (!missing_b && current_b < mnar_threshold) missing_b = true;
     }
     
     // Handle Negative Values & Final Assignment
     // If mnar0 logic is active, negatives below threshold become NA handled by 'missing_x' flag
     // Otherwise, clamp negatives at 0.
     if (!missing_a) {
       if (!use_mnar0_logic && current_a < 0.0) {
         current_a = 0.0;
       }
       MP_A[current_idx] = current_a;
       if (!AR) sum_a += current_a; // Add to sum only if not missing
       count_a++;                  // Count non-missing replicates
     } else {
       MP_A[current_idx] = NA_REAL;
     }
     
     if (!missing_b) {
       if (!use_mnar0_logic && current_b < 0.0) {
         current_b = 0.0;
       }
       MP_B[current_idx] = current_b;
       if (!AR) sum_b += current_b; // Add to sum only if not missing
       count_b++;                  // Count non-missing replicates
     } else {
       MP_B[current_idx] = NA_REAL;
     }
     
   } // end replicate loop (r)
   
   // Calculate Mean of Replicates if AR = FALSE
   if (!AR) {
     SampleID_mor[i] = i + 1;
     MP_A_mor[i] = (count_a > 0) ? (sum_a / static_cast<double>(count_a)) : NA_REAL;
     MP_B_mor[i] = (count_b > 0) ? (sum_b / static_cast<double>(count_b)) : NA_REAL;
   }
   
 } // end sample loop (i)
 
 // --- Prepare Output ---
 List sim_data;
 if (AR) {
   sim_data = List::create(Named("SampleID") = SampleID,
                           Named("ReplicateID") = ReplicateID,
                           Named("MP_A") = MP_A,
                           Named("MP_B") = MP_B);
 } else {
   sim_data = List::create(Named("SampleID") = SampleID_mor,
                           Named("MP_A") = MP_A_mor,
                           Named("MP_B") = MP_B_mor);
 }
 
 if (include_parameters) {
   List used_parameters = List::create(
     Named("n") = n,
     Named("R") = R,
     Named("cvx") = cvx,
     Named("cvy") = cvy,
     Named("base_sdx") = base_x, // Base SD used
     Named("base_sdy") = base_y, // Base SD used
     Named("cil") = cil,
     Named("ciu") = ciu,
     Named("dist") = dist,
     Named("df_tau") = (dist == "lst" ? df_tau : NA_REAL),
     Named("eta") = eta,
     Named("eta0") = eta0,
     Named("prop") = prop,
     Named("qpos") = (use_qran_nonselectivity ? qpos : NA_INTEGER),
     Named("qran") = (use_qran_nonselectivity ? qran : NA_REAL),
     Named("qdir") = (use_qran_nonselectivity ? qdir : NA_INTEGER),
     // Named("qlim_lower") = (use_qran_nonselectivity ? qlim_lower : NA_REAL), // Maybe too much detail?
     // Named("qlim_upper") = (use_qran_nonselectivity ? qlim_upper : NA_REAL), // Maybe too much detail?
     Named("mmax") = mmax,
     Named("b0") = b0,
     Named("b1") = b1,
     Named("c0") = c0,
     Named("c1") = c1,
     Named("error_dist") = error_dist,
     Named("dfx") = (error_dist == "lt" ? dfx : NA_REAL),
     Named("dfy") = (error_dist == "lt" ? dfy : NA_REAL),
     Named("md_method") = md_method,
     Named("mar_prob") = (use_mar ? mar_prob : NA_REAL),
     Named("mnar_threshold") = (use_mnar ? mnar_threshold : NA_REAL),
     Named("g_provided") = g_exists,
     Named("type") = (g_exists ? NA_INTEGER : type),
     Named("tau_provided") = obs_tau_exists
   );
   // Add tau values if they were generated/provided and parameter output is requested
   // used_parameters["tau_values"] = tau; // Can make output large
   
   List out = List::create(Named("simulated_data") = sim_data,
                           Named("parameters") = used_parameters);
   return out;
 } else {
   return sim_data;
 }
}

// Ensure Rcpp knows about the custom function helper (if needed elsewhere,
// otherwise being static inline or just defined above main function is fine)
// You might not need this explicit export if it's just a helper.
// NumericVector custom_function(NumericVector x, Rcpp::Function func){
//   NumericVector output = func(x);
//   return output;
// }

