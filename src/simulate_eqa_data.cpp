#include <Rcpp.h>
using namespace Rcpp;

//' Simulation of EQA data based on study design and potential differences in non-selectivity
//'
//' @title Simulation of EQA data based on study design and potential differences in non-selectivity
//' @name simulate_eqa_data
//' 
//' @param parameters A \code{list} of parameters used to simulate the EQA data. You must at least specify one parameter for this function to run. Except that one mandatory parameter, you may optionally choose the remaining of the parameters. These are the optimal parameters that you may include into the list:
//' \itemize{
//'   \item{\code{n: }}{  The number of clinical samples}
//'   \item{\code{R: }}{The number of replicates on each sample}
//'   \item{\code{cvx: }}{The analytical CV of x measurements}
//'   \item{\code{cvy: }}{The analytical CV of y measurements}
//'   \item{\code{cil: }}{The lower range of the concentration interval}
//'   \item{\code{ciu: }}{The upper range of the concentration interval}
//'   \item{\code{dist: }}{The distribution to simulate latent variables from. Possbile choices include 'unif' (uniform distribution, default), 'norm' (normal distribution) , 'lst' (location-scale t-distribution), 'lnorm' (log-normal distribution)}
//'   \item{\code{df_tau: }}{The degrees of freedom for the 'lst' distribution if the distribution of latent variables are location-scale t-distributed ('lst'). Defaults to 5 if not otherwise is specified.}
//'   \item{\code{eta: }}{The heteroscedasticity factor}
//'   \item{\code{eta0: }}{The proportion of base MS standard deviations}
//'   \item{\code{qpos: }}{Position of systematic differences in non-selectivity. 0 signify lower range and 1 signify upper range}
//'   \item{\code{qran: }}{Interquantile range where systematic differences in non-selectivity should have its effect}
//'   \item{\code{prop: }}{average proportion of clinical samples affected by random differences in non-selectivity}
//'   \item{\code{mmax: }}{The maximum relocation magnitude in number of analytical SDs of y measurements. This assumes either prop or qpos and qran to be specified as well}
//'   \item{\code{b0: }}{For systematic linear DINS between IVD-MDs. Intercept. Defaults to 0.}
//'   \item{\code{b1: }}{For systematic linear DINS between IVD-MDs. Slope. Defaults to 1.}
//'   \item{\code{c0: }}{For systematic linear non-selectivity in IVD-MD 1. Intercept. Defaults to 0.}
//'   \item{\code{c1: }}{For systematic linear non-selectivity in IVD-MD 1. Slope. Defaults to 1.}
//'   \item{\code{error_dist: }}{The distribution to simulate measurement error components from. Possible choices include 'norm' (normal distribution, default) and 'lt' (location t-distribution)}
//'   \item{\code{dfx: }}{The degrees of freedom for the measurement error components in IVD-MD 1 if error_dist = 'lt'. Defaults to 5 if not otherwise is specified.}
//'   \item{\code{dfy: }}{The degrees of freedom for the measurement error components in IVD-MD 2 if error_dist = 'lt'. Defaults to 5 if not otherwise is specified.}
//'   \item{\code{md_method: }}{Method for simulation missing data. Possible choices include 'none' (no missing data is simulated, default), 'mar' (missing at random), 'mnar' (missing not at random) and 'marmnar' (missing both at random and not at random)}
//'   \item{\code{mar_prob: }}{The probability (value between 0 and 1) of having a missing measurement. Only relevant if \code{md_method} is 'mar' or 'marmnar'. If not specified, but \code{md_method} = 'mar' or \code{md_method} = 'marmnar', it defaults to 0.05.}
//'   \item{\code{mnar_threshold: }}{The lower bound threshold (a real value) for when a measurement should be missing. Only relevant if \code{md_method} is 'mnar' or 'marmnar'. If not specified, but \code{md_method} = 'mnar' or \code{md_method} = 'marmnar', it defaults to \code{cil}. Alternatively, if not specified, but \code{md_method} = 'mnar0' or \code{md_method} = 'marmnar0', it defaults to 0.}
//'   
//' }
//' @param silence \code{Integer} should temporary calculation results be printed? This may be useful for debugging or strange curiosity. \code{silence = 1} signify that printing will be suppressed, which is the default. \code{silence = 0} allows such printing 
//' 
//' @description Simulates a data set with n x R rows, and four columns. The two first columns are the base ID columns (\code{SampleID} and \code{ReplicateID}). The remaining columns are numeric columns holding measurement results from the two IVD-MDs in comparison (denoted 'MP_A' (y) and 'MP_B' (x)).
//'              The form of the simulated data depends greatly on which parameters are specified in the the \code{parameters} argument. 
//' @details Simulates method comparison data for \code{n} samples (e.g., clinical samples, pools, external quality assessment samples, reference material samples), where each sample is measured \code{R} times (replicated measurements). In theory, we simulate from (x_ir, y_ir) where x_ir = f(tau_i) + h_ir and y_ir = g(f(tau_i)) + v_ir.
//'          The form of f is specified through parameters \code{c0} and \code{c1}, whereas g is specified through numerous parameters such as \code{b0}, \code{b1}, \code{qpos}, \code{qran}, \code{mmax}, \code{prop}. tau_i is modelled through \code{cil}, \code{ciu}, \code{dist} and \code{df_tau}.
//'          h_ir and v_ir are measurement error components modelled through \code{cvx}, \code{cvy}. \code{cvx}, \code{cvy} can be functions of \code{error_dist}, \code{dfx}, \code{dfy}, \code{eta} and \code{eta0}. 
//'          In order to convert the outputted list to a table, use either \code{as.data.frame()}, \code{as.data.table()}, \code{as.tibble()}. The most efficient way to convert is \code{setDT()} from the \code{data.table} package.
//'
//' @return A list where each list element is a column of the generated EQA data
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
      unsorted_tau[i] = mu_tau + sigma_tau + R::rt(df_tau);
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

