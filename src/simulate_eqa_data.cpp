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
//'   \item{\code{eta: }}{The heteroscedasticity factor}
//'   \item{\code{eta0: }}{The proportion of base MS standard deviations}
//'   \item{\code{qpos: }}{Position of systematic differences in non-selectivity. 0 signify lower range and 1 signify upper range}
//'   \item{\code{qran: }}{Interquantile range where systematic differences in non-selectivity should have its effect}
//'   \item{\code{prop: }}{average proportion of clinical samples affected by random differences in non-selectivity}
//'   \item{\code{mmax: }}{The maximum relocation magnitude in number of analytical SDs of y measurements. This assumes either prop or qpos and qran to be specified as well}
//'   \item{\code{b0: }}{For systematic linear differences in non-selectivity. Intercept}
//'   \item{\code{b1: }}{For systematic linear differences in non-selectivity. Slope}
//' }
//' @param silence \code{Integer} should temporary calculation results be printed? This may be useful for debugging or strange curiosity. \code{silence = 1} signify that printing will be suppressed, which is the default. \code{silence = 0} allows such printing 
//' 
//' @description Simulates a data set with n x R rows, and four columns. The two first columns are the base ID columns (\code{SampleID} and \code{ReplicateID}). The remaining columns are numeric columns holding measurement results from the two MSs in comparison (MS_A and MS_B) denoted \code{MP_A} and \code{MP_B}.
//' @details In order to convert the outputted list to a table, use either \code{as.data.frame()}, \code{as.data.table()}, \code{as.tibble()}. The most efficient way to convert is \code{setDT()} from the \code{data.table} package.
//'
//' @return A list where each list element is a column of the generated EQA data
//'
//' @examples \dontrun{
//'   simulated_data <- simulate_eqa_data(parameters = list(n = 25, R = 3, prop = 0.1, mmax = 5))
//'   simulated_data <- setDT(simulated_data)
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
  int eta_exists = 0;
  int eta0_exists = 0;
  int qpos_exists = 0;
  int qran_exists = 0;
  int prop_exists = 0;
  int mmax_exists = 0;
  int b0_exists = 0;
  int b1_exists = 0;
  
  // Checks which of the parameters found in parameters argument
  CharacterVector given_parameters = parameters.names();
  // Rcout << "given parameters : " << given_parameters << "\n";
  int number_given_parameters = given_parameters.size();
  // Rcout << "number_given_parameters : " << number_given_parameters << "\n";
  for(int i = 0; i < number_given_parameters; ++i){
    
    if(silence == 0){
      Rcout << "ith given parameter : " << given_parameters[i] << "\n";  
    }
    
    CharacterVector candidate_param(14);
    candidate_param[0] = "n";
    candidate_param[1] = "R";
    candidate_param[2] = "cvx";
    candidate_param[3] = "cvy";
    candidate_param[4] = "cil";
    candidate_param[5] = "ciu";
    candidate_param[6] = "eta";
    candidate_param[7] = "eta0";
    candidate_param[8] = "qpos";
    candidate_param[9] = "qran";
    candidate_param[10] = "prop";
    candidate_param[11] = "mmax";
    candidate_param[12] = "b0";
    candidate_param[13] = "b1";
    
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
      eta_exists = 1;
    }
    else if(candidate_param[7] == given_parameters[i]){
      eta0_exists = 1;
    }
    else if(candidate_param[8] == given_parameters[i]){
      qpos_exists = 1;
    }
    else if(candidate_param[9] == given_parameters[i]){
      qran_exists = 1;
    }
    else if(candidate_param[10] == given_parameters[i]){
      prop_exists = 1;
    }
    else if(candidate_param[11] == given_parameters[i]){
      mmax_exists = 1;
    }
    else if(candidate_param[12] == given_parameters[i]){
      b0_exists = 1;
    }
    else if(candidate_param[13] == given_parameters[i]){
      b1_exists = 1;
    }
  }
  if(silence == 0){
    Rcout << "Does n exists? 1 for yes and 0 for no : " << n_exists << "\n";
    Rcout << "Does R exists? 1 for yes and 0 for no : " << R_exists << "\n";
    Rcout << "Does cvx exists? 1 for yes and 0 for no : " << cvx_exists << "\n";
    Rcout << "Does cvy exists? 1 for yes and 0 for no : " << cvy_exists << "\n";
    Rcout << "Does cil exists? 1 for yes and 0 for no : " << cil_exists << "\n";
    Rcout << "Does ciu exists? 1 for yes and 0 for no : " << ciu_exists << "\n";
    Rcout << "Does eta exists? 1 for yes and 0 for no : " << eta_exists << "\n";
    Rcout << "Does eta0 exists? 1 for yes and 0 for no : " << eta0_exists << "\n";
    Rcout << "Does qpos exists? 1 for yes and 0 for no : " << qpos_exists << "\n";
    Rcout << "Does qran exists? 1 for yes and 0 for no : " << qran_exists << "\n";
    Rcout << "Does prop exists? 1 for yes and 0 for no : " << prop_exists << "\n";
    Rcout << "Does mmax exists? 1 for yes and 0 for no : " << mmax_exists << "\n";
    Rcout << "Does b0 exists? 1 for yes and 0 for no : " << b0_exists << "\n";
    Rcout << "Does b1 exists? 1 for yes and 0 for no : " << b1_exists << "\n";
  }

  int n = 0;
  int R = 0;
  float cvx = 0;
  float cvy = 0;
  float cil = 0;
  float ciu = 0;
  float eta = 0;
  float eta0 = 0;
  int qpos = 0;
  float qran = 0;
  float prop = 0;
  float mmax = 0;
  float b0 = 0;
  float b1 = 1;
  
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
    cvx += R::rbeta(2, 5) / 10;
  }
  if(cvy_exists == 1){
    float reg_cvy = parameters["cvy"];
    cvy += reg_cvy;
  }
  else{
    cvy += R::rbeta(2, 5) / 10;
  }
  if(cil_exists == 1){
    float reg_cil = parameters["cil"];
    cil += reg_cil;
  }
  else{
    cil += R::rf(1.057057, 8.15) * 44;
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

  NumericVector unsorted_tau(n);
  NumericVector tau(n);
  int nR = n * R;
  IntegerVector SampleID(nR);
  IntegerVector ReplicateID(nR);
  NumericVector MP_A(nR);
  NumericVector MP_B(nR);
  
  float average = (cil + ciu) / 2;
  
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
  
  for(int i = 0; i < n; ++i){
    unsorted_tau[i] = R::runif(cil, ciu);  
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
      //Rcout << "sdy for i = " << i << " is " << sdy << "\n";  
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
    for(int r = 0; r < R; ++r){
      SampleID[idx[r]] = i + 1;
      ReplicateID[idx[r]] = r + 1;
      MP_A[idx[r]] = tau[i] + R::rnorm(0, sdy) + relocate_sample_i * relocating_magnitude * sdy;
      MP_B[idx[r]] = b0 + tau[i] * b1 + R::rnorm(0, sdx);
    }
  }
  
  List out = List::create(Named("SampleID") = SampleID, Named("ReplicateID") = ReplicateID, Named("MP_A") = round(MP_A, 3), Named("MP_B") = round(MP_B, 3));
  return out;
}

