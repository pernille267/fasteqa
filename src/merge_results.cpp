#include <Rcpp.h>
using namespace Rcpp;

// Helper function to create an index map for quick lookup
std::map<String, int> create_index_map(const CharacterVector& keys) {
  std::map<String, int> index_map;
  for (int i = 0; i < keys.size(); ++i) {
    index_map[keys[i]] = i;
  }
  return index_map;
}

//' Merge all commutability evaluation computations efficiently into one object
//' 
//' @title Merge all computations efficiently into one object
//' @name merge_results
//' @param pb_data A \code{list} or \code{data.table}. This should contain the
//'        prediction band data. Estimated prediction band data is a collection
//'        of estimated pointwise prediction intervals evaluated along a given
//'        grid of predictor values. Must contain:
//'        \itemize{
//'           \item \code{comparison: } A \code{character} vector. The
//'                                     comparison identifiers. Typically on
//'                                     the form MP_A - MP_B.
//'           \item \code{predictor: } A \code{numeric} vector. The \code{MP_B}
//'                                    values used to predict the \code{MP_A}
//'                                    values based on the model.
//'           \item \code{prediction: } A \code{numeric} vector. The predicted
//'                                     \code{MP_A} values based on the
//'                                     predictor values \code{MP_B}.
//'           \item \code{lwr: } A \code{numeric} vector. The lower bound of the
//'                              pointwise prediction intervals.
//'           \item \code{upr: } A \code{numeric} vector. The upper bound of the
//'                              pointwise prediction intervals.
//'        }
//' @param ce_data A \code{list} or \code{data.table}. This should contain the
//'                estimated prediction intervals for evaluated data along with
//'                observed values of \code{MP_B} and \code{MP_B}.
//'                Must contain:
//'                \itemize{
//'                   \item \code{comparison: } A \code{character} vector. The
//'                                             comparison identifiers.
//'                                             Typically on the form
//'                                             \code{'MP_A - MP_B'}.
//'                   \item \code{SampleID: } A \code{character} vector. The
//'                                           sample identifiers for the
//'                                           evaluated samples.
//'                   \item \code{MP_A: } A \code{numeric} vector. The observed
//'                                       measurement results from IVD-MD
//'                                       \code{MP_A} (response).
//'                   \item \code{MP_B: } A \code{numeric} vector. The observed
//'                                       measurement results from IVD-MD
//'                                       \code{MP_B} (predictor).
//'                   \item \code{prediction: } A \code{numeric} vector. The
//'                                             predicted values of \code{MP_A}
//'                                             based on \code{MP_B}.
//'                   \item \code{lwr: } A \code{numeric} vector. The lower
//'                                      bound of the prediction intervals of
//'                                      \code{MP_A}.
//'                   \item \code{upr: } A \code{numeric} vector. The upper
//'                                      bound of the prediction intervals of
//'                                      \code{MP_A}.
//'                   \item \code{inside: } An \code{integer} vector. Should be
//'                                         \code{1} if \code{MP_A}
//'                                         \eqn{\in [\mathrm{lwr},
//'                                         \mathrm{upr}]}, and \code{0}
//'                                         otherwise.
//'                }
//' @param zeta_data A \code{list} or \code{data.table}. This should contain
//'                  relevant statistics of \eqn{\hat{\zeta}}. Must contain:
//'                  \itemize{
//'                     \item \code{comparison: } A \code{character} vector.
//'                                               The comparison identifiers.
//'                                               Typically on the form
//'                                               \code{'MP_A - MP_B'}.
//'                     \item \code{zeta: } A \code{numeric} vector. The
//'                                         point estimates of \eqn{\zeta}.
//'                                         Calculated using
//'                                         \code{estimate_zeta()} or
//'                                         \code{estimate_zeta_ss()}
//'                                         (\code{smooth.commutability}
//'                                         package).
//'                     \item \code{lwr: } A \code{numeric} vector. The lower
//'                                        bound of an estimated confidence
//'                                        interval for \eqn{\hat{\zeta}}.
//'                     \item \code{upr: } A \code{numeric} vector. The upper
//'                                        bound of an estimated confidence
//'                                        interval for \eqn{\hat{\zeta}}.
//'                     \item \code{zeta_critical: } A \code{numeric} vector.
//'                                                  The largest value
//'                                                  \eqn{\hat{\zeta}} can take
//'                                                  to conclude that a
//'                                                  particular IVD-MD
//'                                                  comparison hav acceptable
//'                                                  differences in 
//'                                                  non-selectivity.
//'                     \item \code{zeta_conclusion: } An \code{integer} vector.
//'                                                    Should be \code{1} if
//'                                                    \code{zeta} \eqn{\geq}
//'                                                    \code{zeta_critical} and
//'                                                    \code{0} otherwise.
//'                  }
//' @param imprecision_data A \code{list} or \code{data.table}. Repeatability
//'                         uncertainty estimates and corresponding confidence
//'                         intervals. Must contain:
//'                         \itemize{
//'                             \item \code{comparison: } A \code{character}
//'                             vector. The comparison identifiers. Typically
//'                             on the form \code{'MP_A - MP_B'}.
//'                             \item \code{CV_A: } A \code{numeric} vector.
//'                             The point estimates of the coefficient of
//'                             variability for the repeatability in
//'                             \code{MP_A}.
//'                             \item \code{CV_A_lwr: } A \code{numeric} vector.
//'                             The lower bounds of a confidence interval of
//'                             the coefficient of variability for the
//'                             repeatability in \code{MP_A}.
//'                             \item \code{CV_A_upr: } A \code{numeric} vector.
//'                             The upper bounds of a confidence interval of
//'                             the coefficient of variability for the
//'                             repeatability in \code{MP_A}.
//'                             \item \code{CV_B: } A \code{numeric} vector.
//'                             The point estimates of the coefficient of
//'                             variability for the repeatability in
//'                             \code{MP_B}.
//'                             \item \code{CV_B_lwr: } A \code{numeric} vector.
//'                             The lower bounds of a confidence interval of
//'                             the coefficient of variability for the
//'                             repeatability in \code{MP_B}.
//'                             \item \code{CV_B_upr: } A \code{numeric} vector.
//'                             The upper bounds of a confidence interval of
//'                             the coefficient of variability for the
//'                             repeatability in \code{MP_B}.
//'                             \item \code{lambda: } A \code{numeric} vector.
//'                             The point estimates of the repeatability
//'                             variance ratio between \code{MP_A} and
//'                             \code{MP_B}.
//'                             \item \code{lambda_lwr: } A \code{numeric}
//'                             vector. The lower bounds of a confidence
//'                             interval of \code{lambda}.
//'                             \item \code{lambda_upr: } A \code{numeric}
//'                             vector. The upper bounds of a confidence
//'                             interval of \code{lambda}.
//'                         }
//'                         This can be omitted.
//' @param rounding An \code{integer}. The maximum number of decimals to use
//'                 for the output \code{numeric} values. 
//' @param include_imprecision_estimates An \code{logical} value. If
//'                                      \code{TRUE}, imprecision estimate are
//'                                      included in the output. In this case,
//'                                      \code{imprecision_data} cannot be
//'                                      omitted.
//' @param silence An \code{integer}. Set to \eqn{\leq} \code{0} to allow
//'                verbose.
//' 
//' @description
//' Merge all the essential components for a complete commutability evaluation
//' data analysis. This merging allows for presenting the commutability
//' evaluation results in both plots and tables. 
//' 
//' @details
//' The merging is done with respect to the \code{comparison} column. Hence,
//' \code{comparison} must be part of all input components (i.e.,
//' \code{pb_data}, \code{ce_data}, \code{zeta_data}, \code{imprecision_data}).
//' 
//' In practice, this function results in two merged datasets.
//' 
//' (1) \code{ pb_data} is merged with \code{zeta_data} and
//' \code{imprecision_data}, using \code{comparison} as the key for the merging
//' process. The output will be a new \code{pb_data} that contain relevant
//' information from \code{zeta_data} and \code{imprecision_data}.
//' 
//' (2) \code{ ce_data} is merged with \code{zeta_data} and
//' \code{imprecision_data}, using \code{comparison} as key for the merging
//' process. The output will be a new \code{ce_data} that contain relevant
//' information from \code{zeta_data} and \code{imprecision_data}.
//' 
//' @return
//' A \code{list} of length two:
//' \itemize{
//'   \item \code{pb_data: } The input \code{pb_data} with additional elements
//'                          added in the merging process with \code{zeta_data}
//'                          and \code{imprecision_data}.
//'   \item \code{ce_data: } The input \code{ce_data} with additional elements
//'                          added in the merging process with \code{zeta_data}
//'                          and \code{imprecision_data}
//' } 
//'
//' @examples \dontrun{
//'   library(fasteqa)
//'   data <- simulate_data_eqa(list(n = 25, R = 3, cvx = 0.06, cvy = 0.04))
//'   estimate_zeta(data)
//' }
//'

 // [[Rcpp::export]]
 List merge_results(List pb_data,
                    List ce_data,
                    List zeta_data,
                    Nullable<List> imprecision_data = R_NilValue,
                    int rounding = 3,
                    bool include_imprecision_estimates = true,
                    int silence = 1) {
   
   CharacterVector comparison_pb = pb_data["comparison"];
   CharacterVector comparison_ce = ce_data["comparison"];
   CharacterVector comparison_zeta = zeta_data["comparison"];
   
   int N_pb = comparison_pb.size();
   int N_ce = comparison_ce.size();
   
   // Create maps for quick lookup
   std::map<String, int> zeta_map = create_index_map(comparison_zeta);
   
   // Extract zeta data
   NumericVector zeta_pe = zeta_data["zeta"];
   NumericVector zeta_lwr = zeta_data["lwr"];
   NumericVector zeta_upr = zeta_data["upr"];
   NumericVector zeta_crit = zeta_data["zeta_critical"];
   IntegerVector zeta_conc = zeta_data["zeta_conclusion"];
   
   // Prepare output vectors for pb_data
   NumericVector pb_zeta(N_pb), pb_zeta_lwr(N_pb), pb_zeta_upr(N_pb), pb_zeta_crit(N_pb);
   IntegerVector pb_zeta_conc(N_pb);
   
   // Fill pb_data with corresponding zeta values
   for (int i = 0; i < N_pb; ++i) {
     int idx = zeta_map[comparison_pb[i]];
     pb_zeta[i]      = zeta_pe[idx];
     pb_zeta_lwr[i]  = zeta_lwr[idx];
     pb_zeta_upr[i]  = zeta_upr[idx];
     pb_zeta_crit[i] = zeta_crit[idx];
     pb_zeta_conc[i] = zeta_conc[idx];
   }
   
   // Prepare output vectors for ce_data
   NumericVector ce_zeta(N_ce), ce_zeta_lwr(N_ce), ce_zeta_upr(N_ce), ce_zeta_crit(N_ce);
   IntegerVector ce_zeta_conc(N_ce);
   
   // Fill ce_data with corresponding zeta values
   for (int i = 0; i < N_ce; ++i) {
     int idx = zeta_map[comparison_ce[i]];
     ce_zeta[i]      = zeta_pe[idx];
     ce_zeta_lwr[i]  = zeta_lwr[idx];
     ce_zeta_upr[i]  = zeta_upr[idx];
     ce_zeta_crit[i] = zeta_crit[idx];
     ce_zeta_conc[i] = zeta_conc[idx];
   }
   
   // Include imprecision data if requested
   bool has_imprecision = include_imprecision_estimates && imprecision_data.isNotNull();
   
   NumericVector CV_A, CV_A_lwr, CV_A_upr, CV_B, CV_B_lwr, CV_B_upr, lambda, lambda_lwr, lambda_upr;
   
   std::map<String, int> imprec_map;
   
   if (has_imprecision) {
     List imp_data(imprecision_data);
     CharacterVector comparison_imp = imp_data["comparison"];
     imprec_map = create_index_map(comparison_imp);
     
     CV_A       = imp_data["CV_A"];
     CV_A_lwr   = imp_data["CV_A_lwr"];
     CV_A_upr   = imp_data["CV_A_upr"];
     CV_B       = imp_data["CV_B"];
     CV_B_lwr   = imp_data["CV_B_lwr"];
     CV_B_upr   = imp_data["CV_B_upr"];
     lambda     = imp_data["lambda"];
     lambda_lwr = imp_data["lambda_lwr"];
     lambda_upr = imp_data["lambda_upr"];
     
     // Expand imprecision data to match pb and ce lengths
  #define EXPAND_IMPRECISION(DATA_VEC, COMP_VEC) \
     ({ NumericVector res(COMP_VEC.size());    \
  for(int k=0;k<COMP_VEC.size();++k) res[k]=DATA_VEC[imprec_map[COMP_VEC[k]]]; res; })

  List merged_pb_data =
    List::create(
      _["comparison"]       = comparison_pb,
      _["zeta"]             = round(pb_zeta, rounding),
      _["zeta_ci_lwr"]      = round(pb_zeta_lwr, rounding),
      _["zeta_ci_upr"]      = round(pb_zeta_upr, rounding),
      _["zeta_upper"]       = round(pb_zeta_crit, rounding),
      _["dins_conclusion"]  = pb_zeta_conc,
      _["predictor"]        = round(as<NumericVector>(pb_data["predictor"]), rounding),
      _["prediction"]       = round(as<NumericVector>(pb_data["prediction"]), rounding),
      _["pi_lwr"]           = round(as<NumericVector>(pb_data["lwr"]), rounding),
      _["pi_upr"]           = round(as<NumericVector>(pb_data["upr"]), rounding),
      _["CV_A"]             = round(EXPAND_IMPRECISION(CV_A, comparison_pb), rounding),
      _["CV_A_lwr"]         = round(EXPAND_IMPRECISION(CV_A_lwr, comparison_pb), rounding),
      _["CV_A_upr"]         = round(EXPAND_IMPRECISION(CV_A_upr, comparison_pb), rounding),
      _["CV_B"]             = round(EXPAND_IMPRECISION(CV_B, comparison_pb), rounding),
      _["CV_B_lwr"]         = round(EXPAND_IMPRECISION(CV_B_lwr, comparison_pb), rounding),
      _["CV_B_upr"]         = round(EXPAND_IMPRECISION(CV_B_upr, comparison_pb), rounding),
      _["lambda"]           = round(EXPAND_IMPRECISION(lambda, comparison_pb), rounding),
      _["lambda_lwr"]       = round(EXPAND_IMPRECISION(lambda_lwr, comparison_pb), rounding),
      _["lambda_upr"]       = round(EXPAND_IMPRECISION(lambda_upr, comparison_pb), rounding)
    );
  
  List merged_ce_data =
    List::create(
      _["comparison"]       = comparison_ce,
      _["SampleID"]         = ce_data["SampleID"],
      _["zeta"]             = round(ce_zeta, rounding),
      _["zeta_ci_lwr"]      = round(ce_zeta_lwr, rounding),
      _["zeta_ci_upr"]      = round(ce_zeta_upr, rounding),
      _["zeta_upper"]       = round(ce_zeta_crit, rounding),
      _["dins_conclusion"]  = ce_zeta_conc,
      _["MP_B"]             = round(as<NumericVector>(ce_data["MP_B"]), rounding),
      _["MP_A"]             = round(as<NumericVector>(ce_data["MP_A"]), rounding),
      _["prediction"]       = round(as<NumericVector>(ce_data["prediction"]), rounding),
      _["pi_lwr"]           = round(as<NumericVector>(ce_data["lwr"]), rounding),
      _["pi_upr"]           = round(as<NumericVector>(ce_data["upr"]), rounding),
      _["pi_inside"]        = ce_data["inside"],
      _["CV_A"]             = round(EXPAND_IMPRECISION(CV_A, comparison_ce), rounding),
      _["CV_A_lwr"]         = round(EXPAND_IMPRECISION(CV_A_lwr, comparison_ce), rounding),
      _["CV_A_upr"]         = round(EXPAND_IMPRECISION(CV_A_upr, comparison_ce), rounding),
      _["CV_B"]             = round(EXPAND_IMPRECISION(CV_B, comparison_ce), rounding),
      _["CV_B_lwr"]         = round(EXPAND_IMPRECISION(CV_B_lwr, comparison_ce), rounding),
      _["CV_B_upr"]         = round(EXPAND_IMPRECISION(CV_B_upr, comparison_ce), rounding),
      _["lambda"]           = round(EXPAND_IMPRECISION(lambda, comparison_ce), rounding),
      _["lambda_lwr"]       = round(EXPAND_IMPRECISION(lambda_lwr, comparison_ce), rounding),
      _["lambda_upr"]       = round(EXPAND_IMPRECISION(lambda_upr, comparison_ce), rounding)
    );
  
 #undef EXPAND_IMPRECISION
  
  return List::create(_["merged_pb_data"]=merged_pb_data,_["merged_ce_data"]=merged_ce_data);
  
  }
  else {
    List merged_pb_data = List::create(
      _["comparison"]      = comparison_pb,
      _["zeta"]            = round(pb_zeta, rounding),
      _["zeta_ci_lwr"]     = round(pb_zeta_lwr, rounding),
      _["zeta_ci_upr"]     = round(pb_zeta_upr, rounding),
      _["zeta_upper"]      = round(pb_zeta_crit, rounding),
      _["dins_conclusion"] = pb_zeta_conc,
      _["predictor"]       = round(as<NumericVector>(pb_data["predictor"]), rounding),
      _["prediction"]      = round(as<NumericVector>(pb_data["prediction"]), rounding),
      _["pi_lwr"]          = round(as<NumericVector>(pb_data["lwr"]), rounding),
      _["pi_upr"]          = round(as<NumericVector>(pb_data["upr"]), rounding)
    );
    
    List merged_ce_data = List::create(
      _["comparison"]      = comparison_ce,
      _["SampleID"]        = ce_data["SampleID"],
      _["zeta"]            = round(ce_zeta, rounding),
      _["zeta_ci_lwr"]     = round(ce_zeta_lwr, rounding),
      _["zeta_ci_upr"]     = round(ce_zeta_upr, rounding),
      _["zeta_upper"]      = round(ce_zeta_crit, rounding),
      _["dins_conclusion"] = ce_zeta_conc,
      _["MP_B"]            = round(as<NumericVector>(ce_data["MP_B"]), rounding),
      _["MP_A"]            = round(as<NumericVector>(ce_data["MP_A"]), rounding),
      _["prediction"]      = round(as<NumericVector>(ce_data["prediction"]), rounding),
      _["pi_lwr"]          = round(as<NumericVector>(ce_data["lwr"]), rounding),
      _["pi_upr"]          = round(as<NumericVector>(ce_data["upr"]), rounding),
      _["pi_inside"]       = ce_data["inside"]
    );
    
    return List::create(
      _["merged_pb_data"] = merged_pb_data,
      _["merged_ce_data"] = merged_ce_data
    );
  }
   
 }
