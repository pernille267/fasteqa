% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{predict_eqa}
\alias{predict_eqa}
\title{Prediction Interval Estimation for EQA Data via Deming or
       OLS Regression}
\arguments{
\item{data}{A \code{list} or \code{data.table}. Must contain:
\itemize{
    \item \code{SampleID: } A \code{character} vector. The
    clinical sample identifiers.
    \item \code{MP_A: } A \code{numeric} vector. The means of
    replicated measurements from IVD-MD \code{MP_A} (response).
    \item \code{MP_B: } A \code{numeric} vector. The means of
    replicated measurements from IVD-MD \code{MP_B} (predictor).
}}

\item{new_data}{A \code{list} or \code{data.table}. Can contain:
\itemize{
    \item \code{SampleID: } A \code{character} vector. The
    evaluated material sample identifiers. Optional.
    \item \code{MP_A: } A \code{numeric} vector. The means
    of replicated measurements from IVD-MD \code{MP_A}.
    Optional.
    \item \code{MP_B: } A \code{numeric} vector. The means
    of replicated measurements from IVD-MD \code{MP_B}.
    Mandatory.
}
The structure of \code{new_data} determines the output
type. For constructing prediction band (PB) data, include
only \code{MP_B}. For constructing commutability
evaluation data, include all three. Inside checks can only
be performed if \code{MP_A} is included.}

\item{imprecision_estimates}{A \code{list} or \code{data.table}. Must
include:
\itemize{
   \item \code{Var_B: } Pooled variance 
         estimate of \eqn{\sigma_h^2}.
   \item \code{lambda: } Estimate of
         \eqn{\sigma_v^2 / \sigma_h^2}. 
}
See details.}

\item{R}{An \code{integer}. The average number of replicates on which
\code{new_data} is based. \code{count_samplewise_replicates()} can
be employed to generate a suitable input for this parameter
directly.}

\item{R_ratio}{A \code{double}. The ratio of the number of replicates used
in \code{new_data} and \code{data}. Only relevant if
\code{method = 'fg'} or \code{method = 'ols'}. Defaults to
\code{1L}.}

\item{method}{A \code{character} string. The desired method for estimating
prediction intervals. Possible prediction estimation
approaches includes:
\itemize{
    \item \code{fg: } Implements standard Deming regression
          with prediction intervals calculated using the
          'Fuller & Gillard' approach. See details.
    \item \code{clsi: } Utilizes standard Deming regression,
          with prediction intervals derived from the CLSI
          method. See details.
    \item \code{ols: } Implements Ordinary Least Squares
          regression.
}}

\item{level}{A \code{double}. Must be between \code{0} and \code{1}. The
nominal confidence level for the estimated prediction
intervals. Defaults to \code{0.99} (\eqn{99\%}).}

\item{allow_reverse_regression}{A \code{logical} value. If \code{TRUE} and
\code{method = 'ols'}, the response and
predictor change roles if
\code{lambda < 0.5}. Defaults to
\code{FALSE}.}

\item{rounding}{An \code{integer}. The desired number of decimal places for
the predictions and prediction intervals. Defaults to
\code{3L}, offering sufficient precision in relevant
applications.}
}
\value{
A \code{list}. The resulting estimated prediction interval data based on
the function inputs. Contain the following elements:
\itemize{
    \item \code{SampleID: } A \code{character} vector. The ID(s) of the
          evaluated material(s).
    \item \code{MP_B: } A \code{numeric} vector. The measurement result(s)
          of IVD-MD \code{MP_B} (predictor).
    \item \code{MP_A: } A \code{numeric} vector. The measurement result(s)
          of IVD-MD \code{MP_A} (response).
    \item \code{prediction: } A \code{numeric} vector. The predicted
          value(s) of \code{MP_A} given the value(s) of \code{MP_B}.
    \item \code{lwr: } A \code{numeric} vector. The lower limit(s) of the
          estimated prediction interval(s) of \code{MP_A} given the
          value(s) of \code{MP_B}.
    \item \code{upr: } A \code{numeric} vector. The upper limit(s) of the
           estimated prediction interval(s) of \code{MP_A} given the
           value(s) of \code{MP_B}.
    \item \code{inside: } An \code{integer} vector. The Inside checks.
          \code{1} if \code{MP_A} is inside the estimated prediction
          interval for \code{MP_A}. Otherwise, \code{0}.
}
Note: the output \code{list} may be converted to \code{data.table} by using
the \code{setDT()} function from the \code{data.table} package in R.
}
\description{
This function estimates prediction intervals for evaluated material data by
applying either Deming (\code{method = 'fg'} or \code{method = 'clsi'}) or
Ordinary Least Squares (\code{method = 'ols'}) regression methodologies.
}
\details{
For commutability evaluation purposes, inclusion of all \code{SampleID},
\code{MP_A}, and \code{MP_B} in \code{new_data} is required. 
If only prediction band data is required, inclusion of \code{MP_B} in
\code{new_data} is the sole requirement.

Imprecision Estimates

Estimating prediction intervals may require certain imprecision estimates,
which can be estimated from the raw data. Which imprecision estimates that
are required depends on which method that is used to construct the
prediction intervals. For \code{method = 'fg'} or \code{method = 'ols'},
\code{lambda} is required. \code{lambda} is calculated using

\eqn{\lambda = \frac{\hat{\sigma}_v^2}{\hat{\sigma}_h^2}},

where \eqn{\hat{\sigma}_v^2} and \eqn{\hat{\sigma}_h^2} are pooled 
variance estimates. If differences in non-selectivity (DINS) or equation
error cannot be ruled out, \eqn{\lambda} will generally underestimate the
true parameter \eqn{\Lambda}. Keep this in mind!

For \code{method = 'clsi'}, both \code{Var_B} and \code{lambda} are
required. \code{lambda} is calculated as before, and \code{Var_B} is just
\eqn{\hat{\sigma}_h^2}.

Note that both \code{lambda} and \code{Var_B} always must be given as
input. However, if one of them is not necessary to estimate prediction
intervals, it may take an arbitrary value. It is nevertheless recommended
to always include actual values \code{lambda} and \code{Var_B},
independently of whether they are actual used or not.

Note also that \code{global_precision_estimates()} can be used to calculate
both \code{lambda} and \code{Var_B}.
}
\examples{
# Required packages
library(fasteqa)
library(data.table)

# Read data and convert to data.table
test_data_example <- as.data.table(test_data)

# Log-transform data
test_data_example[, MP_A := log(MP_A)]
test_data_example[, MP_B := log(MP_B)]

# Use one of the clinical samples as a fictive evaluated material sample
test_cs_data <- test_data_example[SampleID != "1"]
test_eq_data <- test_data_example[SampleID == "1"]

# Estimate repeatability uncertainty statistics
impr_data <- global_precision_estimates(test_cs_data)

# Calculate mean-of-replicates data
test_cs_data <- test_cs_data[, fun_of_replicates(.SD)]
test_eq_data <- test_eq_data[, fun_of_replicates(.SD)]

# Calculate 95\% OLS commutability evaluation data
ols_pi <- predict_eqa(data = test_cs_data,
                      new_data = test_eq_data,
                      imprecision_estimates = impr_data,
                      method = "ols",
                      level = 0.95,
                      allow_reverse_regression = FALSE,
                      rounding = 3L)

# Calculate 95\% F-G Deming commutability evaluation data
fg_pi <- predict_eqa(data = test_cs_data,
                     new_data = test_eq_data,
                     imprecision_estimates = impr_data,
                     method = "fg",
                     level = 0.95,
                     allow_reverse_regression = FALSE,
                     rounding = 3L)

# Calculate 95\% CLSI EP14 Deming commutability evaluation data
clsi_pi <- predict_eqa(data = test_cs_data,
                      new_data = test_eq_data,
                      imprecision_estimates = impr_data,
                      method = "clsi",
                      level = 0.95,
                      allow_reverse_regression = FALSE,
                      rounding = 3L)

# Convert to data.table objects
lapply(X = list(ols_pi,
                fg_pi,
                clsi_pi),
       FUN = setDT)

# Gather into one data.table
pis <- rbindlist(list("ols" = ols_pi,
                      "fg" = fg_pi,
                      "clsi" = clsi_pi),
                  idcol = "method")

# The result                  
print(pis)

}
