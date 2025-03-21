% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{global_precision_estimates}
\alias{global_precision_estimates}
\title{Calculate Repeatability Variability Components for a Particular IVD-MD Comparison}
\arguments{
\item{data}{A \code{list} or \code{data.table}. Must contain the following
variables:
\itemize{
  \item \code{SampleID: } A \code{character} vector. The sample
                          identifiers.
  \item \code{ReplicateID: } A \code{character} vector. The
                             replicate measurement identifiers.
  \item \code{MP_A: } A \code{numeric} vector. The measurements
                      from IVD-MD \code{MP_A} (response).
  \item \code{MP_B: } A \code{numeric} vector. The measurements
                      from IVD-MD \code{MP_B} (predictor).
}}
}
\value{
A \code{list} of length five. Each element contains the an estimated
statistic. See details for information on each of them.
}
\description{
Calculate various repeatability components, such as variance, coefficient
of variability and the ratio of variances.
}
\details{
Five statistics are estimated. The repeatability variances of \code{MP_A}
and \code{MP_B} are estimated by calculating the \eqn{n} (number of
unique values in \code{SampleID}) sample variances using the \code{numeric}
values in \code{MP_A} and \code{MP_B}. Then we take the mean of the sample
variances to obtain pooled estimated variances for the true repeatability
variances. We denote these estimated pooled variances by \code{Var_A} and
\code{Var_B}, respectively. Using the grand mean of \code{MP_A} and
\code{MP_B} we can then calculate three other statistics, \code{CV_A},
\code{CV_B} and \code{lambda}. Here is a summary the different estimated
statistics:
\itemize{
  \item \code{Var_A: } Pooled variance of all samplewise variances based on
                       IVD-MD \code{MP_A}. An estimator for
                       \eqn{\sigma_v^2}. Denote the estimator
                       \eqn{\hat{\sigma}_v^2}
  \item \code{Var_B: } Pooled variance of all samplewise variances based on
                       IVD-MD \code{MP_B}. An estimator for
                       \eqn{\sigma_h^2}. Denote the estimator
                       \eqn{\hat{\sigma}_h^2}
  \item \code{CV_A: } Estimated repeatbility coefficient of variation.
                      calculated by \eqn{\hat{\sigma}_v / \overline{y}},
                      where \eqn{\overline{y}} is the grand sample mean of
                      the values in \code{MP_A}.
  \item \code{CV_B: } Estimated repeatbility coefficient of variation.
                      calculated by \eqn{\hat{\sigma}_h / \overline{x}},
                      where \eqn{\overline{x}} is the grand sample mean of
                      the values in \code{MP_B}.
  \item \code{lambda: } Estimated repeatbility variance ratio. 
                        calculated by
                        \eqn{\hat{\sigma}_v^2 / \hat{\sigma}_h^2}.
}

By default, \code{CV_A} and \code{CV_B} are represented as a decimal
number. These values can also be represented as percentages, and to
covert to percentages, one should multiply their raw values with \code{100}.

Note: If one uses log-transformed \code{data}, the interpretation of
\code{CV_A} and \code{CV_B} may change, depending on the application.
In the log-transformation case, the square-root of \code{Var_A} and
\code{Var_B} have a similar interpretation as \code{CV_A} and \code{CV_B}
calculated using raw \code{data}.
}
\examples{
library(fasteqa)
# Calculate imprecision estimates
repeatability_statistics <- global_precision_estimates(test_data)

# Output
print(repeatability_statistics)

# Convert CV_A and CV_B to percentages
repeatability_statistics$CV_A <- repeatability_statistics$CV_A * 100
repeatability_statistics$CV_B <- repeatability_statistics$CV_B * 100

# Round results (two decimals)
repeatability_statistics <- lapply(X = repeatability_statistics,
                                   FUN = round,
                                   digits = 2L)

# Convert to data.frame
repeatability_statistics <- as.data.frame(repeatability_statistics)

print(repeatability_statistics)

}
