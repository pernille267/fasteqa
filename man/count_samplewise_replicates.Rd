% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{count_samplewise_replicates}
\alias{count_samplewise_replicates}
\title{Count the Number of Replicated for Each Sample}
\arguments{
\item{data}{A \code{list} or \code{data.table} representing CS data or EQAM data
with the \code{list} elements or \code{data.table} columns:
\code{SampleID}, \code{MP_A}, and \code{MP_B}.}

\item{summary}{A \code{character} specifying the summary statistic of the sample-wise
numbers of replicates. Default is \code{'mode'}. Possible summary statistics include:
\itemize{
   \item \code{none}: Returns an \code{integer} vector containing every
   sample-wise number of replicatees for \code{data}.
   \item \code{mode}: Returns the sample mode.
   \item \code{median}: Returns the sample median.
   \item \code{mean}: Returns the sample mean.
   \item \code{ceiling}: Returns the rounded up (to nearest integer) sample mean.
   \item \code{floor}: Returns the rounded down (to nearest integer) sample mean.
   \item \code{round}: Returns the rounded (to nearest integer) sample mean.
}}

\item{invalid_NA}{A \code{logical} that determines the behavior of the function in response
to invalid \code{data} or \code{summary} input.
If \code{TRUE}, the function will return an \code{NA}-value when encountering
invalid input or computation errors, rather than throwing an error.
While generally not recommended due to potential masking of issues,
this may be useful in certain scenarios where uninterrupted execution is desired.}

\item{silence}{An \code{integer} that dictates the verbosity level of the console output.
If \code{silence} is set to a value less than 1, various verbose output,
including debugging reports, will be displayed depending on the specific
value of \code{silence}.}
}
\value{
Returns a \code{list} that contains a single element, \code{R_i}.
The type of \code{R_i} depends on the \code{summary} parameter.
It is typically an single integer value.
However, if \code{summary} is set to \code{'none'}, \code{R_i} is an
integer vector because it list the number of replicates for each sample.
When \code{summary} is set to \code{'mean'},
\code{R_i} is returned as a \code{double}.
}
\description{
This function counts the number of replicated measurements done on each sample
within an IVD-MD comparison for either clinical sample (CS) data or
external quality asessment material (EQAM) data. Alternatively, a summary
statistic can be returned to summarize the numbers of replicates across all
samples in the \code{data}.
}
\details{
The function \code{predict_eqa()} hinges on the count of replicates
performed for its \code{data} parameter. When the \code{method} is either
\code{fg} or \code{ols}, it is essential to tally the replicates for both
\code{data} and \code{new_data} inputs. This handy function offers a
streamlined solution. It allows for the direct usage of \code{predict_eqa()}
without requiring manual counting of the sample-wise replicates,
thereby enhancing efficiency and ease of use.

Note: If the count of replicates for a particular sample is 0, it will be
removed from the output.
}
\examples{
\dontrun{
library(fasteqa)
# Simulation parameters for clinical sample data
cs_parameters <- list(n = 25, R = 3,
                      cvx = 0.01, cvy = 0.015,
                      cil = 10, ciu = 70)
# Use the simulation parameters to simulate toy clinical sample data                      
cs_data <- simulate_eqa_data(cs_parameters)
# Calculate mode of sample-wise number of replicates
mode_R <- count_samplewise_replicates(cs_data, summary = 'mode')$R_i
# Count sample-wise number of replicates
samplewise_R <- count_samplewise_replicates(cs_data, summary = 'none')$R_i                                    
}
}
