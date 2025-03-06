\name{fasteqa-package}
\alias{fasteqa}
\docType{package}
\title{
  Utility Functions for External Quality Assessment Material Data Analysis
}
\description{
  Streamline commutability evaluation data analysis by providing functions
  that effectively calculates relevant statistics. It also provides functions
  that allow standardization of output structure. This package can be viewed
  as backbone of the more general \code{commutability} package.
}
\details{
  Key functions:
  \itemize{
    \item \code{predict_eqa()}: Model data and calculate estimated
    prediction intervals.
    \item \code{estimate_zeta_ols()}: Estimate \eqn{zeta} assuming a linear
    relationship between the IVD-MDs in comparison.
    \item \code{global_precision_estimates()}: Estimates repeatability
    statistics for both IVD-MDs in comparison.
    \item \code{bootstrap_ci()}: Estimating confidence intervals for
    repeatability statistics and \eqn{\zeta}.
    \item \code{merge_results()}: Merge individual components into a
    single object.
  }
  
}
\author{
Pernille Kjeilen Fauskanger

Maintainer: Pernille Kjeilen Fauskanger <pernille.fauskanger@noklus.no>
}
\references{
  Sandberg S. et al. (2023) for the IFCC Working Group on Commutability in Metrological Traceability, Recommendations for Setting a Criterion and Assessing Commutability of Sample Materials Used in External Quality Assessment/Proficiency Testing Schemes, \emph{Clinical Chemistry}. 69(11): 1227-1237.
  \url{https://doi.org/10.1093/clinchem/hvad135}
  
  Fauskanger P.K., et al. (2025) Quantification of Difference in Nonselectivity Between In Vitro Diagnostic Medical Devices. \emph{Biometrical Journal}. 67: e70032.
  \url{https://doi.org/10.1002/bimj.70032}
}
\keyword{ package }
\seealso{
  Useful links:
}
\examples{
  # 1
  print(1)
}
