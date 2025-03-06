#' Fictive Clinical Sample Data
#'
#' @format ## `test_data`
#' A data frame with 60 rows and 4 columns:
#' \describe{
#'   \item{SampleID}{ID of clinical sample}
#'   \item{ReplicateID}{ID of replicated measurements}
#'   \item{MP_A}{First IVD-MD in the comparison}
#'   \item{MP_B}{Second IVD-MD in the comparison}
#' }
#' @source clsi
"test_data"

#' Glucose Serum Clinical Sample Meausurements
#'
#' @format ## `test_data`
#' A data frame with 75 rows and 6 columns:
#' \describe{
#'   \item{SampleID}{ID of clinical sample}
#'   \item{ReplicateID}{ID of replicated measurements}
#'   \item{Advia}{Measurement results from an IVD-MD named Advia.}
#'   \item{Alinity}{Measurement results from an IVD-MD named Alinity.}
#'   \item{Cobas}{Measurement results from an IVD-MD named Cobas.}
#'   \item{Vitros}{Measurement results from an IVD-MD named Vitros.}
#' }
#' @source Noklus
"glucose_cs_data"

#' Glucose Serum External Quality Assessment Material Sample Meausurements
#'
#' @format ## `test_data`
#' A data frame with 9 rows and 6 columns:
#' \describe{
#'   \item{SampleID}{ID of EQA material sample}
#'   \item{ReplicateID}{ID of replicated measurements}
#'   \item{Advia}{Measurement results from an IVD-MD named Advia.}
#'   \item{Alinity}{Measurement results from an IVD-MD named Alinity.}
#'   \item{Cobas}{Measurement results from an IVD-MD named Cobas.}
#'   \item{Vitros}{Measurement results from an IVD-MD named Vitros.}
#' }
#' @source Noklus
"glucose_eqam_data"
