# Description for the excel file
# multiculturalSortingSpices.xlsx
# Hervé Abdi: July 7, 2018.
# for SPISE 2018 Workshop
#
#  multiculturalSortingSpices Preambule ----
#' @title An excel file:
#' 61 assessors from 5 countries sorted 16 spice samples
#'
#' @description
#'\code{multiculturalSortingSpices}. An excel file
#' storing the data from a sorting task:
#' 62 participants from 5 different countries
#' (USA, France, India, Spain, and Vietnam)
#' sorted 16 different spices
#' (including 6 mixtures of spices).
#' This excel file can be accessed with the \code{R}-command
#' \code{system.file()}, and can be
#' read by the function
#' \code{DistatisR::read.df.excel} (note that the package
#' \code{DistatisR} needs to have been installed first
#' from \code{Github} as 
#' \code{devtools::install_github('HerveAbdi/DistatisR')}).
#'
#'@details
#'  The data are stored in an excel file containing 3 sheets.
#' The first sheet, called \code{DataSort} contains the
#' sorting data: the rows are the products and the
#' columns are the assessors. The data are organized and
#' coded as indicated in the help for \code{distatis}
#' (see also the Abdi et al.' 2007 Distatis paper).
#' Some spices are mixtures of several spices and the names are
#' short names to make graphs easier to read (the long names
#' are given in the sheet \code{legend4Spices}).
#' The first letter of the name of the assessors gives
#' their nationalities: S for Spanish (10),
#' V for Vietnamese (6), I for Indian (15),
#' F for French (21), and A for American (9).
#' @docType data
#' @references
#' Part of these data (i.e., the French sample) is described
#' and analyzed in
#'  Chollet, S., Valentin, D., & Abdi, H. (2014).
#'  Free sorting task. In P.V. Tomasco & G. Ares (Eds),
#'  \emph{Novel Techniques in Sensory Characterization
#'  and Consumer Profiling}.
#'  Boca Raton: Taylor and Francis. pp 207-227.
#'
#'  The format of the sorting task is the same as the beer example in
#'  Abdi, H., Valentin, D., Chollet, S., & Chrea, C. (2007).
#'  Analyzing assessors and products in sorting tasks:
#'  DISTATIS, theory and applications.
#'  \emph{Food Quality and Preference, 18}, 627–640.
#' @keywords datasets DistatisR
#' @author  Chollet, S., Valentin, D., & Abdi, H.
#' @name multiculturalSortingSpices
#' @section FileName: multiculturalSortingSpices.xlsx
#' @section ReadingTheData:
#' To fetch this dataset use \code{system.file()} (see example below).
#' @examples
#' \dontrun{
#' # Note we need to have DistatisR installed for the example to run
#'path2file <- system.file("extdata",
#'        "multiculturalSortingSpices.xlsx", package = "R4SPISE2018")
#'spiceDataSort <- DistatisR::read.df.excel(path = path2file,
#'                                          sheet = "DataSort")
#'           }
NULL
# End of  multiculturalSortingSpices----
#_____________________________________________________________________



