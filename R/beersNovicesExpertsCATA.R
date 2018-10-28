# Description for the excel file
# beersNovicesExpertsCATA.xlsx
# Hervé Abdi: July 8, 2018.
# for SPISE 2018 Workshop
#
#  beersNovicesExpertsCATA Preambule ----
#'@title An excel file storing the data from
#' a CATA task: Novices and Experts describe 9 beers.
#'
#'@description \code{beersNovicesExpertsCATA}.
#' An excel file
#' storing the data from a CATA task:
#' 80 Novices and 12 Experts evaluated 9 beers
#' described by 57 descriptors
#' that belong
#' to one of three blocks (Gustative, Emotion, or Context).
#' This excel file can be accessed with the \code{R}-command
#' \code{system.file()}, and can be
#' read by the function
#' \code{PTCA4CATA::read.xls.CATA} (note that the package
#' \code{PTCA4CATA} needs to have been installed first
#' from \code{Github} with the
#' command \code{devtools::install_github('HerveAbdi/PTCA4CATA')}).
#'
#'@details
#'  The data are stored in an excel file containing 5  sheets.
#' The first sheet, called \code{DataCATA} contains the
#' CATA data: the rows are the products and the
#' columns are the descriptors. In this sheet
#' the data are organized and
#' coded as indicated in the help for \code{PTCA4CATA::read.xls.CATA}
#' and for  \code{PTCA4CATA::OrangeJuiceSortingRawData}.
#'  ..
#' The second sheet, called \code{DescriptionJudges} contains the
#' description of the Assessors; The explanation of the description
#' is given in the  third sheet called
#' \code{Legend4DescriptionJudges}.
#' The fourth sheet, called \code{Categories4Descriptors}
#' gives the category (C/E/G) of the descriptors;
#' The explanation of the description
#' is given in the fifth sheet called
#' \code{Legend4Categories} (C = Context,
#' E = Emotion, G = Gustative).
#'
#'
#' @references
#'  These data were used for a keynote address given
#'  at the 2014 meeting of the Sensometric society
#'  (see \url{http://www.sensometric.org/page-1757209})
#'   Abdi, H., Chollet, S., Valentin, D., & Lelievre, M.
#'  \emph{ 	Partial Triadic Correspondence Analysis of BeTA
#'      (Best that Apply) Data}. Chicago, 2014.
#' @keywords datasets DistatisR
#' @author  Chollet, S., Valentin, D., Lelievre, M., & Abdi, H.
#' @name beersNovicesExpertsCATA
#' @section FileName: beersNovicesExpertsCATA.xlsx
#' @section ReadingTheData:
#' To fetch this dataset use \code{system.file()} (see example below).
#' @examples
#' \dontrun{
#' # Note we need to have DistatisR installed for the example to run
#'path2file <- system.file("extdata",
#'        "beersNovicesExpertsCATA.xlsx", package = "R4SPISE2018")
#'beersCATA <- PTCA4CATA::read.xls.CATA(path = path2file,
#'                                          sheet = "DataCATA")
#'           }
NULL
# End of  beersNovicesExpertsCATA----
#_____________________________________________________________________



