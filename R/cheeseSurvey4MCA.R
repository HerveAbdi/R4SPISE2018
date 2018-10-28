# Description for the excel file
# cheeseSurvey4MCA.xlsx
# Hervé Abdi: July 12, 2018.
# for SPISE 2018 Workshop
#
#  cheeseSurvey4MCA Preambule ----
#'@title An excel file storing the answers
#' a survey questionnaire to be analyzed by MCA.
#'
#'@description \code{cheeseSurvey4MCA}:
#' An excel file
#' storing the data from an attitude survey:
#' 305 (French) consumers
#' described their attitude, consumption, and knowledge about
#' a 
#'  cheese (Maroille) made in northen France.
#'  These data are analyzed with multiple correspondence analysis
#'  with the packages \code{ExPosition} and \code{PTCA4CATA}
#' (note that the package
#' \code{PTCA4CATA} needs to have been installed
#' from \code{Github} with the
#' command \code{devtools::install_github('HerveAbdi/PTCA4CATA')}).
#'
#'@details
#'  The data are stored in an excel file containing 5  sheets.
#' The first sheet, called \code{Feuil1} contains the
#'  data for MCA: the rows are consumers and
#' columns are variables.
#' The other sheets (called \code{Feuil2} to \code{Feuil5})
#' provide additional information.
#' The questionnaire is provided as
#' \code{cheeseSurveyQuestionnaire}.
#'
#' @references
#'  These data were used in a dissertation on
#'  Cheese consumption in France written by  Menouar Nacef
#'  under the direction of Professor Sylvie Chollet.
#' @keywords datasets R4SPISE2014
#' @author  Nacef, M., & Chollet, S.
#' @name cheeseSurvey4MCA
#' @section FileName: cheeseSurvey4MCA.xlsx
#' @section ReadingTheData:
#' To fetch this dataset use \code{system.file()} (see example below).
#' @examples
#' \dontrun{
#' # Note we need to have DistatisR installed for the example to run
#'path2file <- system.file("extdata",
#'        "cheeseSurvey4MCA.xlsx", package = "R4SPISE2018")
#'cheeseDATA <- PTCA4CATA::read.xls.CATA(path = path2file,
#'                                          sheet = "Feuil1")
#'           }
NULL
# End of  cheeseSurvey4MCA ----
#_____________________________________________________________________



