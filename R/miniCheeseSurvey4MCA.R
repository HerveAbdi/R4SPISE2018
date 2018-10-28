# Description for the excel file
# miniCheeseSurvey4MCA.xlsx
# Hervé Abdi: July 12, 2018.
# for SPISE 2018 Workshop
#
#  miniCheeseSurvey4MCA Preambule ----
#'@title An excel file storing the data of
#' a subset of questions from a
#' survey questionnaire to be analyzed by MCA.
#'
#'@description \code{miniCheeseSurvey4MCA}:
#' An excel file
#' storing the data from an attitude survey:
#' 305 (French) consumers
#' described their attitude, consumption, and knowledge about
#'  cheese (Maroille).
#'  These data are analyzed with multiple correspondence analysis
#'  with the packages \code{ExPosition} and \code{PTCA4CATA}
#' (note that the package
#' \code{PTCA4CATA} needs to have been installed
#' from github with the
#' command \code{devtools::install_github('HerveAbdi/PTCA4CATA')}).
#'
#'@details
#'  The data are stored in an excel file containing 5 sheets.
#' The first sheet, called \code{Feuil1} contains the
#'  data for MCA: the rows are consumers and
#' columns are variables.
#' The other sheets (called \code{Feuil2} to \code{Feuil5})
#' provide additional information.
#' The questionnaire is provided as
#' \code{cheeseSurveyQuestionnaire}.
#' The questions kept in the mini version concern
#' \emph{Knowledge} (8 questions), \emph{Beliefs} (24 questions)
#' along with descriptors of the participants: Age,
#' Sex, CSP, size of the family.
#'
#' @references
#'  These data were used in a dissertation on
#'  Cheese consumption in France written by Menouar Nacef
#'  under the direction of Professor Sylvie Chollet.
#' @keywords datasets DistatisR
#' @author  Nacef, M., & Chollet, S.
#' @name cheeseSurvey4MCA
#' @section FileName: miniCheeseSurvey4MCA.xlsx
#' @section ReadingTheData:
#' To fetch this dataset use \code{system.file()} (see example below).
#' @examples
#' \dontrun{
#' # Note we need to have DistatisR installed for the example to run
#'path2file <- system.file("extdata",
#'        "miniCheeseSurvey4MCA.xlsx", package = "R4SPISE2018")
#'cheeseDATA <- PTCA4CATA::read.xls.CATA(path = path2file,
#'                                          sheet = "Data")
#'           }
NULL
# End of  miniCheeseSurvey4MCA----
#_____________________________________________________________________



