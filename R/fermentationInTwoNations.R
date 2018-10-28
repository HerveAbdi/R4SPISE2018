# Description for the excel file
# questionsCheeseSurvey.docx
# Hervé Abdi: July 12, 2018.
# for SPISE 2018 Workshop
#
#  fermentationInTwoNations Preamble ----
#'@title An \code{.xlsx} and a \code{.pdf} file
#'giving the answers (Excel file) to  a
#'questionnaire (pdf file) about fermented products.
#'
#' @description \code{fermentationInTwoNations}:
#' An \code{.xlsx} file and a \code{.pdf} file
#'giving the results (Excel file) to a
#'questionnaire (pdf file) about fermented products.
#' These results come from
#' a study about the
#' social representations of fermented products
#' in France and Vietnam.
#' The excel file name is 
#' \code{fermentedFoodSurvey.xlsx}
#' (the data are stored in the sheet" \code{data4MCA}),
#' and the questionnaire is
#' in the file \code{questionnaireFermentation.pdf}.
#' 
#'@details
#' The excel file includes data from 373 participants 
#' (220 French and 153, Vietnamese participants) to the
#' original study plus 30 "new" participants (students
#' from the 2018 SPISE Workshop).
#' A detailled analysis of these data is given in the vignette
#' \code{fermentationIn2NationsMCA}
#' from the \code{R4SPISE2018} package.
#'
#' @references
#'  These data were used in a project directed
#'  by Professor Sylvie Chollet.
#' @keywords datasets DistatisR
#' @author   Chollet, S.
#' @name fermentationInTwoNations
#' @section FileNames: \code{fermentedFoodSurvey.xlsx} and
#' \code{questionnaireFermentation.pdf}
#' @section ReadingTheData:
#' To fetch this dataset use \code{system.file()} (see example below).
#' @examples
#' \dontrun{
#' # to get the path of the files:
#' path2file.pdf <- system.file("extdata",
#'        "questionnaireFermentation.pdf", package = "R4SPISE2018")
#' path2file.xlsx <- system.file("extdata",
#'        "fermentedFoodSurvey.xlsx", package = "R4SPISE2018")
#'           }
NULL
# End of  fermentationInTwoNations----
#_____________________________________________________________________



