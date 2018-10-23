## ----note, include = FALSE, ECHO = FALSE, eval = FALSE-------------------
#  **NOTE:**
#  
#  This `pdf` was generated from the vignette
#  `R4SPISE2018::cheeseMCA` from the `R` Package
#  `R4SPISE2018`. Check the help for the
#  very last version of this document.

## ---- include = TRUE, echo = TRUE----------------------------------------
rm(list = ls())
graphics.off()

## ----setup, include = FALSE, ECHO = FALSE--------------------------------
# Important: Remember 
#     build the vignettes with devtools::build_vignettes()
knitr::opts_chunk$set( collapse = TRUE, comment = "#>")

## ---- eval = FALSE,ECHO = FALSE , include = FALSE------------------------
#  knitr::opts_knit$get()

## ----loadPackages--------------------------------------------------------
# Decomment all/some these lines if the packages are not installed
# devtools::install_github('HerveAbdi/PTCA4CATA')
# devtools::install_github('HerveAbdi/DistatisR')
# devtools::install_github('HerveAbdi/R4SPISE2018') # of course!
#  install.packages(prettyGraphs)
#  install.packages('Matrix')
#  install.packages('dplyr')
#  install.packages('gridExtra')
#  install.packages('grid')
#  install.packages('gtable')
#  install.packages('stringi')
#  load the libraries that we will need
suppressMessages(library(Matrix))
suppressMessages(library(DistatisR))
suppressMessages(library(PTCA4CATA))
suppressMessages(library(prettyGraphs))
suppressMessages(library(ExPosition))
suppressMessages(library(dplyr))
suppressMessages(library(gridExtra)) # to save a table as a graph
suppressMessages(library(grid))      # that will be saved in the
suppressMessages(library(gtable))    # powerpoint with the figures


## ----findDataPath--------------------------------------------------------
path2file <- system.file("extdata",
       "miniCheeseSurvey4MCA.xlsx", package = "R4SPISE2018")

## ----xls.datafile, echo=FALSE, fig.cap="The Data Excel File \\label{fig:spicesxl}", fig.height=3, fig.width=4, include=TRUE, out.width='70%'----

knitr::include_graphics('../man/figures/maroilesMCA.png')


## ----resdSortingData-----------------------------------------------------
rawData <- read.df.excel(path = path2file, sheet = 'Data')$df.data

## ----savexls-------------------------------------------------------------
saveFile <- file.copy(from = path2file, to = '~/Downloads/myDataFile.xlsx')

## ----maketmpData---------------------------------------------------------
temp.df <- dplyr::select(rawData, Q09_Sex:Q11_City, Q01_Know1:Q08_Know8, 
                   Q15_C01:Q38_C24)

## ----recodeSex-----------------------------------------------------------
# We recode sex as m/f
sex  =  car::recode(temp.df[,'Q09_Sex'], " 1 = 'm';2 = 'f'")

## ---- message = FALSE, warning = FALSE-----------------------------------
# correct responses
correct.Know = c(3,3,1,2,1,1,1,3)
# where are the knowledge questions
knowQuestions = (substr( colnames(temp.df),5,8)   == 'Know') 
# Get the number of correct amswers
correct.answers =  rowSums(matrix(correct.Know, nrow = nrow(temp.df), 
                   ncol = length(correct.Know), byrow = TRUE) ==
           temp.df[,which(knowQuestions)])

## ----histKnow, fig.width = 6, fig.height= 5------------------------------
hist(correct.answers, main = '', breaks = 0:8,  xlab = 'Number of Correct Answers')

## ----rescaleKnow---------------------------------------------------------
# Make only three categories use recode from package car
know =  car::recode(correct.answers, "0 = 1; 1 = 1; 2 = 1; 3 = 1; 4 = 2; 5 = 3; 6 = 3")

## ----rescaleAge----------------------------------------------------------
oriAge <- temp.df$Q10_Age
newAge  <- car::recode(oriAge, "1 = 1; 2 = 1; 3 = 2; 4 = 3; 5 = 4; 6 = 4")

## ---- fig.cap = 'Response distribution for Likert scales \\label{tab:Likert}', outwidth = '100%'----
knitr::kable( t(apply(rawData[,32:55],2,function(x){summary(as.factor(x))})))

## ----Lick2Keep-----------------------------------------------------------
Lick2Keep = c('Q15_C01','Q19_C05','Q21_C07',
               'Q24_C10','Q25_C11','Q32_C18')

## ----binLick-------------------------------------------------------------

Lickert <- temp.df[, colnames(temp.df) %in% Lick2Keep]
recl <- function(x){y = car::recode(x, "1 = 1; 2 = 1; 3 = 2; 4 = 2") }
rec.Lickert <- apply(Lickert,2,recl)

## ----cleanData-----------------------------------------------------------
cleanData <- cbind(sex,newAge,rawData$Q11_City,know,rec.Lickert)
colnames(cleanData) <- c('Sex','Age','City','Know',
                          'C01','C05','C07','C10','C11','C18')
# make sure that there is no NA left
cleanData <- cleanData[complete.cases(cleanData),]

## ----runMCA--------------------------------------------------------------
resMCA <- epMCA(cleanData, graphs = FALSE) 

## ----subset.L------------------------------------------------------------
 cleanData.Lille = cleanData[cleanData[,'City'] == 'Lille',]
 cleanData.Lille <- cleanData.Lille[, (colnames(cleanData.Lille) != 'City') ]

## ----runMCA.lille--------------------------------------------------------
resMCA.Lille <- epMCA(cleanData.Lille, graphs = FALSE) 

## ----subset.A------------------------------------------------------------
 cleanData.Angers <- cleanData[cleanData[,'City'] == 'Angers',]
 cleanData.Angers <- cleanData.Angers[, (colnames(cleanData.Angers) != 'City') ]

## ----runMCA.angers-------------------------------------------------------
 resMCA.Angers <- epMCA(cleanData.Angers, graphs = FALSE) 

## ----colors--------------------------------------------------------------
cJ <- resMCA$ExPosition.Data$cj
color4Var <- prettyGraphs::prettyGraphsColorSelection(ncol(cleanData))

## ----ctrVar--------------------------------------------------------------
# Extract the root names before the "."
varNames <- stringi::stri_extract(rownames(cJ),regex = '[^.]*')
varCtr.tmp <- aggregate(cJ ~ varNames, (cbind(varNames,cJ)),sum)
varCtr <- varCtr.tmp[,-1]
rownames(varCtr)    <- varCtr.tmp[,1]
rownames(color4Var) <- varCtr.tmp[,1]

## ----ctrVar.Tab----------------------------------------------------------
nFact <- min(5, ncol(cJ) - 1)
# knitr::kable(round( varCtr[,1:nFact]*1000 ) )
# save table as a graph
ctrTable <- tableGrob(round(varCtr[,1:nFact]*1000))
h <- grobHeight(ctrTable)
w <- grobWidth(ctrTable)
title <- textGrob("Variable Contributions",
                  y = unit(0.5,"npc") + 1.2*h, 
                  vjust = 0,
                  gp = gpar(fontsize=15))
TableWithTitle <- gTree(children = gList(ctrTable, title))
grid.draw(TableWithTitle)
a00.1.ctrTable  <- recordPlot()

## ------------------------------------------------------------------------
nM   <- nrow(cJ)
nVar <- nrow(color4Var)
col4Labels <- rep("",nM)

for (i in 1:nVar){
  lindex <- varNames %in% rownames(color4Var)[i]
  col4Labels[lindex] <- color4Var[i]
}

## ----screeMCA, fig.height=4, fig.width= 7--------------------------------
# 5.A. A scree plot for the RV coef. Using standard plot (PTCA4CATA)
scree.mca <- PlotScree(ev = resMCA$ExPosition.Data$eigs, 
                   p.ev = NULL, max.ev = NULL, alpha = 0.05,
                   col.ns = "#006D2C", col.sig = "#54278F",
                   title = "MCA. Explained Variance per Dimension", 
                   plotKaiser = FALSE,
                   color4Kaiser = "darkorchid4", 
                   lwd4Kaiser = 2.5)
a1.Scree <- recordPlot() # Save the plot

## ----createFjMap---------------------------------------------------------
axis1 = 1
axis2 = 2
Fj <- resMCA$ExPosition.Data$fj
# generate the set of maps
BaseMap.Fj <- createFactorMap(X = Fj , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              constraints = NULL,
                              title = 'MCA. Variables', 
                              col.points = col4Labels,
                              display.points = TRUE,
                              pch = 19, cex = 1,
                              display.labels = TRUE,
                              col.labels = col4Labels,
                              text.cex = 2.5, font.face = "bold",
                              font.family = "sans",
                              col.axes = "darkorchid",
                              alpha.axes = 0.2,
                              width.axes = 1.1,
                              col.background = adjustcolor("lavender",
                                                       alpha.f = 0.2),
                              force = 2, segment.size = 0)
# add labels
labels4MCA <- createxyLabels.gen(x_axis = axis1,
                                   y_axis = axis2,
               lambda = resMCA$ExPosition.Data$eigs,
               tau = resMCA$ExPosition.Data$t)
# make the maps
aa.1.BaseMap.Fj <- BaseMap.Fj$zeMap +  labels4MCA 
aa.2.BaseMapNoDot.Fj  <- BaseMap.Fj$zeMap_background +
                          BaseMap.Fj$zeMap_text + labels4MCA 

## ----plotaMap, fig.width= 8 , fig_width = '100%'-------------------------
print(aa.1.BaseMap.Fj)

## ----createFiMap---------------------------------------------------------
Fi <- resMCA$ExPosition.Data$fi
colCity <- c('darkblue', 'red4')
nI <- nrow(Fi)
col4I.City <- rep("",nI)
for (i in 1:length(colCity) ){
  lindex <- cleanData[,'City'] %in% unique(cleanData[,'City'])[i]
  col4I.City[lindex] <- colCity[i]
}
# generate the set of maps
BaseMap.Fi <- createFactorMap(X = Fi , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              constraints = NULL,
                              title = 'MCA. Variables', 
                              col.points = col4I.City,
                              alpha.points = .2,
                              display.points = TRUE,
                              pch = 19, cex = .8,
                              display.labels = TRUE,
                              col.labels = col4I.City,
                              text.cex = 2.5, font.face = "bold",
                              font.family = "sans",
                              col.axes = "darkorchid",
                              alpha.axes = 0.2,
                              width.axes = 1.1,
                              col.background = adjustcolor("lavender",
                                                       alpha.f = 0.2),
                              force = 2, segment.size = 0)
# make the maps
aa.5.BaseMapNoLabels.Fi  <- BaseMap.Fi$zeMap_background +
                          BaseMap.Fi$zeMap_dots + labels4MCA 

## ----plotaMapi, fig.width= 8---------------------------------------------
print(aa.5.BaseMapNoLabels.Fi)

## ----createFjMap.A-------------------------------------------------------
col4Labels.sub <- col4Labels[varNames != 'City']
axis1 = 1
axis2 = 2
Fj.Angers <- resMCA.Angers$ExPosition.Data$fj
# generate the set of maps
BaseMap.Fj.Angers <- createFactorMap(X = Fj.Angers , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              constraints = NULL,
                              title = 'Angers. MCA. Variables', 
                              col.points = col4Labels.sub,
                              display.points = TRUE,
                              pch = 19, cex = 1,
                              display.labels = TRUE,
                              col.labels = col4Labels.sub,
                              text.cex = 2.5, font.face = "bold",
                              font.family = "sans",
                              col.axes = "darkorchid",
                              alpha.axes = 0.2,
                              width.axes = 1.1,
                              col.background = adjustcolor("lavender",
                                                       alpha.f = 0.2),
                              force = 2, segment.size = 0)
# add labels
labels4MCA.Angers <- createxyLabels.gen(x_axis = axis1,
                                   y_axis = axis2,
               lambda = resMCA.Angers$ExPosition.Data$eigs,
               tau = resMCA.Angers$ExPosition.Data$t)
# make the maps
ca.1.BaseMap.Fj.Angers       <- BaseMap.Fj.Angers$zeMap +  labels4MCA.Angers 
ca.2.BaseMapNoDot.Fj.Angers  <- BaseMap.Fj.Angers$zeMap_background +
                          BaseMap.Fj.Angers$zeMap_text + labels4MCA.Angers 

## ----plotaMap.A, fig.width = 10, fig.height = 8--------------------------
print(ca.1.BaseMap.Fj.Angers)

## ----createFjMap.L-------------------------------------------------------
col4Labels.sub <- col4Labels[varNames != 'City']
axis1 = 1
axis2 = 2
Fj.Lille <- resMCA.Lille$ExPosition.Data$fj
# generate the set of maps
BaseMap.Fj.Lille <- createFactorMap(X = Fj.Lille , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              constraints = NULL,
                              title = 'Lille. MCA. Variables', 
                              col.points = col4Labels.sub,
                              display.points = TRUE,
                              pch = 19, cex = 1,
                              display.labels = TRUE,
                              col.labels = col4Labels.sub,
                              text.cex = 2.5, font.face = "bold",
                              font.family = "sans",
                              col.axes = "darkorchid",
                              alpha.axes = 0.2,
                              width.axes = 1.1,
                              col.background = adjustcolor("lavender",
                                                       alpha.f = 0.2),
                              force = 2, segment.size = 0)
# add labels
labels4MCA.Lille <- createxyLabels.gen(x_axis = axis1,
                                   y_axis = axis2,
               lambda = resMCA.Lille$ExPosition.Data$eigs,
               tau    = resMCA.Lille$ExPosition.Data$t)
# make the maps
cb.1.BaseMap.Fj.Lille       <- BaseMap.Fj.Lille$zeMap + labels4MCA.Lille 
cb.2.BaseMapNoDot.Fj.Lille  <- BaseMap.Fj.Lille$zeMap_background +
                          BaseMap.Fj.Lille$zeMap_text + labels4MCA.Lille 

## ----plotaMap.L, fig.width = 10, fig.height=8----------------------------
print(cb.1.BaseMap.Fj.Lille)

## ----saveGraphs, message = FALSE, warning = FALSE, error = FALSE, eval = FALSE----
#  name4Graphs = 'cheeseFrom2Cities.pptx'
#  list2Graphs <- PTCA4CATA::saveGraph2pptx(file2Save.pptx = name4Graphs,
#                   title = 'Attitudes toward Cheese',
#                   addGraphNames = TRUE)

## ----powerpoint,  message = FALSE, warning = FALSE, error = FALSE, eval = FALSE----
#  output:
#        powerpoint_presentation:
#             slide_level: 4

## ----vignettes,  message = FALSE, warning = FALSE, error = FALSE, eval = FALSE----
#  output:
#         rmarkdown::html_vignette:
#            toc: true
#            number_sections: true

