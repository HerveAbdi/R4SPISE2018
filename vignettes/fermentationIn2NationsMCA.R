## ---- eval = TRUE, ECHO = FALSE , include = FALSE------------------------
knitr::opts_chunk$set(collapse = TRUE, fig.width = 9, comment = "#>")
#knitr::opts_chunk$set(fig.path = '../test4Spise2018/lesFigs4Ferment/') 
#knitr::opts_chunk$set(tex.path = '../test4Spise2018/lesFigs4Ferment/') 
# Knitr options here
knitr::opts_knit$get()

## ----pdfIt, echo = FALSE, warning = FALSE, eval = FALSE------------------
#  # This code can be used to generate pdf or words
#  rmarkdown::render('../test4Spise2018/fermentationMCA.Rmd',
#                     c('bookdown::pdf_document2','bookdown::html_document2') ) ,
#                     intermediates_dir = '../../test4Spise2018/stuff4Ferment/')
#  rmarkdown::render('../test4Spise2018/fermentationMCA.Rmd',
#                     c('bookdown::pdf_document2','bookdown::html_document2'),
#                     output_dir = '../../test4Spise2018/stuff4Ferment/')
#  # Previous command cannot pass options such as keep_tex
#  # to do so run only one format at a time such as:
#  rmarkdown::render("fermentationIn2NationsMCA.Rmd",
#                     bookdown::pdf_document2(keep_tex = TRUE,
#                     toc_depth = 3)  ,
#                     output_dir = '../../test4Spise2018/stuff4Ferment/')
#  # rmarkdown::render('../test4Spise2018/fermentationMCA.Rmd','tufte::tufte_html')
#  # Not tufte works only with 2 levels of section
#  # rmarkdown::render('../test4Spise2018/fermentationMCA.Rmd',tufte::tufte_handout(toc_depth = 2))
#  # rmarkdown::render('../test4Spise2018/fermentationMCA.Rmd','html_notebook')
#  # better
#  # xaringan::inf_mr()
#  # pretty format
#  rmarkdown::render('../../test4Spise2018/fermentationMCA.Rmd',
#                     prettydoc::html_pretty(keep_tex = "hpstr",
#                     toc_depth = 3) )
#  rmarkdown::render('../../test4Spise2018/fermentationMCA.Rmd',
#                     rmdformats::html_clean(keep_tex = TRUE
#                      ) )

## ---- include = TRUE, echo = TRUE----------------------------------------
rm(list = ls())
graphics.off()

## ----setup, include = FALSE, ECHO = FALSE--------------------------------
# Important: Remember 
#     build the vignettes with devtools::build_vignettes()
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 9,
  fig.width =  8
)

## ---- eval = FALSE,ECHO = FALSE , include = FALSE------------------------
#  knitr::opts_knit$get()

## ----loadPackages--------------------------------------------------------
# Decomment all/some these lines if the packages are not installed
#_____________________________________________________________________
# # October 19, 2018. Temporary fix for an Rstudio problem
# # if this error message is generated: 
# build_site()
# Error in eval(substitute(expr), data, enclos = parent.frame()) :
# is.character(repos) is not TRUE
# # from Rstudio team: as a temporary fix use these two lines
# repos <- getOption("repos")
# options(repos = setNames(as.character(repos), names(repos)))
#_____________________________________________________________________
#
# devtools::install_github('HerveAbdi/PTCA4CATA')
# devtools::install_github('HerveAbdi/DistatisR')
# devtools::install_github('HerveAbdi/data4PCCAR')
# devtools::install_github('HerveAbdi/R4SPISE2018') # of course!
#  install.packages('prettyGraphs')
#  install.packages('Matrix')
#  install.packages('dplyr')
#  install.packages('gridExtra')
#  install.packages('grid')
#  install.packages('gtable')
#  install.packages('stringi')
#  install.packages('printr')
#  install.packages('kableExtra')
#  load the libraries that we will need
suppressMessages(library(Matrix))
suppressMessages(library(prettyGraphs))
suppressMessages(library(ExPosition))
suppressMessages(library(InPosition))
suppressMessages(library(DistatisR))
suppressMessages(library(PTCA4CATA))
suppressMessages(library(data4PCCAR))
suppressMessages(library(dplyr))
suppressMessages(library(gridExtra))    # to save a table as a graph
suppressMessages(library(grid))         # that will be saved in the
suppressMessages(library(gtable))       # powerpoint with the figures
# suppressMessages(library(printr))     # To pretty print tables 
# suppressMessages(library(kableExtra)) # To pretty print tables 


## ----name4pptx-----------------------------------------------------------
name4Graphs = 'fermentationFrom2Countries_withJsup.pptx'

## ----filename, eval = TRUE-----------------------------------------------
file2read.name     <- 'fermentedFoodSurvey.xlsx' # xls data file name
path2file <- system.file("extdata", file2read.name, package = "R4SPISE2018")

## ----sheetName-----------------------------------------------------------
sheetName4Data     <- 'dataMCA' # sheet name for the data

## ----xls.datafile, echo=FALSE, fig.cap="The Data Excel File \\label{fig:spicesxl}", fig.height=3, fig.width=4, include=TRUE, out.width='70%'----

knitr::include_graphics('../man/figures/fermentationXlsFile.png')


## ----resdSortingData-----------------------------------------------------
rawData.tmp <- DistatisR::read.df.excel(path = path2file, 
                             sheet = sheetName4Data)$df.data
# 

## ------------------------------------------------------------------------
# Transform the data into factors
rawData <- apply(rawData.tmp,2,as.factor)

## ----summary-------------------------------------------------------------
 knitr::kable( t( summary(rawData) ) )

## ----maketmpData---------------------------------------------------------
temp.df <- dplyr::select(as.data.frame(rawData), 
                         supVar,
                         sex , age, occupation, nation, frequency, 
                         Wpreservation, Wquality, Whealth, Wtaste, 
                         fish, cereal,  vegetable, fruit, meat,
                         quality, industrial, preservation, health,  
                         expensive, taste, trust, clear, innovation)

## ----recodeAge-----------------------------------------------------------
temp.df[,'age'] <- plyr::mapvalues(temp.df[,'age'], 
                          from = c("<25", ">25"), to = c("Y", "O"))
temp.df[,'occupation'] <- plyr::mapvalues(temp.df[,'occupation'], 
                          from = c("FM", "other"), to = c("F", "O"))
temp.df[,'nation'] <- plyr::mapvalues(temp.df[,'nation'], 
                          from = c("F", "VN"), to = c("F", "V"))

## ----cleanData-----------------------------------------------------------
cleanData.tmp <- temp.df
cleanData.tmp <- cleanData.tmp[complete.cases(cleanData.tmp),]
cleanData.allVar <- cleanData.tmp[cleanData.tmp[,1] == 'A', 2:ncol(cleanData.tmp)]
cleanData.varSup <- cleanData.allVar[,1:4]
cleanData        <- cleanData.allVar[,5:ncol(cleanData.allVar)]
cleanData.sup <- cleanData.tmp[cleanData.tmp[,1] == 'S', 6:ncol(cleanData.tmp)]

## ----runMCA--------------------------------------------------------------
resMCA <- epMCA(cleanData, graphs = FALSE) 

## ----runMCA.sup----------------------------------------------------------
# recode the factors as set of 0/1 variables
testclean <- makeNominalData(rbind(cleanData,cleanData.sup))
clean.Sup <-  testclean[cleanData.tmp[,1] == 'S',]
# barycentric code for nation
#clean.Sup[,(colnames(testclean) %in% 'nation.F')] <- .5
#clean.Sup[,(colnames(testclean) %in% 'nation.V')] <- .5
#
resMCA.sup <- supplementaryRows(SUP.DATA = clean.Sup, res = resMCA)
colnames(resMCA.sup$fii) <- paste0('Dimension ',
                                1:ncol(resMCA.sup$fii))

## ----runMCA.varsup-------------------------------------------------------
#
resMCA.varSup <- supplementaryCols(SUP.DATA = makeNominalData(cleanData.varSup),
                                        res = resMCA)
colnames(resMCA.varSup$fjj) <- paste0('Dimension ',
                                1:ncol(resMCA.varSup$fjj))

## ----inferences, message = FALSE, warning = FALSE, results= FALSE--------
resMCA.inf <- epMCA.inference.battery(cleanData, graphs = FALSE)

## ----screeMCA, fig.height=4, fig.width= 7--------------------------------
scree.mca <- PlotScree(ev = resMCA$ExPosition.Data$eigs, 
                    title = "MCA. Explained Variance per Dimension")
b0001a.Scree <- recordPlot() # Save the plot

## ----screeMCA.inf, fig.height=4, fig.width= 7----------------------------
scree.mca <- PlotScree(ev = resMCA$ExPosition.Data$eigs, 
               p.ev = resMCA.inf$Inference.Data$components$p.vals, 
               plotKaiser = TRUE,
               title = "MCA. Explained Variance per Dimension")
b0001b.Scree <- recordPlot() # Save the plot

## ----colors--------------------------------------------------------------
cJ <- resMCA$ExPosition.Data$cj
color4Var <- prettyGraphs::prettyGraphsColorSelection(ncol(cleanData))

## ----phi2----------------------------------------------------------------
# Pseudo Heat Map. Correlation ----
# We need correlation to compare with PCA
corrMatBurt.list <- phi2Mat4BurtTable(cleanData)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corr4MCA.r <- corrplot::corrplot(
         as.matrix(corrMatBurt.list$phi2.mat^(1/2)),# to get correlations 
         method="color", col=col(200),  
         type="upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = color4Var, 
         tl.cex = .9,
         tl.srt = 45, #Text label color and rotation
         number.cex = .8,
         diag = TRUE # needed to have the color of variables correct
         )
# dev.new()
a0000.corMat.phi <- recordPlot()

## ----ctrVar--------------------------------------------------------------
varCtr <- data4PCCAR::ctr4Variables(cJ) 
rownames(color4Var) <- rownames(varCtr)

## ----ctrVar.Tab----------------------------------------------------------
nFact <- min(5, ncol(cJ) - 1)
#knitr::kable(round( varCtr[,1:nFact]*1000 ) )
# save table as a graph
ctrTable <- tableGrob(round(varCtr[,1:nFact]*1000))
h <- grobHeight(ctrTable)
w <- grobWidth(ctrTable)
title <- textGrob("Variable Contributions",
                 y = unit(0.5,"npc") + 0.92*h, 
                 # fine tune the position of the title 
                  just = "centre",
                  gp = gpar(fontsize = 14))
TableWithTitle <- gTree(children = gList(ctrTable, title))

## ---- fig.height= 10, fig.width= 6, fig.cap = 'Variable Contributions (per mille). \\label{fig:varCtr}', collapse = TRUE, fig.margin = TRUE----
# Note: Potential problems with grid.draw(). If it does not plot
# recordPlot() will fail and the graph will not be saved in the powerpoint
# and will generate a strange error message
grid.draw(TableWithTitle)
#a0000.2.ctrTable  <- recordPlot()

## ----printr,  echo = FALSE, message = FALSE------------------------------
# As an alternative we print the contributions with a combination
#of `kable` and `printr` as:
laTable <- round(varCtr[,1:nFact]*1000)
# knitr::kable(round(varCtr[,1:nFact]*1000), caption = 'Variable Contributions')
#    %>%
#   kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
#  add_header_above(c(" ", "Dimensions" = nFact))
  

## ----color4Levels--------------------------------------------------------
col4Levels <- data4PCCAR::coloringLevels(
                       rownames(resMCA$ExPosition.Data$fj), color4Var)
col4Labels <- col4Levels$color4Levels

## ----ctrV1---------------------------------------------------------------
varCtr1 <- varCtr[,1]
names(varCtr1) <- rownames(varCtr)
a0005.Var.ctr1  <- PrettyBarPlot2(varCtr1,
                    main = 'Variable Contributions: Dimension 1',
                                ylim = c(-.05, 1.2*max(varCtr1)),
                                font.size = 5,
                                threshold = 1 / nrow(varCtr),
                                color4bar = gplots::col2hex(color4Var)
)
print(a0005.Var.ctr1)

## ----ctrV2---------------------------------------------------------------
varCtr2 <- varCtr[,2]
names(varCtr2) <- rownames(varCtr)
a0006.Var.ctr2  <- PrettyBarPlot2(varCtr2,
                    main = 'Variable Contributions: Dimension 2',
                                ylim = c(-.05, 1.2*max(varCtr2)),
                                threshold = 1 / nrow(varCtr),
                                font.size = 5,
                                color4bar = gplots::col2hex(color4Var)
)
print(a0006.Var.ctr2)

## ----ctrV3---------------------------------------------------------------
varCtr3 <- varCtr[,3]
names(varCtr3) <- rownames(varCtr)
a0006.Var.ctr3  <- PrettyBarPlot2(varCtr3,
                    main = 'Variable Contributions: Dimension 3',
                                ylim = c(-.05, 1.2*max(varCtr2)),
                                threshold = 1 / nrow(varCtr),
                                font.size = 5,
                                color4bar = gplots::col2hex(color4Var)
)
print(a0006.Var.ctr3)

## ----ctrV12--------------------------------------------------------------
ctrV12 <- PTCA4CATA::createFactorMap(X =  varCtr, 
                        title = "Variable Contributions", 
                        col.points = color4Var,
                        col.labels = color4Var,
                        alpha.points = 0.5,
                        cex = 2.5, 
                        alpha.labels = 1, 
                        text.cex = 4,
                        font.face = "plain", 
                        font.family = "sans")

ctr.labels <- createxyLabels.gen(
  1,2, lambda = resMCA$ExPosition.Data$eigs,
  tau = resMCA$ExPosition.Data$t
)
a0007.Var.ctr12  <- ctrV12$zeMap  + ctr.labels
#
print(a0007.Var.ctr12)


## ----getCtr12------------------------------------------------------------
absCtrVar <- as.matrix(varCtr) %*% diag(resMCA$ExPosition.Data$eigs)
varCtr12  <- (absCtrVar[,1] + absCtrVar[,2]) / 
   (resMCA$ExPosition.Data$eigs[1] + resMCA$ExPosition.Data$eigs[2])
importantVar <- (varCtr12 >=  1 / length(varCtr12))
col4ImportantVar <- color4Var
col4NS <- 'gray90' 
col4ImportantVar[!importantVar] <- col4NS

## ----ctrV12.ns-----------------------------------------------------------
ctrV12.imp <- PTCA4CATA::createFactorMap(X =  varCtr, 
                        title = "Important Variables: Contributions", 
                        col.points = col4ImportantVar,
                        col.labels = col4ImportantVar,
                        alpha.points = 0.5,
                        cex = 2.5, 
                        alpha.labels = 1, 
                        text.cex = 4,
                        font.face = "plain", 
                        font.family = "sans")
a0008.Var.ctr12.imp  <- ctrV12.imp$zeMap  + ctr.labels
#
print(a0008.Var.ctr12.imp)


## ----getCtr23------------------------------------------------------------
#absCtrVar <- as.matrix(varCtr) %*% diag(resMCA$ExPosition.Data$eigs)
varCtr23  <- (absCtrVar[,3] + absCtrVar[,2]) / 
   (resMCA$ExPosition.Data$eigs[3] + resMCA$ExPosition.Data$eigs[2])
importantVar23 <- (varCtr23 >=  1 / length(varCtr23))
col4ImportantVar23 <- color4Var
col4NS <- 'gray90' 
col4ImportantVar23[!importantVar23] <- col4NS

## ----ctrV23.ns-----------------------------------------------------------
ctrV23.imp <- PTCA4CATA::createFactorMap(X =  varCtr,
                                         axis1 = 3, axis2 = 2,
                        title = "Important Variables: Contributions 3 * 2", 
                        col.points = col4ImportantVar23,
                        col.labels = col4ImportantVar23,
                        alpha.points = 0.5,
                        cex = 2.5, 
                        alpha.labels = 1, 
                        text.cex = 4,
                        font.face = "plain", 
                        font.family = "sans")
ctr.labels23 <- createxyLabels.gen(
  3,2, lambda = resMCA$ExPosition.Data$eigs,
  tau = resMCA$ExPosition.Data$t
)
a0009.Var.ctr23.imp  <- ctrV23.imp$zeMap  + ctr.labels23
#
print(a0009.Var.ctr23.imp)


## ----BR4var--------------------------------------------------------------
# Get the pseudo Bootstrap Rqtios
BrLevels <- resMCA.inf$Inference.Data$fj.boots$tests$boot.ratios
wJ       <- 1 / resMCA.inf$Fixed.Data$ExPosition.Data$W
nIter    <- 1000
Br4Variables <- data4PCCAR::BR4varMCA(BrLevels, wJ, nIter) 

## ----BR41----------------------------------------------------------------
VarBR1 <- Br4Variables$pseudoBR.pos[,1]
c0010.Var.br1  <- PrettyBarPlot2(VarBR1,
                    main = 'Variable Pseudo Bootstrap Ratios: Dimension 1',
                               ylim = 2,
                                threshold = 2,
                                font.size = 5,
                                color4bar = gplots::col2hex(color4Var)
)
print(c0010.Var.br1)

## ----BR42----------------------------------------------------------------
VarBR2 <- Br4Variables$pseudoBR.pos[,2]
c0011.Var.br2  <- PrettyBarPlot2(VarBR2,
                    main = 'Variable Pseudo Bootstrap Ratios: Dimension 2',
                               ylim = 2,
                                threshold = 2,
                                font.size = 5,
                                color4bar = gplots::col2hex(color4Var)
)
print(c0011.Var.br2)

## ----BR43----------------------------------------------------------------
VarBR3 <- Br4Variables$pseudoBR.pos[,3]
c0012.Var.br3  <- PrettyBarPlot2(VarBR3,
                    main = 'Variable Pseudo Bootstrap Ratios: Dimension 3',
                               ylim = 2,
                               threshold = 2,
                               font.size = 5,
                               color4bar = gplots::col2hex(color4Var)
)
print(c0012.Var.br3)

## ----createFjMap---------------------------------------------------------
axis1 = 1
axis2 = 2
Fj <- resMCA$ExPosition.Data$fj
# generate the set of maps
BaseMap.Fj <- createFactorMap(X = Fj , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              title = 'MCA. Variables', 
                              col.points = col4Labels, cex = 1,
                              col.labels = col4Labels, text.cex = 2.5,
                              force = 2)
# add labels
labels4MCA <- createxyLabels.gen(x_axis = axis1, y_axis = axis2,
               lambda = resMCA$ExPosition.Data$eigs,
               tau = resMCA$ExPosition.Data$t)
# make the maps
b0002.BaseMap.Fj <- BaseMap.Fj$zeMap + labels4MCA 
b0003.BaseMapNoDot.Fj  <- BaseMap.Fj$zeMap_background +
                          BaseMap.Fj$zeMap_text + labels4MCA 

## ----plotaMap, fig.width= 8 , fig_width = '100%'-------------------------
print(b0002.BaseMap.Fj)

## ----mapJ-grey-----------------------------------------------------------
col4Levels.imp <- data4PCCAR::coloringLevels(rownames(Fj),
                                             col4ImportantVar)
BaseMap.Fj.imp <- createFactorMap(X = Fj , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              title = 'MCA. Important Variables', 
                      col.points = col4Levels.imp$color4Levels, 
                              cex = 1,
                      col.labels = col4Levels.imp$color4Levels, 
                              text.cex = 2.5,
                              force = 2)
b0010.BaseMap.Fj <- BaseMap.Fj.imp$zeMap + labels4MCA 
print(b0010.BaseMap.Fj)

## ----adLines-------------------------------------------------------------
lines4J <- addLines4MCA(Fj, col4Var = col4Levels.imp$color4Variables, size = .7)
 b0020.BaseMap.Fj <-  b0010.BaseMap.Fj + lines4J
 print( b0020.BaseMap.Fj)

## ----someLines-----------------------------------------------------------
zeNames          <- getVarNames(rownames(Fj)) 
importantsLabels <- zeNames$stripedNames %in% zeNames$variableNames[importantVar]
Fj.imp <- Fj[importantsLabels,]
lines4J.imp <- addLines4MCA(Fj.imp, 
                            col4Var = col4Levels$color4Variables[which(importantVar)], 
                            size = .9, linetype = 3, alpha = .5)
 b0021.BaseMap.Fj <-  b0020.BaseMap.Fj + lines4J.imp
 print( b0021.BaseMap.Fj)

## ----mapJ23-grey---------------------------------------------------------
col4Levels23.imp <- data4PCCAR::coloringLevels(rownames(Fj),
                                             col4ImportantVar23)
axis3 = 3
BaseMap.Fj23.imp <- createFactorMap(X = Fj , # resMCA$ExPosition.Data$fj,
                              axis1 = axis3, axis2 = axis2,
                              title = 'MCA. Important Variables. Dimensions 2 & 3', 
                      col.points = col4Levels23.imp$color4Levels, 
                              cex = 1,
                      col.labels = col4Levels23.imp$color4Levels, 
                              text.cex = 2.5,
                              force = 2)
labels4MCA23 <- createxyLabels.gen(x_axis = axis1, y_axis = axis2,
               lambda = resMCA$ExPosition.Data$eigs,
               tau = resMCA$ExPosition.Data$t)
b0030.BaseMap.Fj23 <- BaseMap.Fj23.imp$zeMap + labels4MCA23 

# zeNames          <- getVarNames(rownames(Fj)) 
importantsLabels23 <- zeNames$stripedNames %in% zeNames$variableNames[importantVar23]
Fj23.imp <- Fj[importantsLabels23,]
lines4J23.imp <- addLines4MCA(Fj23.imp, 
                    col4Var = col4Levels$color4Variables[
                               which(importantVar23)],
                    axis_h = axis3,
                    axis_v = axis2,
                    size = .9, linetype = 3, alpha = .5)
 b0031.BaseMap.Fj23 <-  b0030.BaseMap.Fj23 + lines4J23.imp
 print( b0031.BaseMap.Fj23)

## ----mapvarSup-----------------------------------------------------------
col4VarSup <- prettyGraphs::prettyGraphsColorSelection(ncol(cleanData.varSup))
Fj.sup <- resMCA.varSup$fjj
col4Levels.sup <- data4PCCAR::coloringLevels(rownames(Fj.sup), col4VarSup)
BaseMap.Fj.sup <- createFactorMap(X = Fj.sup , # resMCA$ExPosition.Data$fj,
                axis1 = axis1, axis2 = axis2,
                constraints  = BaseMap.Fj$constraints, # to get same size
                title = 'MCA. Supplementary and Important Variables', 
                col.points = col4Levels.sup$color4Levels, 
                              cex = 1,
                col.labels = col4Levels.sup$color4Levels, 
                text.cex = 2.5,
                force = 2)
lines4J.sup <- addLines4MCA(Fj.sup, 
                  col4Var = col4Levels.sup$color4Variables, size = .7)
b0030.Sup.Fj <- BaseMap.Fj.sup$zeMap + 
                     BaseMap.Fj.imp$zeMap_dots + 
                     BaseMap.Fj.imp$zeMap_text +
                     labels4MCA + 
                     lines4J + lines4J.sup
print(b0030.Sup.Fj)

## ----mapvarSup.only------------------------------------------------------
b0031.Sup.Fj.only <- BaseMap.Fj.sup$zeMap + 
                     BaseMap.Fj.imp$zeMap_dots + 
                     labels4MCA + 
                      lines4J.sup
print(b0031.Sup.Fj.only)

## ----BR1-----------------------------------------------------------------

c0001.Levels.BR  <- PrettyBarPlot2(
      resMCA.inf$Inference.Data$fj.boots$tests$boot.ratios[,1], # BR
                    main = 'Bootstrap Ratios for Columns : Dimension 1',
                             threshold = 2,
                             color4bar = gplots::col2hex(col4Labels)
)
print(c0001.Levels.BR)

## ----createFiMap---------------------------------------------------------
Fi <- resMCA$ExPosition.Data$fi
colCity <- c('darkblue', 'red4')
nI <- nrow(Fi)
col4I.City <- rep("",nI)

for (i in 1:length(colCity) ){
  lindex <- cleanData.allVar[,'nation'] %in% unique(cleanData.allVar[,'nation'])[i]
  col4I.City[lindex] <- colCity[i]
}
# generate the set of maps
BaseMap.Fi <- createFactorMap(X = Fi , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              title = 'MCA. Observations (by nation)', 
                              col.points = col4I.City,
                              alpha.points = .4, cex = .9,
                              col.labels = col4I.City,
                              text.cex = 2.5, 
                              force = 2)
# make the maps
d0001.BaseMapNoLabels.Fi  <- BaseMap.Fi$zeMap_background +
                                  BaseMap.Fi$zeMap_dots + labels4MCA 

## ----plotaMapi, fig.width= 8---------------------------------------------
print(d0001.BaseMapNoLabels.Fi)

## ----Boot4CI-------------------------------------------------------------
# Bootstrap for CI:
BootCube.Gr <- PTCA4CATA::Boot4Mean(resMCA$ExPosition.Data$fi, 
                                 design = cleanData.allVar$nation,
                                 niter = 100,
                                 suppressProgressBar = TRUE)
nationsMeans <- PTCA4CATA::getMeans(resMCA$ExPosition.Data$fi, cleanData.allVar$nation)
# colCity <- c('darkblue', 'red4')
MapGroup <- PTCA4CATA::createFactorMap(nationsMeans,
                            # use the constraint from the main map
                            constraints = BaseMap.Fi$constraints,
                            col.points = colCity,
                            cex = 7,  # size of the dot (bigger)
                            col.labels = colCity,
                            text.cex = 6)
d002.Map.I.withMeans <- d0001.BaseMapNoLabels.Fi  +
                          MapGroup$zeMap_dots + MapGroup$zeMap_text
print(d002.Map.I.withMeans)

## ----graphElli-----------------------------------------------------------
GraphElli <- PTCA4CATA::MakeCIEllipses(BootCube.Gr$BootCube[,1:2,],
                            names.of.factors = c("Dimension 1","Dimension 2"),
                            col = colCity,
                            p.level = .95)
d003.Map.I.withCI <-  d0001.BaseMapNoLabels.Fi + 
                          MapGroup$zeMap_text +  GraphElli
print(d003.Map.I.withCI)

## ----TI------------------------------------------------------------------
GraphTI.Hull <- PTCA4CATA::MakeToleranceIntervals(resMCA$ExPosition.Data$fi,
                            design = as.factor(cleanData.allVar$nation),
                            # line below is needed
                            names.of.factors =  c("Dim1","Dim2"), # needed 
                            col = colCity,
                            line.size = .50, 
                            line.type = 3,
                            alpha.ellipse = .2,
                            alpha.line    = .4,
                            p.level       = .75)
#_____________________________________________________________________
# Create the map:
d005.Map.I.withTIHull <- d002.Map.I.withMeans  +
                           GraphTI.Hull + MapGroup$zeMap_dots +
                           MapGroup$zeMap_text + MapGroup$zeMap_dots
#_____________________________________________________________________
# plot it
# dev.new()
print(d005.Map.I.withTIHull)

## ----createFiMap.sup-----------------------------------------------------
Fi.sup  <- resMCA.sup$fii
col     <- 'green'  
nI.sup <- nrow(Fi.sup)
col4I.sup <- rep("",nI.sup)
# generate the set of maps
BaseMap.Fi.sup <- createFactorMap(X = Fi.sup , # resMCA$ExPosition.Data$fj,
                              axis1 = axis1, axis2 = axis2,
                              constraints = BaseMap.Fi$constraints,
                              title = '', 
                              col.points = 'green',
                              alpha.points = .2, cex = 1.2,
                              col.labels = 'green',
                              text.cex = 2.5,
                              force = 2)
# make the maps
e0001.BaseMapNoLabels.Fi.sup  <- BaseMap.Fi$zeMap_background +
                                   BaseMap.Fi.sup$zeMap_dots + 
                                   BaseMap.Fi.sup$zeMap_text + 
                                   BaseMap.Fi$zeMap_dots +
  ggplot2::ggtitle('MCA Active and Supplementary Observations') +
                                   labels4MCA 

## ----plotaMapi.sup, fig.width= 8-----------------------------------------
 print(e0001.BaseMapNoLabels.Fi.sup)

## ----saveGraphs, message = FALSE, warning = FALSE, error = FALSE, eval = FALSE----
#  list2Graphs <- PTCA4CATA::saveGraph2pptx(file2Save.pptx = name4Graphs,
#                   title = 'Attitudes toward fermented products',
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

