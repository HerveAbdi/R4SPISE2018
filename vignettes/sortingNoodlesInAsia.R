## -----------------------------------------------------------------------------
# A clean start
rm(list = ls())
graphics.off()

## ----setup, include = FALSE, ECHO = FALSE-------------------------------------
# Important: Remember 
#     build the vignettes with devtools::build_vignettes()
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width = 9,
  comment = "#>"
)

## ----pdfIt, echo = FALSE, warning = FALSE, eval = FALSE-----------------------
#  # To get the vignette as a pdf file
#  rmarkdown::render("sortingNoodlesInAsia.Rmd",
#           bookdown::pdf_document2(keep_tex = TRUE,
#           toc_depth = 3)  ,
#           output_dir = '../../test4Spise2018/sortingNoodles/')

## ---- eval = FALSE,ECHO = FALSE , include = FALSE-----------------------------
#  # Knitr options here
#  knitr::opts_knit$get()

## ----nsubject, echo = FALSE---------------------------------------------------
nK = 30

## ----zeNoodles, echo=FALSE, fig.cap="The Ramen Noodles Pictures \\label{fig:ramenPictures}", fig.height=4, fig.width=6, include=TRUE, out.width='85%'----
knitr::include_graphics('../man/figures/ramenNoodlesImages.png')

## ---- include = TRUE, echo = TRUE---------------------------------------------
rm(list = ls())
graphics.off()

## ----loadPackages-------------------------------------------------------------
# Decomment all/some these lines if the packages are not installed
# devtools::install_github('HerveAbdi/PTCA4CATA')
# devtools::install_github('HerveAbdi/DistatisR')
#  install.packages(prettyGraphs)
#  install.packages('Matrix')
#  install.packages('factoextra')
#  install.packages('ExPosition')
#  install.packages('pander') # nice pdf-tables
#
#  load the libraries that we will need
suppressMessages(library(Matrix))
suppressMessages(library(factoextra))
suppressMessages(library(DistatisR))
suppressMessages(library(PTCA4CATA))
suppressMessages(library(prettyGraphs))
suppressMessages(library(ExPosition))




## ----filenome-----------------------------------------------------------------
file2read.name <- 'dataSortingRamen.xlsx'
path2file <- system.file("extdata", file2read.name, package = "R4SPISE2018")
sheetName4Data       <- 'DataSort'
sheetName4Vocabulary <- "Vocabulary"
sheetName4Judges     <- "JudgesDescription"

## ----findDataPath-------------------------------------------------------------
path2file <- system.file("extdata",
       "dataSortingRamen.xlsx", package = "R4SPISE2018")

## ----xls.datafile, echo=FALSE, fig.cap="The Data Excel File \\label{fig:spicesxl}", fig.height=3, fig.width=4, include=TRUE, out.width='70%'----
# label.spicesxl = capFigNo
# ![Toto ](Screen Shot 2018-07-09 at 13.51.37.png)
# decomment if needed 
knitr::include_graphics('../man/figures/imageOfSortingNoodles.png')
#knitr::include_graphics('imageOfSortingNoodles.png')
# copy the file 
# file.copy(path2file, 'ramen.xlsx')


## ----resdSortingData----------------------------------------------------------
# read the sortinig data and the vocabulary
multiSort.list <- read.df.excel(path = path2file, 
                    sheet = sheetName4Data,
                    voc.sheet = sheetName4Vocabulary)
multiSort  <- multiSort.list$df.data
vocabulary <- multiSort.list$df.voc

## ----savexls, echo = FALSE, eval = FALSE, include = FALSE---------------------
#  # saveFile <- file.copy(from = path2file, to = '~/Downloads/myDataFile.xlsx')

## ----peekASort----------------------------------------------------------------
library(pander)
# knitr::kable(multiSort[1:5,1:10])
pander::pander(multiSort[1:5,1:10])

## ----readJudges---------------------------------------------------------------
# read the sortinig data and the vocabulary
judgesDescription <- read.df.excel(path = path2file, 
                                  sheet = sheetName4Judges)$df.data
nVarJudges  <- ncol(judgesDescription)

## ----k4Judges-----------------------------------------------------------------
k <- 1 # this is the first descriptor for the judges

## ----descJudges---------------------------------------------------------------
descJudges <- judgesDescription[,k ]

## ----colJudges----------------------------------------------------------------
# Create a 0/1 group matrix with ExPosition::makeNominalData()
nominal.Judges <- makeNominalData(as.data.frame(descJudges))
# get the colors
color4Judges.list <- prettyGraphs::createColorVectorsByDesign(nominal.Judges)
# color4Judges.list

## ----getCube------------------------------------------------------------------
DistanceCube <- DistanceFromSort(multiSort)

## ----runDistatis--------------------------------------------------------------
resDistatis <- distatis(DistanceCube)

## ----rvGroups-----------------------------------------------------------------
# Get the factors from the Cmat analysis
G <- resDistatis$res4Cmat$G 
# Compute the mean by groups of HJudges
JudgesMeans.tmp <- aggregate(G, list(descJudges), mean) 
JudgesMeans <- JudgesMeans.tmp[,2:ncol(JudgesMeans.tmp )] 
rownames(JudgesMeans) <- JudgesMeans.tmp[,1]
# Get the bootstrap estimates
BootCube <- PTCA4CATA::Boot4Mean(G, design = descJudges,
                       niter = 100,
                       suppressProgressBar = TRUE)
# head(BootCube)

## ----computeSk----------------------------------------------------------------
F_j     <- resDistatis$res4Splus$PartialF
alpha_j <- resDistatis$res4Cmat$alpha
# create the groups of Judges
#groupsOfJudges <- substr(names(alpha_j),1,1)
groupsOfJudges <- descJudges
code4Groups <- unique(groupsOfJudges)
nK <- length(code4Groups)
# initialize F_K and alpha_k
F_k <- array(0, dim = c(dim(F_j)[[1]], dim(F_j)[[2]],nK))
dimnames(F_k) <- list(dimnames(F_j)[[1]], 
                         dimnames(F_j)[[2]], code4Groups)
alpha_k <- rep(0, nK)
names(alpha_k) <- code4Groups
Fa_j <- F_j
# A horrible loop
for (j in 1:dim(F_j)[[3]]){ Fa_j[,,j]  <- F_j[,,j] * alpha_j[j] }
# Another horrible loop
for (k in 1:nK){
  lindex <- groupsOfJudges == code4Groups[k]
  alpha_k[k] <- sum(alpha_j[lindex])
  F_k[,,k] <- (1/alpha_k[k])*apply(Fa_j[,,lindex],c(1,2),sum)
}


## ----projVoc------------------------------------------------------------------
F4Voc <- projectVoc(multiSort.list$df.voc, resDistatis$res4Splus$F)

## ----RV.scree.MapPlain, fig.height=4, fig.width= 7----------------------------
# 5.A. A scree plot for the RV coef. Using standard plot (PTCA4CATA)
scree.rv.out <- PlotScree(ev = resDistatis$res4Cmat$eigValues, 
                   title = "RV-map: Explained Variance per Dimension")
a1.Scree.RV <- recordPlot() # Save the plot

## ----RVGplot------------------------------------------------------------------
# Create the layers of the map
gg.rv.graph.out <- createFactorMap(X = resDistatis$res4Cmat$G, 
                            axis1 = 1, axis2 = 2, 
                            title = "Judges: RVMap", 
                            col.points = color4Judges.list$oc, 
                            col.labels = color4Judges.list$oc)
# create the labels for the dimensions of the RV map
labels4RV <- createxyLabels.gen(
                  lambda = resDistatis$res4Cmat$eigValues , 
                  tau    = resDistatis$res4Cmat$tau,
                  axisName = "Dimension ")
# # Create the map from the layers
# Here with lables and dots
a2a.gg.RVmap <- gg.rv.graph.out$zeMap + labels4RV
# Here with colored dots only
a2b.gg.RVmap <- gg.rv.graph.out$zeMap_background +
                gg.rv.graph.out$zeMap_dots + labels4RV

## ----mapa2a, fig.height=6, fig.width= 9---------------------------------------
print(a2a.gg.RVmap )

## ----RVwithCI-----------------------------------------------------------------
# First the means
# A tweak for colors
in.tmp    <- sort(rownames(color4Judges.list$gc), index.return = TRUE)$ix
col4Group <- color4Judges.list$gc[in.tmp]
#
gg.rv.means <- PTCA4CATA::createFactorMap(JudgesMeans,
                      axis1 = 1, axis2 = 2, 
                      constraints = gg.rv.graph.out$constraints,
                      col.points =  col4Group ,
                      alpha.points = 1, # no transparency
                      col.labels = col4Group)
#
 dimnames(BootCube$BootCube)[[2]] <- 
                    paste0('dim ',1: dim(BootCube$BootCube)[[2]])
  #c('Dim1','Dim2') 
GraphElli.rv <- MakeCIEllipses(BootCube$BootCube[,1:2,],
                 names.of.factors = c("dim 1","dim 2"), 
                 col = col4Group, 
                 p.level = .95)
a2d.gg.RVMap.CI <- a2b.gg.RVmap + gg.rv.means$zeMap_dots + GraphElli.rv 

## ----meansRV------------------------------------------------------------------
knitr::kable(JudgesMeans[,1:3])


## ----mapa2d, fig.height=6, fig.width= 9---------------------------------------
print(a2d.gg.RVMap.CI )

## ----HCA----------------------------------------------------------------------
 D <- dist(resDistatis$res4Cmat$G, method = "euclidean")
 fit <- hclust(D, method = "ward.D2")
 a05.tree4participants <- fviz_dend(fit,  k = 1, 
                        k_colors = 'burlywood4', 
                        label_cols = color4Judges.list$oc[fit$order],
                        cex = .7, xlab = 'Participants',
                        main = 'Cluster Analysis: Participants') 

## ----plothca, fig.height = 9, fig.width = 9-----------------------------------
 print(a05.tree4participants)

## ----kmeans-------------------------------------------------------------------
# First plain k-means
set.seed(42)
participants.kMeans <- kmeans(x = G , centers = 4)
#_____________________________________________________________________
# Now to get a map by cluster:
col4Clusters  <- createColorVectorsByDesign(
              makeNominalData(
              as.data.frame(participants.kMeans$cluster)  ))

## ----kmeans.graph, fig.height=6, fig.width= 9---------------------------------
baseMap.i.km <- PTCA4CATA::createFactorMap(G,
                             title = "RV map. k-means 4 groups",        
                                  col.points = col4Clusters$oc,
                                  col.labels = col4Clusters$oc,
                 constraints =    gg.rv.graph.out$constraints,
                                  alpha.points =  .4)
a06.aggMap.i.km <- baseMap.i.km$zeMap_background +
  baseMap.i.km$zeMap_dots + baseMap.i.km$zeMap_text +  labels4RV
print(a06.aggMap.i.km)

## ----scree4S, fig.height=4, fig.width=7---------------------------------------
#---------------------------------------------------------------------
# A scree plot for the Compromise.
scree.S.out <- PlotScree(
              ev = resDistatis$res4Splus$eigValues, 
              title = "Compromise: Explained Variance per Dimension")
b1.Scree.S <- recordPlot()
#---------------------------------------------------------------------

## ----createGr4S, echo=TRUE, error=FALSE, warning=FALSE,message=FALSE,results=FALSE----
# 4.1 Get the bootstrap factor scores (with default 1000 iterations)
BootF <- BootFactorScores(resDistatis$res4Splus$PartialF)
# 5.2 a compromise plot
# General title for the compromise factor plots:
genTitle4Compromise = 'Compromise.'
# To get graphs with axes 1 and 2:
h_axis = 1
v_axis = 2
# To get graphs with say 2 and 3 
# change the values of v_axis and h_axis
color4Products <- #  Create color for the Products from prettyGraph
 prettyGraphsColorSelection(n.colors = nrow(resDistatis$res4Splus$F))
gg.compromise.graph.out <- createFactorMap(resDistatis$res4Splus$F,
                                    axis1 = h_axis, 
                                    axis2 = v_axis,
                                    title = genTitle4Compromise,
                                    col.points = color4Products ,
                                    col.labels = color4Products) 
# NB for the lines below You need DISTATIS version > 1.0.0
#  to get the eigen values and tau for the compromise
label4S <- createxyLabels.gen(
            x_axis   = h_axis, y_axis = v_axis,
            lambda   = resDistatis$res4Splus$eigValues , 
            tau      = resDistatis$res4Splus$tau,
            axisName = "Dimension ")
b2.gg.Smap <-  gg.compromise.graph.out$zeMap + label4S 
#  
# 5.4 a bootstrap confidence interval plot 
# 5.3  create the ellipses
gg.boot.graph.out.elli <- MakeCIEllipses(
                              data = BootF[,c(h_axis,v_axis),],
                              names.of.factors = 
                                c(paste0('Factor ',h_axis),
                                  paste0('Factor ',v_axis)),
                              col = color4Products,
)  
# Add ellipses to compromise graph
b3.gg.map.elli <- gg.compromise.graph.out$zeMap + gg.boot.graph.out.elli + label4S 
#

## ----plot4S, fig.height=6, fig.width= 9---------------------------------------
print(b2.gg.Smap)

## ----cluster4Prod-------------------------------------------------------------
nFac4Prod = 3
D4Prod <- dist(resDistatis$res4Splus$F[,1:nFac4Prod], method = "euclidean")
 fit4Prod <- hclust(D4Prod, method = "ward.D2")
 b3.tree4Product <- fviz_dend(fit4Prod,  k = 1, 
                        k_colors = 'burlywood4', 
                        label_cols = color4Products[fit4Prod$order],
                        cex = .7, xlab = 'Products',
                        main = 'Cluster Analysis: Products') 

## ----plothcaProd, fig.height = 9, fig.width = 9-------------------------------
 print(b3.tree4Product)

## ----PartialFS----------------------------------------------------------------
# get the partial map
map4PFS <- createPartialFactorScoresMap(
                         factorScores = resDistatis$res4Splus$F,      
                          partialFactorScores = F_k,  
                          axis1 = 1, axis2 = 2,
                          colors4Items = as.vector(color4Products), 
                          names4Partial = dimnames(F_k)[[3]], # 
                          font.labels = 'bold')

d1.partialFS.map.byProducts <- gg.compromise.graph.out$zeMap + 
                                  map4PFS$mapColByItems + label4S 
d2.partialFS.map.byCategories  <- gg.compromise.graph.out$zeMap + 
                                  map4PFS$mapColByBlocks + label4S 

## ----SwithCategories.1, fig.height=6, fig.width= 9, message = FALSE, warning = FALSE, error = FALSE----
print(d1.partialFS.map.byProducts )

## ----SwithCategories.2, fig.height=6, fig.width= 9, message = FALSE, warning = FALSE, error = FALSE----
print(d2.partialFS.map.byCategories)

## ----graphVoc-----------------------------------------------------------------
# 5.5. Vocabulary
# 5.5.2 CA-like Barycentric (same Inertia as products)
gg.voc.bary <- createFactorMap(F4Voc$Fvoca.bary,
                    title = 'Vocabulary',
                    col.points = 'red4',
                    col.labels = 'red4',
                    display.points = FALSE,
                    constraints = gg.compromise.graph.out$constraints)
#
e1.gg.voc.bary.gr <- gg.voc.bary$zeMap + label4S 

#print(e1.gg.voc.bary.gr)
b5.gg.voc.bary.dots.gr <- gg.compromise.graph.out$zeMap_background +
                          gg.compromise.graph.out$zeMap_dots + 
                          gg.voc.bary$zeMap_text + label4S 
#print(gg.voc.bary.dots.gr)

## ----vocbary, fig.height=6, fig.width = 9-------------------------------------
print(e1.gg.voc.bary.gr)

## ----vocbaryProd, fig.height=6, fig.width = 9---------------------------------
print(b5.gg.voc.bary.dots.gr)

## ----saveGraphs, message = FALSE, warning = FALSE, error = FALSE, eval = FALSE----
#  toto <- PTCA4CATA::saveGraph2pptx(file2Save.pptx = name4Graphs,
#                   title = '30 (mostly Asian) participants sort pictures of 20 Ramen noodles',
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

