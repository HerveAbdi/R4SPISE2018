---
title: "Package: R4SPISE2018"
author: "Hervé Abdi"
date: "`October 23, 2018`"
---

# R4SPISE2018 .0.1.0

 A set of data sets, scripts, and vignettes 
 used for the *Advanced Workshop of the 2018 SPISE Meeting*
 that took place in Da Nang, Vietnam: July 20 to July 30, 2018. 

## Workshop: *individual differences* in sensory evaluation

## Packages to install prior to installing this package

Prior to installing the `R4SPISE2018` package, 
we need to install some other packages (mostly from `Github`).

```{r}
# decomment the line below if devtools is not yet installed
# install.packages('devtools') 
install.packages('factoextra')
devtools::install_github('HerveAbdi/PTCA4CATA')
devtools::install_github('HerveAbdi/DistatisR')
devtools::install_github('HerveAbdi/data4PCCAR')
```

## How to install the package

To install the package `R4SPISE2018`, use the following `R`-command:
```{r}
# install.packages('devtools') # decomment this line 
if devtools is not yet installed
devtools::install_github('HerveAbdi/R4SPISE2018', 
       dependencies = TRUE, # use it first time only, comment after
       build_vignettes = TRUE) # to get the vignettes
```
The installation of the "dependent" packages can take some time and
building the vignettes can take also a bit of time to build:
So, be patient!

## How to build the vignettes

If you have forgotten to build the vignettes, you can build them with
the command:
```{r}
devtools::build_vignettes()
```
Again, be prepared to wait till the vignettes are built.

## Main statistical techniques and R-packages used

1. Distatis: package `DistatisR`
2. Partial Triadic Correspondence Analysis (PTCA): package `PTCA4CATA`
3. Multiple Correspondence Analysis (MCA): package `ExPosition`

### Where to download the packages

    * `ExPosition` can be downloaded from `CRAN`
    * `DistatisR` and `PTCA4CATA` should be downloaded from `HerveAbdi/Github'

## Current Vignettes


### beersCATA		

*Correspondence Analysis for CATA:*

Novices and experts evaluated 9 beers.

### cheeseMCA	

*Multiple Correspondence Analysis:*

Analyzing a survey about cheese

### fermentationIn2NationsMCA

*MCA with supplementary observations:*

Analysis of a survey about attitudes towards fermented food answered by 373 participants (220 French and 183 Vietnamese). In addition, 30 participants from the SPISE2018 advanced workshop answered the questionnaire and are considered as supplementary observations.

This example also shows how to create


### multiculturalSortingSpices	

*DISTATIS:*

Analyzing a sorting task with different groups of assessors.

### sortingNoodlesInAsia

*DISTATIS (with vocabulary and barycentric projections):*

Analyzing a sorting task of pictures of Ramen Noodles
with verbal description of the groups made by  the participants.

