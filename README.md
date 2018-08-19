# R4SPISE2018 .0.1.0

## R4SPISE2018: a set of data sets, scripts, and vignettes used for the Advanced Workshop of the 2018 SPISE meeting. Da Nang, Vietnam: July 20 to July 30, 2018. 

## Theme of the workshop: *individual differences* in sensory evaluation

## Packages to install prioir to installingn this package

Prior to installing the `R4SPISE2018` package, we need to install some other pacakges (msotly from `Github`).

```{r}
# install.packages('devtools') # decomment this line if devtools is not yet installed
install.packages('factoextra')
devtools::install_github('HerveAbdi/PTCA4CATA')
devtools::install_github('HerveAbdi/DistatisR')
```

## How to install the package

To install the packages, use the following `R`-command:
```{r}
# install.packages('devtools') # decomment this line if devtools is not yet installed
devtools::install_github('HerveAbdi/R4SPISE2018', 
      dependencies = TRUE, # use it first time only comment after
      build_vignettes = TRUE) # to get the vignettes
```



## Main statistical techniques  and R-packages used

1. Distatis: package `DistatisR`
2. Partial Triadic Correspondence Analysis: package `PTCA4CATA`
3.  Multiple Correspondence Analysis: package `ExPosition`
