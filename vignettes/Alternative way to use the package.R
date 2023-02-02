#### Author: Xiaodong Feng. 2023 ####
# This script provides alternative way to use the R package DynamicTol
# In case the installation of DynamicTol is failed, you can load all the packages 
# and functions at once by using the following scripts
#### Install and load libraries ####
# .libPaths()
library(blandr)
library(CAMERA)
library(cowplot)
library(ChemmineOB)
library(ChemmineR)
# library(DBI)
library(data.table)
library(dplyr)
library(dendextend) #color_labels
library(fitdistrplus)
library(fmcsR)
library(ggpubr)
library(ggplot2)
# library(grid)
# library(import)
# library(mzR)
library(MSnbase)
# library(multtest)
# library(magrittr)
# library(metfRag)
# library(msdata)
library(msPurity)
# library(msPurityData)
# library(pander)
library(plyr)
library(pROC) # to plot the receiver operating characteristic (ROC)
# library(RSQLite)
# library(rJava)
library(Rcpp)
library(readxl)
library(reshape2)
library(reticulate) # used to load the python functions
# library(RColorBrewer)
library(snow)
# library(styler)
# library(stringr)
library(Spectra)
library(writexl)
library(xlsx)
library(MSnbase)
library(xcms)
#### Load functions ####
# The functions are separated in the General.R, Plot.R and SpectralMatching.R
# We just need to load these three files via the source function.
setwd('d:/github/DynamicTol')
source('R/General.R')
source('R/Plot.R')
source('R/SpectralMatching.R')
