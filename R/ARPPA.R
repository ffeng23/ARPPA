##################################################
# A R Software Package for Protein Array analysis
#    -1.use to read/import the data
#	 -2.background correct the raw data
#    -3.normalize between arrays, by a linear model on log-transformed data
#    -4.statistically identify the differentially bound proteins between arrays
#    	corrected for multiple comparison
#	 -5.test the shifted distribution/pattern of Ab binding in comparison with the control array
#            Kolmogorov-Smirnov Test or Chi-square test
#
#      developed by Feng Feng @ Boston University. 
#      All right reserved
#          2/16/2016
####################################################

##the following is the code to install necessary libs manually
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLitr("PAA")
#####################

##########Dependency###########
##  library(limma)
####################

#' @include fileIO.R
# # ' @include data.R

#'
#' @title APRRA allows to preprocess and analyze protoarray data
#'
#' @description An R package to read, preprocess and analyze protoarray data
#'
#' @details The functions you're likely to need from \pkg{ARPPA} are
#' \code{\link{importTextData}},...to be added???. 
#' Otherwise refer to the vignettes to see
#' how to format the documentation.
"_PACKAGE"

