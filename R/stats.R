#####Description for this file runStats.R#####
# this file is used to define the functions to run statistical analysis 
# in order to identify the proteins differentially binding to the testing Abs
############### 

# the following section is used to act like a holder for
# importing libraryies
# package limma is necessary for EList
#' @title an empty file to import necessary libraries 
#' @import limma 
#' @name _empty
NULL

#first define the top most function to drive the analysis
#the overall work flow
#'@title S3 function to identify features with significant interaction
#'@description \code{identifyFeatures} to identify features/proteins that interact 
#'	siginificantly with the testing Abs
#'@details This function takes in an elist object holding the correctly
#'	 normalized protoarray data, and then run statistical analysis 
#' 	 using an linear regression model with interaction terms. The results 
#'	 are also corrected for multiple analysis using FDR
#'	(see details about the model and analysis \code{\link{..}})
#'	Please see the example below for function usage.
#'
#'@param dataFilePath string path to the file containing the text formatted input data.
#'				This is a tab-delimited file exported from the software reading
#'               the arrays. It has both gene data and control data in one.
#'@param TargetFile string path to a text file holding the meta data for the array. 
#'				Please check the template
#'@param start.data numeric the line number indicating where the data section starts
#'@param nrows.data numeric number of rows for the data section
#'@param start.control numeric the line number indicating where the control data section starts
#'@param nrow.control numeric the number of rows for the control section
#'@param aggregation string indicating which type of duplicate aggregation should be performed. 
#			If "min" is chosen, the value for the corresponding feature will be the minimum of both. 
#			If "arithMean" is chosen, the arithmetic mean will be computed. 
#			If "genMean" is chosen, the geometric mean will be computed.
#			The default is "min" (optional).
#'@param as.is boolean TRUE by default,for reading the text file
#'@param header boolean TRUE by default, for reading the text file
#'@param sep char "\t" by default, for reading the text file
#'@param na.strings string "" by default, for reading the text file
#'@param quote char "" (no quote) by default,  for reading the text file
#'@return an ELISTlist object containing both the data and control
#@seealso
#'@examples
#'	datapath<-system.file("extdata", package="ARPPA")
#'	targets <- list.files(system.file("extdata", package="ARPPA"),
#'		 pattern = "targets_text_Batch1", full.names=TRUE) 

#'	elist2<-importTextData(dataFilePath=datapath, targetFile=targets, start.data=51, nrows.data=18803-53,
#'				start.control=18803, nrows.control=23286,aggregation="geoMean",
#'				as.is=TRUE, header=TRUE,sep="\t",na.strings="", quote="")
#'
#'@export
identifyFeatures<-function(elist, expressionCutOff=50 )
		{
		}