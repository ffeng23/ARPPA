
#' @title an empty file to import necessary libraries 
#' @import limma MASS
#' @name _empty
NULL


#'@title S3 function to read/import the data from a text file
#'@description \code{importTextData} reads/imports the data from a text file
#'@details This function take a text format file as input.
#'	 The file has a fixed format. It assumes the input file is
#'	generated by the protoarray software (prospector).
#'	Please see the sample data.
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
importTextData<-function(dataFilePath=NULL, targetFile=NULL, start.data=NULL, nrows.data=NULL,
				start.control=NULL, nrows.control=NULL,aggregation="geoMean",
				as.is=TRUE, header=TRUE,sep="\t",na.strings="", quote="")
		{
			 if (is.null(dataFilePath) || is.null(targetFile)) {
				stop("ERROR: Not all mandatory arguments have been defined!")
				}
			if (start.data<=0|| nrows.data<=0||start.control<=0||nrows.control<=0) {
				stop("ERROR: Not all mandatory arguments have been defined!")
				}
			if(aggregation!="geoMean"&&aggregation!="arithMean"&&aggregation!="min")
			{
				cat("Warning:Unknow aggregation type, use default (min)!\n")
				aggregation="min";
			}
			#now read the data first
			cat("reading target files for \"", targetFile, "\"...\n")
			flush.console();
			targetRd<-read.table(targetFile, sep="\t", header=TRUE, as.is=TRUE, na.strings="", quote="");
			cat("Done!\n")
			original_path<-getwd()
			setwd(dataFilePath);#,
			fileNames<-targetRd[,"FileName"]
			
			for(i in c(1:length(fileNames)))
			{
				cat("reading data files for ", fileNames[i], "...\n")
				flush.console();
				dataRd<-read.table(fileNames[i], sep="\t", header=TRUE, skip=51, nrows=18803-53,as.is=TRUE, na.strings="", quote="");
				ctrDataRd<-read.table(fileNames[i], sep="\t", header=TRUE, skip=18803, nrows=-1,as.is=TRUE, na.strings="", quote="");

				
				#put them together
				if(i==1)
				{
					#control data
					ctrData_signal<-data.frame(ctrDataRd[,'Signal']);
					ctrData_background<-data.frame(ctrDataRd[,'Background']);
					#sample data
					dataAll_signal<-data.frame(dataRd[,'Signal']);
					dataAll_background<-data.frame(dataRd[,'Background']);
					
				}
				else
				{
					#control data
					ctrData_signal<-cbind(ctrData_signal,ctrDataRd[,'Signal']);
					ctrData_background<-cbind(ctrData_background,ctrDataRd[,'Background']);
					
					#sample data
					dataAll_signal<-cbind(dataAll_signal, dataRd[,'Signal']);
					dataAll_background<-cbind(dataAll_background, dataRd[,'Background']);
					
				}				
			}#end of read data
			cat("Done!\n")
			sampleNames<-targetRd[,"ArrayID"]
			colnames(dataAll_signal)<- sampleNames;
			colnames(dataAll_background)<- sampleNames;
			colnames(ctrData_signal)<- sampleNames;
			colnames(ctrData_background)<- sampleNames;

			##aggregation,summary/average the duplicates
			sum_index<-seq(1,length(dataRd[,1]),2);
			#dataAll_signal<-c();
			#dataAll_background<-c();
			csum_index<-seq(1,length(ctrDataRd[,1]),1);
			if(aggregation=="geoMean")
			{
				#cat("geonMean\n");
				dataAll_signal<-exp((log(dataAll_signal[sum_index,])+log(dataAll_signal[sum_index+1,]))/2)
				dataAll_background<-exp((log(dataAll_background[sum_index,])+log(dataAll_background[sum_index+1,]))/2)
				#cat("length dataAll:", dim(dataAll_signal
			}
			else if(aggregation=="arithMean")
			{
				cat("arithMean\n");
				dataAll_signal<-(dataAll_signal[sum_index,]+dataAll_signal[sum_index+1,])/2
				dataAll_background<-(dataAll_background[sum_index,]+dataAll_background[sum_index+1,])/2
			}
			else #default min
			{
				cat("min\n");
				dataAll_temp1<-(dataAll_signal[sum_index,])
				dataAll_temp2<-(dataAll_signal[sum_index+1,])
				dataAll_temp1[dataAll_temp1>dataAll_temp2]<-dataAll_temp2[dataAll_temp1>dataAll_temp2]
				dataAll_signal<-dataAll_temp1
				
				dataAll_temp1<-(dataAll_background[sum_index,])
				dataAll_temp2<-(dataAll_background[sum_index+1,])
				dataAll_temp1[dataAll_temp1>dataAll_temp2]<-dataAll_temp2[dataAll_temp1>dataAll_temp2]
				dataAll_background<-dataAll_temp1
			}
			
			#now use the last readin file to create gene info frame
			genes<-dataRd[sum_index,c('Block','Column','Row','Protein.Amount','Description')];

			genes$Name<-paste("Hs~",dataRd[sum_index,'Database.ID'],"~uORF:",dataRd[sum_index,'Ultimate.ORF.ID'],sum_index,sep="");
			genes$ID<-paste("HA20251~",dataRd[sum_index,'Array.ID'],sep="");


			cgenes<-ctrDataRd[csum_index,c('Block','Column','Row')];
			cgenes$Description<-rep('Control',length(csum_index))
			cgenes$Name<-ctrDataRd[csum_index,'ControlGroup'];
			cgenes$ID<-paste("HA20251~",ctrDataRd[csum_index,'ID'],sep="");
			#cgenes$
			dataAll_signal<-as.matrix(dataAll_signal)
			dataAll_background<-as.matrix(dataAll_background)
			ctrData_signal<-as.matrix(ctrData_signal);
			ctrData_background<-as.matrix(ctrData_background);

			#create a target object dataframe
			#targets<-data.frame(ArrayID=sampleNames,FileName=fileNames,Group=group, Batch=rep('B1',length(sampleNames)),
			#					Data=rep('11.21.2015',length(sampleNames)),Array=c(1:length(sampleNames)),SerumID=sampleNames);
			#this is part of information for organization of section in the array
			printer<-list(ngrid.r=12, ngrid.c=4,nspot.r=22, nspot.r=22);					
								
			el<-list(E=dataAll_signal, Eb=dataAll_background, targets=targetRd, genes=genes,
					source="genepix.median", printer=printer, C=ctrData_signal, Cb=ctrData_background, cgenes=cgenes);
			setwd(original_path);		
			elist2<-new("EListRaw", el)
			##........now we have object, do the jobs
		}#end of function###########
		
##############testing area###########
