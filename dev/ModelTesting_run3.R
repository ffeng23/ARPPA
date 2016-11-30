##########################################
##the code here is used to run simulation tests for 
##checking the model and methods. 
##1.testing whether the simulated data can be used to
##	detecting the nonZero interaction effects. ----it should
##2.checking for the heterscedasticity cases vs. homoscedasticity
##3.checking for validity of the modeling without control data
##4.all above first using the smaller number of genes---- to keep 
##	the memory usage low so as to using the R built-in linear 
##	function for testing
##
##	Oct 1st, 2016 Feng@BU
##	see ModelTesting_run1.R for equal variance cases
##
##@Oct 23, 2016 by Feng
##ModelTesting code file #3
##this one is used to read the file and run analysis for unequal variance cases
##here, we need to do the transformation!!!
##since clearly the unequal variances have lower true discovery rate.
## To transform we need to calculate the expected variance using bayesian approach
## (smyth 2004), then transform to stabalize the variance by 
##	 var_hat=(d0S0^2+diSi^2)/(d0+di)
##		Y_t=(Y-Y_bar)/sqrt(var_hat)+Y_bar
##by doing this, we kind of stablized the variance, but did not 
##changed the variance
##
##########################################
library(ARPPA)
##first checking for homoscedastical dataset
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon")
nGene<-2000
nTreatment<-2   #number of different beta
sampleSize<-5 ##this is the number of repeats for each group

#================analysis=======
###start doing the summary statistics of the repeated runs
#
#loading the data
repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_NegativeCon");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
i<-1
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave
	rm(lstToSave)
	
	####data transformation###
	##First, estimate the variance prior, bayesian
	dataExp<-lstData[[1]]
	dataExp<-dataExp[[1]]
	svList<-sampleVariance(dataExp,nTreatment, sampleSize)
	#reformat it into a vector
	s2<-as.vector(as.matrix(svList))
	df<-as.vector(sampleSizes(dataExp,nTreatment, sampleSize))
	r2<-sampleVariancePrior(s2,df)
	s0<-r2[[2]]
	df0<-r2[[1]]

	##
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex<-which(gammaTrue[,2]!=0)
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j])
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="xy")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.1,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]

########=================Isotype control data analysis
#loading the data
repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_IsotypeCon");

#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave_EI
	rm(lstToSave_EI)
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex<-which(gammaTrue[,2]!=0)
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j])
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="xy")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.02,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]
 
 
 #######====================No control++++++++++

#loading the data
repeats<-100
setwd("H:\\feng\\LAB\\hg\\proteinArray_Masa\\ARPPA\\data\\UnEqualVar_UnControl");
#setwd("~/Desktop/arppr/EqualVar_NegativeCon/")
i<-1
portions<-seq(0.001, 0.5, 0.001)
stats_df<-data.frame("portions"=portions);
for(i in c(1:repeats))
{
	cat("reading the ",i, "/", repeats," data sets.....\n")
	load(paste("result_", i,".RData",sep=""))
	cat("count the correct ones.....\n")
	#depending on which data,
	lstData<-lstToSave_ENC
	rm(lstToSave_ENC)
	
	#try to know which one is nonZero as the true values
	gammaTrue<-lstData$data$gamma
	gammaTrueIndex<-which(gammaTrue[,2]!=0)
	
	#now need to sort the p-value of the 
	#and pick top n genes to be the significant ones as the analysis results
	#
	lr_coe<-lstData$lregCoef
	lr_coe_names<-rownames(lr_coe)
	lr_coe_gamma_index<-which(grepl("gene\\d*:group2",lr_coe_names));
	lr_coe_int<-lr_coe[lr_coe_gamma_index,]
	lr_coe_int_names<-lr_coe_names[lr_coe_gamma_index];

	#sort the array according to the prob
	lr_coe_int_sort_order<-order(lr_coe_int[,4]) #NOTE:the names is one more than its index, eg. 811 is gene812:group2
	#lr_coe_int_sort<-lr_coe_int[lr_coe_int_sort_order,]
	

	#now we have everything, just need to collect statistics
	p_correct<-portions;#initialize the vector
	for(j in c(1:length(portions)))
	{
		#for each portion, we need to check what is percentage to be correct
		numToPick<-floor(lstData$data$params[1]*portions[j])
		#picking from the order array
		genesToPick<-lr_coe_int_sort_order[c(1:numToPick)]+1

		#now we got the genes, just need to check whether the genes so far are the TRUE ones
		numOfCorrect<-sum(is.element(genesToPick, gammaTrueIndex));
		p_correct[j]<-numOfCorrect/numToPick
	}
	stats_df[,paste("repeat_",i)]<-p_correct;
}
#statistics with the stats_df array
mean_correct<-portions
max_correct<-portions
min_correct<-portions
std_correct<-portions
m_stats_df<-as.matrix(stats_df)
for(k in c(1:length(portions)))
{
	mean_correct[k]<-mean(m_stats_df[k,c(-1)])
	max_correct[k]<-max(m_stats_df[k,c(-1)])
	min_correct[k]<-min(m_stats_df[k,c(-1)])
	std_correct[k]<-sqrt(var(m_stats_df[k,c(-1)]))
}

stats_df[,"mean"]<-mean_correct
stats_df[,"min"]<-min_correct
stats_df[,"max"]<-max_correct
stats_df[,"std"]<-std_correct

plot(c(0.001,0.51),c(0.02, 1.1), type="n", main="true discovery rate for protein selection", 
	xlab="portion of protein selected", ylab="portion of true discovery", log="xy")
lines(stats_df[,1], stats_df[,"mean"], col=1, lty=1, lwd=2)
lines(stats_df[,1], stats_df[,"min"], col="grey", lty=2, lwd=1)
lines(stats_df[,1], stats_df[,"max"], col="grey", lty=2, lwd=1)
lines(c(0.01, 0.01), c(0.02,1.2), col=2, lty=3)

#show stats
stats_df[c(1:30),c("portions", "mean", "min", "max", "std")]
 