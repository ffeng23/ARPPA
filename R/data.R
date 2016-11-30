###this is the module to define the data simulation and related functions for expression
#####this is r developing code for running simulation to generate 
##the DNA/protein microarray data
## the model (ref: smyth 2004, Bayesian hierarchical model and linear model)
##  linear model:
##             Y_ijk=alpha_i+beta_j+ gamma_ij+ epsilon_ijk
##		see the doc named: "Microarray data simulation.doc"
##    feng@BU   09/03/2016
####################################

#note: all the effects are fixed, although the observed ones are "random", since there is 
# 	 observation error or disturbance. the error or disturbance made the estimated effects 
#    "random". 
##S3 function

#comment for Oxygen2 helper page
#'@title S3 function to simulate expression data
#'@description \code{simulateExpression} simulates the expression data assuming
#'	a specific/non-specific effect model (biologically) and a Bayesian 
#'  hierarchical model (statistically)
#'@details This function is used to simulate the data assuming statistically a 
#'  a linear model with Bayesian and hierarchical structure according to Smyth 2004.
#'  Biologically, it parses the observed expression data into effects from different 
#'	sources, mainly nonspecific and specific ones. The nonspecific effects comes 
#'  due to the nonspecifi binding by either proteins or treatment agents. The specific one 
#'  is cause by the specific interaction due to protein and treatment agents (here antibodies).
#'  Linear model:  >>>>>>>>>To be continue>>>>>>>>>>>>>>
#'	About the control groups in the simulation, it is always true that there are 
#'	gene-wise control groups. Under the real condition, these are the genes that are
#'	printed at very little amount on the array. In case of the controls in the treatment group,
#'	there are 3 cases. First is the isotype control group, in which the isotype 
#'	antibody shows no reactivity to genes, but only nonspecific effects. Second 
#'	one is the negative control, meaning there is no antibody at all but only solution.
#'	in this case, there is no nonspecific or interaction effects. The third case
#'	has no control. "Therefore, we need to compare the every treatment group with
#'	the grand mean to estimate the interaction. The assumption here is that 
#'	non-zero interaction is a rare event and grand mean should be a good estimation
#'	of the nonspecific effects of the gene."
#'@param nGene numeric the number of genes to be simulated
#'@param nTreatment number of treatment groups, eg number of antibodies
#'@param sampleSize number of replicates for each gene by treatment group
#'@param control.isotype logic indicating whether for the treatment group
#'	there is isotype controls.
#'@param control.negative logic indicating whether for the treatment group
#'	there is negative controls, meaning the vehicle control and no antibody.
#'	control.negative and control.isotype are exclusive. if both control.negative
#'	and control.isotype are set to be true, then the negative case is assumed.
#'	When alpha and beta are set, all the control arguments are ignored.
#'	Again, there are always gene control groups and we always set first one
#'	to be the control group by gene.
#'
#'@param alpha numeric vector specifying the non-specific gene effects with 
#'  length of nGene
#'@param alpha.mean and alpha.sigma numeric the parameters used to randomly specify
#	the non-specific alpha effects following the normal distribution with
#	mean and standard deviation. When parameter alpha is specified by
#	the user. These two parameters are ignored
#'@param beta numeric vector specifying the non-specific treatment/antibody
#'  effects with length of nTreatment
#'@param beta.mean and beta.sigma numeric the parameters used to randomly specify
#	the non-specific beta effects following the normal distribution with
#	mean and standard deviation. When parameter alpha is specified by
#	the user. These two parameters are ignored
#'@param gamma numeric The specific effect between gene and treatment.
#'	Under the null hypothesis, this is zero assuming no interacting effects.
#' 	The observed expressions for many of the protein and treatment are
#'	caused mainly because of non-specific effects with disturbance. 
#'@param prob.nonZero numeric this is true partion for those proteins
#' 	having the non-specific interaction with the treatment agents.
#'@param gamma.sigma numeric the parameter specific the distribution from which
#'	the non-zero gamma drawing values. Here we assuming a Gaussian 
#'	distribution with mean of zero and standard deviation of 
#'	gamm.sigma. If gamma is specified by the user, then gamma.sigma and 
#'	prob.nonZero will be ignored.
#'@param epsilon.si numeric matrix the standard deviations for Gaussian
#'	distributed errors/disturbances.	
#'@param epsilon.d0 and epsilon.s0 numeric The hyperparameter specifying
#'	the prior distribution of for the variance for each inidividual group.
#'	Here we assume a heteroscedastical error model. 
#'		episilon ~ N(0, si2)
#'		d0*s02/si2~chiSquare(d0,di)
#'		Yijk=Yijk_hat+epsilon 
#'		di is the df of si2
#'
#'@return a list containing 1)a numeric data matrix with a format of gene by treatment/antibody. 
#'	The columns are arrangend to have repeats; 2)parameters alpha for genes;
#'	3)parameters beta for treatment; 4)parameters gamma for interaction;
#'	5)vector of other parameters (named), number of Genes, number of Treatment/Group, 
#'	probability for nonzero interactions;
#'	6)vector of data variance; 7)sample.size (balanced or unbalanced design)
#'	 
#'@examples 
#' 
#' nGene<-100
#' nTreatment<-2   #number of different beta
#' sampleSize<-5 ##this is the number of repeats for each group

#' alpha.mean<-0  #variance for alpha prior
#' alpha.sigma<-3
#' beta.mean<-0
#' beta.sigma<-2   #variance for beta prior
#'
#' gammaN0.sigma<-10  # the unscaled factor for gamma given gamma<=>0
#' p_diff<-0.01

#' #priors for variance distribution
#' d0<-5
#' s0<-2

#' #call it
#' set.seed(2004);
#' dataExp<-simulateExpression(nGene, nTreatment, sampleSize,
#'					alpha.mean=alpha.mean, alpha.sigma=alpha.sigma,
#'					beta.mean=beta.mean, beta.sigma=beta.sigma,
#'					prob.nonZero=0.01, gamma.sigma=gammaN0.sigma,
#'					epsilon.d0=d0, epsilon.s0=s0
#'					)
#'
#'@export

simulateExpression<-function(nGene, nTreatment, sampleSize=2,
			alpha=NULL, alpha.mean=NULL, alpha.sigma=NULL, 
			control.isotype=TRUE, control.negative=TRUE, 
			control.index=c(1),
			beta=NULL,beta.mean=NULL, beta.sigma=NULL, 
			gamma=NULL, prob.nonZero=NULL, gamma.sigma=NULL,
			epsilon.si=NULL, epsilon.d0=NULL, epsilon.s0=NULL
			)
			{
				#first need to check for the correct input
				#assuming the nGene and nTreatment and sampleSize
				#set.seed(2004)
				if(sampleSize<2)
				{
					stop("ERROR: sample size has to be bigger than 2!");
				}
				if(nTreatment<1)
				{
					stop("ERROR: number of treatment groups has to be bigger than 1!");
				}
				if(nGene<1)
				{
					stop("ERROR: number of genes has to be bigger than 1!");
				}
				if(nGene<100)
				{
					cat("WARNING: number of gene is small. The estimation might be inconsistent!\n")
				}
				
				if(is.null(alpha))
				{
					if(is.null(alpha.mean)||is.null(alpha.sigma))
					{
						stop("ERROR: no values has been set for alpha (the non-specific effects for genes!");
					}
					#now set alpha according to the distribution
					alpha<-rnorm(nGene,alpha.mean, alpha.sigma)
					alpha[1]<-0;
				}
				if(length(alpha)==1)
				{
					alpha<-rep(alpha, nGene)
				}
				if(length(alpha)!=nGene)
				{
					stop("ERROR:the length of input argument 'alpha' is not correct.")
				}
				#cat("alpha:", alpha,"\n")
				if(is.null(beta))
				{
					if(is.null(beta.mean)||is.null(beta.sigma))
					{
						stop("ERROR: no values has been set for beta (the non-specific effects for treatment group!");
					}
					#now set alpha according to the distribution
					beta<-rnorm(nTreatment,beta.mean, beta.sigma)
					if(control.negative)
					{
					#cat("yes\n")
						beta[control.index]=0;
					}
					if(control.isotype)
					{
					#cat("yes2\n")
						#beta[control.index]
						#here we don't have to do anything.
					}
				}
				if(length(beta)==1)
				{
					beta<-rep(beta,nTreatment);
				}
				if(length(beta)!=nTreatment)
				{
					stop("ERROR:the length of input argument 'beta' is not correct.")
				}
				#cat("beta:",beta,"\n")
				if(is.null(gamma))
				{
					if(is.null(prob.nonZero)||is.null(gamma.sigma))
					{
						stop("ERROR: no values has been set for alpha (the non-specific effects for genes!");
					}
					#now set alpha according to the distribution
					gamma<-matrix(rep(0,nGene*nTreatment),nrow=nGene, byrow=T)
					##in this following code snip, we randomly distribute the differential gamma into
					##different positions across different treatment group
					##other part is 
					numGene_diff<-floor(nGene*prob.nonZero)
					#cat("numGene_diff:",numGene_diff,"\n")
					antibodyIndex<-c(1:nTreatment)
					if(control.isotype||control.negative)
					{
						
						antibodyIndex<-antibodyIndex[-1*control.index]
					}
					cat("antibodyIndex:", antibodyIndex, "\n");
					for(i in antibodyIndex)
					{
						pos_diff<-sample(c(2:nGene), size=numGene_diff, replace=F)
						#cat("pos_diff:",pos_diff,"\n")
						#cat("i:",i,"\n")
						#cat("gamma sub:", gamma[pos_diff,i],"\n")
						#we want to be replace false, since otherwise we will have the collision (same number drawn multiple time)
						gamma[pos_diff,i]<-(rnorm(numGene_diff,0,gamma.sigma))
					}
					#gamma[0,]<-0;
					#if(control.isotype||control.negative)
					#{
					#	gamma[,index]<-0;
					#}
				}
				#cat("gamma:\n")
				#print(gamma)
				#cat("\n")
				if(length(gamma)==1)
				{
					gamma <-matrix(rep(gamma, nGene*nTreatment), nrow=nGene, byrow=T);
				}
				if(dim(gamma)[1]!=nGene||dim(gamma)[2]!=nTreatment)
				{
					stop("ERROR:the length of input argument 'gamma' is not correct.")
				}
				if(is.null(epsilon.si))
				{
					if(is.null(epsilon.d0)||is.null(epsilon.s0))
					{
						stop("ERROR: no distribution parameters has been set!");
					}
					##now ready to generate variance
					epsilon.si<-rScaledChisq(nGene*nTreatment,epsilon.d0,epsilon.s0*epsilon.s0)
					#epsilon.si<-matrix(rchisq(nGene*nTreatment, df=epsilon.d0),nrow=nGene, byrow=T)
					#epsilon.si<-1/(epsilon.d0*epsilon.s0*epsilon.s0)*epsilon.si
					#epsilon.si<-1/epsilon.si
					epsilon.si<-sqrt(epsilon.si)
					epsilon.si<-matrix(epsilon.si,nrow=nGene, byrow=T)
				
				}
				if(length(epsilon.si)==1)
				{
					epsilon.si<-matrix(rep(epsilon.si,nGene*nTreatment), nrow=nGene, byrow=T);
				}
				if(dim(epsilon.si)[1]!=nGene||dim(epsilon.si)[2]!=nTreatment)
				{
					stop("ERROR:the length of input argument 'epsilon.si' is not correct.")
				}
				##with the variance we can generate the data
				#now we have everything ready, do the observation values
				#put them into matrix first
				groupInfo<-rep(0,sampleSize*nTreatment)
				Y_ijk<-matrix(rep(0,nGene*nTreatment*sampleSize),nrow=nGene,byrow=T)
				for(j in c(1:nTreatment)) #for different treatment
				{
					samplePos<-c(1:sampleSize)+(j-1)*sampleSize;
					groupInfo[samplePos]<-j;
					for(k in c(1:nGene))
					{
						#write the meta data first
						Y_ijk[k,samplePos]<-alpha[k]+beta[j]+gamma[k,j]
						Y_ijk[k,samplePos]<-Y_ijk[k,samplePos]+rnorm(sampleSize,0,epsilon.si[k,j])
					}
				}

				#for now, return the matrix and will try other later
				retList<-list( "exp"=Y_ijk, "alpha"=alpha, "beta"=beta, 
						"gamma"=gamma, 
						"params"=c("nGene"=nGene, "nGroup"=nTreatment, "prob.nonZero"=prob.nonZero),
						"sample.Size"=sampleSize,  "epsilon.sd"=epsilon.si
						);
			}
			
#accessary function for converting the matrix to data.frame
#
#Oxygen 2 comment
#'@title convert a data matrix into data frame
#'@description takes in data matrix into a data frame with necessary 
#'	categorical information
#'@details It assumes a data matrix with a format of gene by groups
#'	and then turn it into a data frame with data fields of expression followed
#'	by group id and gene id.
#'@param x matrix containing the data with a format of genes as row and group as column
#'@param nGroup numeric the number of groups
#'@param sampleSize vector the number of repeats in each group. If it is a balanced 
#'	design, then a scalar value is enough
#'@return a data frame holding the data and categorical information
#'@export
matrix2dframe<-function(x, nGroup, sampleSize)
{
	if(missing(x))
	{
		stop("ERROR:undefined function argument \'x\'!")
	}
	if(missing(nGroup))
	{
		stop("ERROR:undefined function argument \'nGroup\'!")
	}
	if(missing(sampleSize))
	{
		stop("ERROR:undefined function argument \'sampleSize\'!")
	}
	
	#now check for the input
	if(!is.matrix(x)&&!is.data.fram(x))
	{
		#not a matrix
		stop("ERROR: the input data is not in a correct fromat!!")
	}
	#make sure the 
	if(length(sampleSize)==1)
	{
		sampleSize<-rep(sampleSize,nGroup)
	}
	#now check the size of group
	if(dim(x)[2]!=sum(sampleSize))
	{
		#something wrong
		stop("ERROR:the input data object doesn't have the correct dimension as specified!.")
	}
	#now , so far so good
	#do the conversion column by column
	dtm<-data.frame()
	gene<-c(1:dim(x)[1])
	index<-0
	for(i in 1:nGroup)
	{
		for(j in 1:sampleSize[i])
		{
			index<-index+1
			exp<-x[,index]
			group<-i
			temp<-data.frame(exp,group, gene)
			if(index==1)
			{
				dtm<-temp
			}
			else
			{
				dtm<-rbind(dtm, temp)
			}
		}
	}
	dtm$gene<-as.factor(dtm$gene)
	dtm$group<-as.factor(dtm$group)
	cat("Done!\n");
	return(dtm);
}


generateGammaMatrix<-function(nGene, nTreatment, prob.nonZero=0,group.differential=NULL)
{
	#now set alpha according to the distribution
	gamma<-matrix(rep(0,nGene*nTreatment),nrow=nGene, byrow=T)
	##in this following code snip, we randomly distribute the differential gamma into
	##different positions across different treatment group
	##other part is 
	numGene_diff<-floor(nGene*prob.nonZero)
	for(i in c(1:nTreatment))
	{
		pos_diff<-sample(nGene, size=numGene_diff, replace=T)
		gamma[pos_diff,i]<-(rnorm(numGene_diff,0,gamma.sigma))
	}
	return(gamma)
}