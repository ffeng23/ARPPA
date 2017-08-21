##=======gainAdjust.R==========
#
## This is the module to do the gainAdjust. The job is to based on the
## readings from multiple gain setting, "summarize" the reading in order
## to 1)minimize the saturation of the high reading points
## as well as 2)increase the low reading points.
##
## The model is based on a 5 parameter logistic function.
## 
## We assume all the readings should be well described a 
## common parameterized logistic function. The different
## protein concentration as well as the different affinity
## together determine/shift the points horizontally
## (on the x-axix). We just need to assume one
## function with count set of parameters (5 of them)
## and but also "align" by fitting to estimate 
## the shift/offset on the x axis
##
## Another thing to note is that we assume the maximum and 
## minimum values of the function is set 0 and 65535 respectively.
## We could allow them to vary and estimate it from the data, but
## proctically this fitting is not ideal. Fixing them makes the 
## fitting more reasonable and easy to converge. (***TODO: allow the   <==============
## maximum one fixed?? allowing minimum value to be varying??Need to
## try!!)
##	a=0.01, d=65535
##
## another thing to mention: the form of the logistic function we assume
## here is the one assume a log transformed x input, since it takes an exponetial
## in the formula. For y, it doesn't have to be log transformed. Either way, it should work,
## but it is said that log transformed is more symetric!!!
## 	We will try both anyway
##
###logistic function 5pl

#for the 
###model 5pl  --- deprecated ---obsolete
##  y<-d+(a-d)/(1+(x/c)^b)^g
### y - log values
### a - smallest y log
### d - hightest y log
### c - might point of x
### b - slope
### g - the asymetric factor
#####

###updated model now, is
###model 5pl
##  y<-d+(a-d)/(1+exp(xmid-x)^scal)^g
### y - log values
### a - smallest y log
### d - hightest y log
### xmid - middle point of x
### scal - slope at the middle point
### g - the asymetric factor
#####
##=======update 
###8/7/2017, move all the residual function to another file residualFunc.R
##
####======================================================


#------------start the code------------

##first define constants
#Maximizing the Dynamic Range of a Scan
#The complete dynamic range of the scanner is being used when you see a range of intensities on the image from 1 to 65535. A pixel with an intensity of 65535 is ¡°saturated¡±. Saturated pixels represent a condition in which there are more photons detected than the PMT can process, or the output of the PMTs exceeds the range of the analog-to-digital converters. A saturated pixel is not an accurate measurement of the signal from the pixel, so it is imperative to set the PMT Gain to avoid saturation.
#from GenePix_6.0_manual.pdf" page 40
pmt_saturated_reading<-65535;
pmt_lowest_reading<-1
a<-pmt_lowest_reading;
d<-(pmt_saturated_reading-5);

#'@include residualFunc.R


#this function is used to prepare the input of x
# in order to call the fitting function,
# since when we call nlsLM, it check for the input y and x,
# it will complain, when lengths of y and x are not identical
# y, dependent variable, a data matrix or data frame, each row is for one set of data
# x, independent variable, a vector of length of ncol of y

# pos, the position of unchanged/reference xs, 1 by default, 
##		this one is used ONLY when we want to move the reference series to the beginning
##		of the data. NOw this is not necessary. The code in the other
##		part is taking in the ref.index for the reference series.

# return value, a list of two vectors of equal length, y and x

##description, in this function, we prepare the input, transform the
## input into vectors (both y and x), and also based on "pos",
## we moved the referenced y and x into the first slot
## so then in the fitting module, we can assume this and 
## add the k, the shifting parameters, to the following 
## xs. 
### y must be either dataframe or matrix in a formate of genes by pmts
### x could be a vector indicating one set of pmts and
###			or it could be matrix or dataframe has same dimension of y
### we also do estimation of variace
#'@export
gainAdjust.prepareInput<-function(y, x, 
		#pos=1, 
		var.log=F
		)
	{
		#need to make sure ncol of y equals to length(x)
		if(class(y)=="numeric") #this is the one set of data, so return everything
		{
			#check for correct input
			if(length(y)!=length(x))
			{
				stop("the input y and x don't have correct format/length, please check!")
			}
			return(list("y"=y,"x"=x, "var"=1))
		}
		
		if(class(y)!="data.frame" && class(y)!="matrix")
		{
			stop("unsupported data format for input y!")
		}
		
	#now we are dealing with data matrix of frame
		
		#check for correct input
		if(class(x)=="matrix"|class(x)=="data.frame")
		{
			x<-as.matrix(x)
			if(dim(y)[2]!=dim(x)[2]|dim(y)[1]!=dim(x)[1])
			{
				stop("ERROR: the input y and x don't have correct format/length, please check!")
			}
		} else if(class(x)=="numeric")
		{
			if(dim(y)[2]!=length(x))
			{
				stop("ERROR: the input y and x don't have correct format/length, please check!")
			}
		}else
		{
			stop("unsupported data format for input x!")
		}
		
		#doing x first
		if(class(x)=="numeric")
		{
			x_t<-rep(x,dim(y)[1])
		} else
		{
			x_t<-rep(0,dim(x)[1]*dim(x)[2])
			for(i in c(1:dim(x)[1]))
			{
				x_t[c(1:dim(x)[2])+(i-1)*dim(x)[2]]<-x[i,]
				
			}
		}
		
		#now doing y
		y_t<-rep(0,dim(y)[1]*dim(y)[2])
		var_t<-rep(0,dim(y)[1]*dim(y)[2]/2)	
		for(i in c(1:(dim(y)[1]/2)))
		{
			y_t[c(1:dim(y)[2])+((i)*2-1-1)*dim(y)[2]]<-y[i*2-1,]
			y_t[c(1:dim(y)[2])+((i)*2-1)*dim(y)[2]]<-y[i*2,]
			if(var.log)
			{
				var_t[c(1:dim(y)[2])+((i-1))*dim(y)[2]]<-(log(y[i*2-1,]/y[i*2,]))^2
			}else {
				var_t[c(1:dim(y)[2])+((i-1))*dim(y)[2]]<-(y[i*2-1,]-y[i*2,])*(y[i*2-1,]-y[i*2,])
			}
		}
		
		##change the parameters
		#if(pos!=1)
		#{
		#	y1<-y_t[c(1:dim(y)[2])]
		#	y_t[c(1:dim(y)[2])]<-y_t[(pos-1)*dim(y)[2]+c(1:dim(y)[2])]
		#	y_t[(pos-1)*dim(y)[2]+c(1:dim(y)[2])]<-y1
		#	# have to do the same anything on x, if x is not a vector
		#	if(class(x)!="numeric")
		#	{
		#		x1<-x_t[c(1:dim(x)[2])]
		#		x_t[c(1:dim(x)[2])]<-x_t[(pos-1)*dim(x)[2]+c(1:dim(x)[2])]
		#		x_t[(pos-1)*dim(x)[2]+c(1:dim(x)[2])]<-x1
		#	}
		#}
		return (list("y"=y_t, "x"=x_t, "var"=var_t))
	}
	
	#------this is the simpliest way to prepare the inits for shifting
	#---- kind of arbitrary and not working well
	#---Depricated!!!
	#-- It works like this.
	#	--- #arbitrarily choose the first one in the data set as the reference
	#	---           all the lines shifted towards this one
	#	--- #to determine the init k values, we simply compare the biggest Y in 
	#	---				in each data series and take log of the quotient scaled by
	#	---				log(100), another aribitrary choice.
	#The accessary function for preparint intial values for nlsLM
	#this is necessary, because we might have many data series.
	#we need this for actually prepare the initial values
	
	#take in the input, before the prepareInput 
	#and make the initial values in order to do nlsLM fitting
	#it returns an array which is one less the total series (row lengths)
	#Always use the first one (row) as the reference
	#all the following rows are shifted left or right
	#input
	#	y is the input dataframe, has NOT been prepared. it has row as genes and col as pmts
	#	aggregated is indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
	#'@export
	gainAdjust.prepareInits<-function(y, aggregated=TRUE)
	{
		#get dimensions
		rlen<-dim(y)[1]
		clen<-dim(y)[2]
		yTemp<-y#y[seq(1,rlen,by=2),]
		####need to check for the aggregated status, in order to "arrange" data
		if(!aggregated)#need to summarize the data
		{	
			yTemp<-y[seq(1,rlen,by=2),]
			if(rlen%%2!=0)
			{
				stop("the number of array rows are not even, can not do aggregation! Stop!! ");
			}
			for(i in seq(1, rlen, by=2))
			{
				yTemp[(i+1)/2,]=(y[i,]+y[i+1,])
			}
		}
		rlen<-dim(yTemp)[1]
		ref<-yTemp[1,clen]
		
		inits<-rep(0,rlen-1)
		for(i in c(2:rlen))
		{
			inits[i-1]<-log(yTemp[i,clen]/yTemp[1,clen])/log(100)
		}
		inits
	}
	#------this is a more complicated way to prepare the inits for shifting
	#---- It takes a local linear model to determine the position of a line in 
	#		the reference data line.
	#
	#-- It works like this.
	#	--- #carefully choose the data set as the reference
	#	---           all the lines shifted towards this one.
	#	---				the reference line could be anywhere in the data
	#	---			To choose, we now use the criterion of MaxF/2 (65553/2). The closest
	#					Ymax will be chosed to be reference.
	#	--- #to determine the init k values, we simply compare the biggest Y in 
	#	---				in each data series. Two cases are possible
	#	---						1) Ymax_i <Ymax_r, nothing to do
	#	---						2) Ymax_i >Ymax_r, find the max Yj_i in the series to be smaller then Ymax_r
	#	---					Now we have a pair (Yj_i, Xj_i).
	#	---				Determine where does this pair belongs to inside the Yr data line
	#						Find (Y(k-1)_r, Yk_r), where Y(k-1)_r<Yj_i<Yk_r. Then 
	#					simply using a linear model to determine Xj_i_pred for Yj_i according
	#						to reference data . It could possible happen that Yj_i is not inside
	#						reference range, meaning Yj_i<Ymin_r. In this case, we simply assuming
	#						Yj_i == Ymin_r, using X(Ymin_r) as the predication and shift the data.
	##						Hopefully, there will not be many cases like this.
	#	--- #to shift, we take the log(Xj_i/Xj_i_pre)
	#The accessary function for preparint intial values for nlsLM
	#this is necessary, because we might have many data series.
	#we need this for actually prepare the initial values
	
	#take in the input, before the prepareInput 
	#and make the initial values in order to do nlsLM fitting
	#it returns an array which has length identifical to the total series (row lengths)
	#
	#input
	#	y is the input dataframe, has NOT been prepared. it has row as genes and col as pmts
	#	aggregated is indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
	#
	#Note: the prepared inits vector always is in a format of aggregated data
	#'@export
	gainAdjust.prepareInitsLM<-function(y, x, aggregated=TRUE)
	{
		#get dimensions
		rlen<-dim(y)[1]
		clen<-dim(y)[2]
		yTemp<-y#y[seq(1,rlen,by=2),]
		####need to check for the aggregated status, in order to "arrange" data
		if(!aggregated)#need to summarize the data
		{	
			yTemp<-y[seq(1,rlen,by=2),]
			if(rlen%%2!=0)
			{
				stop("the number of array rows are not even, can not do aggregation! Stop!! ");
			}
			for(i in seq(1, rlen, by=2))
			{
				yTemp[(i+1)/2,]=(y[i,]+y[i+1,])
			}
		}
		rlen<-dim(yTemp)[1]
		#ref<-yTemp[1,clen]
		
		#sort x in case X is not in order
		idx.srtX<-order(x)  #increasing order
		
		inits<-rep(0,rlen)
		Y<-yTemp[,idx.srtX]
		X<-x[idx.srtX]
		idx.ref<-which.min(abs(pmt_saturated_reading-yTemp[,clen]))
		#idxs.arr<-yTemp[-1])
		ref<-yTemp[idx.ref,idx.srtX]
		for(i in c(1:rlen))
		{
			if(i==idx.ref) #this is the reference array, skip to the next,
			{
				next
			}
			#do the job to calculate the shifting ks
				#first compare the Ymax with Ymax_r
			#want to find a Yi that is smaller than ref_max
			Yi<-Y[i,]
			
			j<-1;
			for(j in length(Yi):1)
			{
				if(Yi[j]<=ref[clen])
				{
					break; ##done,
				}
				
			}
			Yi_j<-Yi[j]
			Xi_j<-X[j]
			if(Yi_j>ref[clen]) #even the smaller Y in the search one is bigger than the ref max, what to do?? 
			{
				Xi_r_pred<-X[clen]
			}
			else
			{
				#now that, we find a pair (Yi, Xi) to compare with ref data
				#next need to determine, which region it locates in ref data
				k<-1
				for(k in (clen-1):1)
				{
					if(Yi_j>=ref[k])
					{
						break; #we found a good one 
					}
				}
				
				if(Yi_j>= ref[1]) #do linear model
				{
					a<-(ref[k]-ref[k+1])/(X[k]-X[k+1])
					b<-((ref[k]/X[k])-(ref[k+1]/X[k+1]))/(1/X[k]-1/X[k+1])
					Xi_r_pred<-(Yi_j-b)/a
				}
				else  #Yi is smaller than smallest ref
				{
					Xi_r_pred<-X[k] #in this case, k=1
				}
				
			}
			inits[i]<-log(Xi_r_pred/Xi_j)
		}
		list("inits"=inits, "ref.index"=idx.ref)
	}
	
	#the five parameter logistic function
	#take in the parameter 5 parameters and get the function values
	#'@export
	f5pl<-function(pars,x)
	{
		if(length(pars)!=5)
			stop("pars are not correctly specified!")
		a<-pars[1]
		d<-pars[2]
		xmid<-pars[3]
		scal<-pars[4]
		g<-pars[5]
		
		pred<-a+(d-a)/(1+exp((xmid-x)/scal))^g
	}
	
	#need to plot as a function
	#y is the original dataframe and has not been "changed"
	#x is the pmt single array
	#model, in this case, we are using 5PL, but with fixed highest end (fix.high), fixed both end (fix.both) and fixed none (fix.none)
	#'@export
	gainAdjust.plot<-function(LMfit,y, x, ylog=TRUE,aggregated=FALSE,model="fix.both",
			a,d, xmid, scal, g, ref.index=1,
			filename=NULL)
	{
		if(missing(LMfit))
		{
			stop("no LM nls fitting object has been specified!")
		}
	
		if(class(LMfit)!="nls.lm"&&class(LMfit)!="list")
		{
			stop("the input is not valid. A nls.lm or list object with par field is needed!")
		}
		#prepare the parameters
		#determine the parameters first
		#define the default case of 3PL
		#xmid<-LMfit$par[1]
		#scal<-LMfit$par[2]
		#g<-LMfit$par[3]
		k_index<-4
		if(missing(a))
			{
				a<-pmt_lowest_reading
			}
		if(missing(d))
			{
				d<-pmt_saturated_reading
			}
		if(ylog)
		{
			a<-log(a)
			d<-log(d)
		}
		switch(model,
				"fix.both"={	
						#this has been taken care of as the default case in above
						xmid<-LMfit$par[1]
						scal<-LMfit$par[2]
						g<-LMfit$par[3]
						#do it again here, because we can not leave it empty
						k_index<-4
					},
				"fix.high"={
						#be careful here, 4PL is not the real 4PL, it is 5PL but fixed d, maximum level
						a<-LMfit$par[1]
						xmid<-LMfit$par[2]
						scal<-LMfit$par[3]
						g<-LMfit$par[4]
						k_index<-5
					},
				"fix.none"={
						#all five parameter
						a<-LMfit$par[1]
						d<-LMfit$par[2]
						xmid<-LMfit$par[3]
						scal<-LMfit$par[4]
						g<-LMfit$par[5]
						k_index<-6
				},
				"fix.all"={
					if(missing(xmid)|missing(scal)|missing(g))
					{
						stop("please specify parameters (xmid/scal/g)!!!")
					}
					k_index<-1
					
				},
				stop("unkown model specified")
			)
		#now doing the plot
		if(!is.null(filename))
		{
			jpeg(filename=filename)
		}
		#op<-par(mfrow=c(2,1), mar = c(3,5,2,1))
		#need to figure out the range of y and x
		#for y, we know it is (0.1, 65535)
		
		#for x, we need to add ks
		xmin<-log(x[1])+min(LMfit$par[c(k_index:length(LMfit$par))])
		xmax<-log(x[length(x)])+max(LMfit$par[c(k_index:length(LMfit$par))])
		
		plot(c(xmin,xmax+0.3), c(min(log(a*0.95),min(y)),log(d)+0.2), bg = "black", cex = 0.5, main="data", type="n")
		#plot(c(xmin,xmax+0.3), c((0.1),(66535)+5000), bg = "black", cex = 0.5, main="data", type="n")
		#points(log(x), log(y[1,]), col="grey", lty=2) 
		params<-LMfit$par[c(k_index:length(LMfit$par))]
		if(ref.index==1){
			params<-c(0,params)
		}else if (ref.index==length(params+1)){
			params<-c(0,params)
		}else{
			params<-c(params[1:(ref.index-1)],0,params[ref.index:length(params)])
		}
		yParams<-y #[-1,]
		if(!aggregated)
		{
			#cat("calling it now for aggregation")
			#yParams<-yParams[-1,];
			#params<-rep(params, rep(2,length(params)));
			yParams<-sqrt(y[seq(1, dim(y)[1],by=2),]*y[seq(2, dim(y)[1],by=2),])#/2
			if(length(params)!=dim(yParams)[1])
			{
				stop("Error: the length of params is not identical to number of y data points!")
			}
			#for(j in c(1:length(params)))
			#{
			#	yParam[j,]<-(y[j*2-1,]+y[j*2,])/2
			#
			#}
			
		}
		for(i in c(1:length(params)))
		{
			#if(i!=304)
			#{
			# next;
			#}
			#cat("\ny:", yParams[i,], "\n")
			#cat("param[i]:", params[i], "\n")
			#points(log(x)+params[i], log(yParams[i,]), col="grey", lty=2)
				if(ylog) #this means that the data is logged, so don't do log again
				{
					lines(log(x)+params[i], yParams[i,], col=i, lty=2) 
				}
				else
				{
					#------>
					#cat("-->long:",log(yParams[i,]), "\n" )
					lines(log(x)+params[i], log(yParams[i,]), col=i, lty=2) 
				}
		}
			
		#now plot the fitted value
							
		xx<-seq(xmin,xmax+0.3, by=(xmax+0.3-xmin)/1000)
		#cat(xx,"\n")
		yy<-a+(d-a)/(1+exp((xmid-xx)/scal))^(g)
		#cat(yy, "\n")
		cat("a:", a, "\td:", d, "\txmid:", xmid, "\tscal:",scal, "\tg:", g, "\n")
		if(!ylog)
		{
			#---->
			yy<-log(yy)
		}
		lines(xx,yy,col=2, lty=1,lwd=2)
		#--plot(1:(LMfit$niter+1), log(LMfit$rsstrace), type="b",
		#--main="log residual sum of squares vs. iteration number",
		#--xlab="iteration", ylab="log residual sum of squares", pch=21,bg=2) 
		#par(op)
		if(!is.null(filename))
		{
			dev.off();
		}
	}
	##do aggregated the duplicated data, either by arithmatic or geometric mean
	#'@export
	gainAdjust.aggregate<-function(y,  mode="geometric")
	{
		if(missing(y))
		{
			stop("data missing, please specify")
		}
		y.aggregated<-c();
		
		switch(class(y),
		"matrix"={
			nrow<-dim(y)[1]
			#ncol<-dim(y)[2]
			
			#assuming the duplicated data are next to each other
			if(nrow!=floor(nrow/2)*2)
			{
				stop("data is corrupted. the number of data rows is not even!")
			}
			idx<-seq(1, nrow, by=2)
			y.aggregated<-y[idx,]
			if(mode=="geometric") #this means that the data is logged, so don't do log again
				{
					y.aggregated<-sqrt(y.aggregated*y[idx+1,]) 
				}
			else if(mode=="arithmatic")
				{
					#------>
					#cat("-->long:",log(yParams[i,]), "\n" )
					y.aggregated<-(y.aggregated+y[idx+1,])/2
				}
			else
				{
					stop("unknow mode, specify either arithmatic or geometric")
				}
		
		}, #case one
		#"numeric"={
		#	ylen<-length(y)
		#	#assuming the duplicated data are next to each other
		#	if(ylen!=floor(ylen/2)*2)
		#	{
		#		stop("data is corrupted. the number of data rows is not even!")
		#	}
		#	idx<-seq(1, ylen, by=2)
		#	y.aggregated<-y[idx]
		#	if(mode=="geometric") #this means that the data is logged, so don't do log again
		#		{
		#			y.aggregated<-sqrt(y.aggregated*y[idx+1]) 
		#		}
		#	else if(mode=="arithmatic")
		#		{
		#			#------>
		#			#cat("-->long:",log(yParams[i,]), "\n" )
		#			y.aggregated<-(y.aggregated+y[idx+1])/2
		#		}
		#	else
		#		{
		#			stop("unknow mode, specify either arithmatic or geometric")
		#		}
			
		#}, #case two
		{
			stop("Not correct data format. data matrix is expected!")
		} #default case
		)
		
		y.aggregated
	}
	
	#this is the function to align the data series into one
	#5pl lines based on the fiting with shifted parameters
	#the input is 
	# LMfit, with the shifted parameters
	#	y in datafame format
	#	x the original x data series in PMTs
	#return value
	#	dataframe for shifted xs, in order of y . And 
	#   this is not log transformed
	#'@export
	gainAdjust.alignData<-function(LMfit,x, aggregated=F,model="fix.both", ref.index=1)
	{
		k_index<-4
		
		switch(model,
				"fix.both"={	
						#do it again here, because we can not leave it empty
						k_index<-4
					},
				"fix.high"={
						#be careful here, 4PL is not the real 4PL, it is 5PL but fixed d, maximum level
						k_index<-5
					},
				"fix.none"={
						#all five parameter
						k_index<-6
				},
				"fix.all"={
						k_index<-1
				},
				stop("unkown model specified")
			)
			params<-c(LMfit$par[c(k_index:length(LMfit$par))])
		if(ref.index==1){
			params<-c(0,params)
		}else if(ref.index==length(params)+1){
			params<-c(params,0)
		}else{
			params<-c(params[1:ref.index-1], 0, params[ref.index:length(params)])
		}
		
		
		if(!aggregated)
		{	
			params<-rep(params, rep(2,length(params)));
				
		}
		params<-rep(params,rep(length(x),length(params)))
		
		params<-matrix(params,ncol=length(x),byrow=T)
		xFrame<-data.frame(matrix(rep(x,dim(params)[1]),ncol=length(x),byrow=T))
		xFrame*exp(params)
		
	}
	
	#'@export
	fitShifts5plSingle<-function(ydata, x, par.5pl, data.aggregated=F, debug=F)
	{
		if(missing(ydata))
		{
			stop("please specify the expression data!")
		}
		if(missing(x))
		{
			stop("please specify the PMT setting data!")
		}
		if(missing(par.5pl))
		{
			stop("please specify the fitting parameters!")
		}
		#if(length(var)!=1&&length(var)!=dim(ydata)[1])
		#{
		#	stop("the var is not set up correctly")
		#}
		initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
				
			#-gainAdjust.fitFinal<-nls.lm(par=initsLM$inits[-initsLM$ref.index], 
			#-	fn=gainAdjust.fnResShift,
			#-	#y=log(yinput.aligned), 
			#-	y=(yinput), 
			#-	x=log(xinput), xlen=xlen,aggregated=data.aggregated,
			#-	ylog=data.ylog, model.weight="log", #in here, this fnResShift doesn't take any model, it hardcodes using log transform
			#-	order=1,d=(d),a=(a), 
			#-	xmid=gainAdjust.fit$par[1],scal=gainAdjust.fit$par[2],g=gainAdjust.fit$par[3],
			#-	shift.index=initsLM$ref.index,
			#-	control = nls.lm.control(nprint=1, maxit=100)
			#-	)
		#the statement block for running fitting with inividual data sets	
		#{
			fitShift<-list();
			fitShift$par<-initsLM$inits;#[-initsLM$ref.index]; #initialize the output
			#--->xsim<-seq(4,7,by=0.01)
			#--->plot(exp(xsim), f5pl(c(a,d, gainAdjust.fit$par[1:3]),xsim),col=2, lwd=2, type="l", log="xy", lty=2)
			cat("Fitting the shifting parameters......\n")
			nprint<-0;
			if(debug)
			{
				nprint<-1;
			}
			for(i in 1:length(initsLM$inits))
			{
				if(debug){
					cat(i,"/", length(initsLM$inits),"...\n")
				}
				if(i==initsLM$ref.index)
				{
					fitShift$par[i]<-0
					next;
				}
				xlen<-length(x)
				#yvar<-(log(ydata[i*2-1,]/ydata[i*2,]))^2
				#idx.zero<-which(yvar==0)
				#if(length(idx.zero)>0)
				#{
				#	yvar[idx.zero]<-min(yvar[-idx.zero]) ##to get rid of zero variance.
				#}
				#yvar<-1
				#do fitting individually
				
				gainAdjust.fitS<-nls.lm(par=fitShift$par[i], fn=gainAdjust.fnResShiftSingle,
					#y=log(yinput.aligned), 
					y=c(ydata[i*2-1,], ydata[i*2,]), 
					x=log(x), #xlen=9,
					aggregated=data.aggregated,
					#ylog=ylog, model.weight="power", #in here, this fnResShift doesn't take any model, it hardcodes using log transform
					#order=1,
					a=par.5pl[1],d=par.5pl[2], 
					xmid=par.5pl[3],scal=par.5pl[4],g=par.5pl[5], mvar=1,#yvar,
					#shift.index=initsLM$ref.index,
					control = nls.lm.control(nprint=nprint, maxit=100))
					
					 #
				if(debug){
					lines(exp(log(x)+gainAdjust.fitS$par[1]),ydata[i*2,],col=i)
					lines(exp(log(x)+gainAdjust.fitS$par[1]),ydata[i*2-1,],col=i)
				}

				fitShift$par[i]<-gainAdjust.fitS$par[1]
			}#end of for loop
		#}
		fitShift$par<-fitShift$par[-initsLM$ref.index]
		fitShift
	}
	
	
	###------------>updated 8/4/2017, feng
	#the function to fit 5PL for parameters and shifts for INDIVIDUAL array with multiple PMT readings
	#In this function, we first 1)do fitting with fixed a and d, with a model of uniform
	#then we 2) do fitting with the 5PL parameters to shifted x's under log-ed data.
	#Then in the third step, 3) estimate 5PL parameters, a, xmid, scal and g (fixing d) with shifted data
	# under sqrt weighted matrix.
	#	If in the last step, the fitted a is negative, we just simply go back to the original fitted 5PL. 
	## NOTE: somehow, the first fitting fixing both has to be uniform weights to get the best fitting.
	## 
	##8******TODO, get a better init for parameters in the 1st 3p fitting
	###Input: y, is the list containing expression and control expression
	##			y$C the data matrix for control expression in array with different PMT readings, 
	#			y$E the data matrix for row is for target expression in array with different PMT readings
	#			y$gene  meta data for gene 
	#			y$cgene meta data for control gene
	### 	x, the vector containing the PMT gain. In the fitting, we always want to do log(xdata)
	###			elements in x must be arranage in a order identical to the column order of PMT settings for array in
	###			
	###		 weight.mode are for advanced users. modify them if you know what you are doing.
	###		data.ylog and data.aggregated are specifying data status. ylog, the y data are in log 
	###			and expression is aggregated, meaning the two replicates are averaged arithmatic and geometric.
	###			(the mode of aggregation was done by user)
	###		fit.aggregated and aggrgated.mode are describing how to do fitting 
	##			either using aggregated or replicated data.  Default is raw non-aggregated data. 
	#			The mode is either arithmatic and geometric. 
	###		fit.mode, to specify whether to use control (control) or control+Expression data (all). By default, only
	#			fit the control data. This is for time issue, since including more data takes too long.
	#		block.size, the size of block. The total block is 48. 
	##			The reason biggest size is 12 control blocks, with about 1000 points (500 unique points)
	##		threshold, the lowest value that maximum PMT gain setting expression should achieve. If data
	#			points with all expression level below this threshold are not included for fitting 5PL
	#
	
	#'@export
	gainAdjust.fit5plArray<-function(y, x, data.ylog=F, data.aggregated=F,
		fit.mode="control",
		fit.aggregated=F, aggregated.mode="geometric", 
		weight.mode, order=1,
		block.size=12,a=1.1, d=65535,threshold=30,
		PMT.gain=NULL,
		debug=F
		)
	{
		if(missing(x))
		{
			stop("input \"x\" missing, please specify!");
		} 
		if(missing(y))
		{
			stop("input \"x\" missing, please specify!");
		}
		if(class(y)!="list"){
			stop("ERROR: the input data of y is not in a correct format. List with expression data needed!!")
		}
		if(missing(weight.mode))
		{
			weight.mode<-c("uniform", "log", "sqrt");
		}
		#check for the correct weight mode
		for(i in 1:length(weight.mode))
		{
			if(!is.element(weight.mode[i], c("log", "uniform", "power","exp", "sqrt")))
			{
				stop("Unknow weight mode specified. Only log/uniform/power/exp supported!");
			}
		}
		#need to get data in below based on input,
		blkNum<-max(y$gene[,"Block"])
		if(block.size<1||block.size>blkNum)
		{
			stop(paste("only a block size between 1 and ", blkNum,sep=""))
		}
		#prepare the output
		fr<-list();#this is the fitting results
		dr<-list("E"=y$E[,1], "C"=y$C[,1], "E.adj"=y$E[,1], "C.adj"=y$C[,1], 
				"E.gain"=y$E[,1], "C.gain"=y$C[,1]); #this is the gain adjust data output, E/C
															#is the array holding the optimal data for this 
															#array, and PMT.gain is the gain setting
															#for the data
																	
		###now good ?? start doing the fitting, based on the control and block size
		for(j in 1:ceiling(blkNum/block.size)){
			cat("Fitting data by blocks: ",j, "/",blkNum/block.size, "...\n")
			#for the array, we need to get the protein data by block, 
			bg<-(j-1)*block.size+1
			ed<-(j-1)*block.size+block.size
			if(ed>blkNum)
			{
				ed<-blkNum
			}
			
			idx.prC<-which(is.element(y$cgene[,"Block"],c(bg:ed)))
			ydata<-y$C[idx.prC,]
			if(fit.mode=="all")
			{
				idx.prE<-which(is.element(y$gene[,"Block"],c(bg:ed)))
				ydata<-rbind(y$E[idx.prE,], ydata)
			}
			#now data is ready, rm the low expressed protein, since they are confusing the fitting
			ydata5PL<-rmLows(ydata, threshold=threshold, index=which.max(x), aggregated=data.aggregated)
			xlen<-length(x)
			
			#now we need to do aggregation
			if(!data.aggregated)
			{
				if(fit.aggregated)
				{
					#do aggregation on data
					if(data.ylog)
					{	
						if(aggregated.mode=="geometric")
						{
							cat("WARNING: call to run aggregation on log transformed data. Could lead to NaN values")
						}
					}
					data.aggregated<-T
					ydata5PL<-gainAdjust.aggregate(ydata5PL, aggregated.mode)
				}
			}else  ##input data already aggregated, then do nothing
			{
				#in this case, we ignore whatever "fit.aggregated" says, since there are no way to revert data back
				#so set data.aggregated as true
				data.aggregated<-T; #redundant
			}
			
			
			#now we have the data, reay, aggregated
			input<-gainAdjust.prepareInput(ydata5PL,x)
			yinput<-input$y
			xinput<-input$x
			
			#get the initial values for fitting, by picking the data series within the middle range 
			initsLM<-gainAdjust.prepareInitsLM(ydata5PL,x, aggregated=data.aggregated)
			#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
			
			#a<-a
			#d<-d
				
			#for three parameters, xmid, scal and g
			parStart<-c(6.5,0.005,0.015,initsLM$inits[-initsLM$ref.index]) #<---3pl, keeping a and d constant
			if(data.ylog)
			{
				a<-log(a)
				d<-log(d)
				parStart<-c(6.3,0.1,0.1,initsLM$inits[-initsLM$ref.index]) 
			}
			cat("\tinitial fitting for shifts and parameters.....\n")
			gainAdjust.fit<-nls.lm(par=parStart, 
				#fn=gainAdjust.fnRes4pShift,  #<----4pl, varying a
				fn=gainAdjust.fnRes3pShift, #<----3pL, keeping a and d constant
				#y=log(yinput), 
				y=yinput,
				x=log(xinput), xlen=xlen, aggregated=data.aggregated, ylog=data.ylog,a=a,
				shift.index=initsLM$ref.index,
				model.weight="uniform", order=1, control = nls.lm.control(nprint=1,maxiter=50)
				)
			if(debug){	
				#sink("debug.txt", append=T)
				cat("summary of block", j, ":\n")
				summary(gainAdjust.fit)
				#sink();
			
			#gainAdjust.plot(gainAdjust.fit,log(y), x,ylog , ref.index=initsLM$ref.index)
			#gainAdjust.plot(gainAdjust.fit,y, x,ylog)
			gainAdjust.plot(gainAdjust.fit,ydata5PL, x, ylog=data.ylog , 
								ref.index=initsLM$ref.index,a=a, aggregated=data.aggregated
								#,filename=paste("fit1_blk_",j,".jpg",sep="")
							)
			#plot 
			x.aligned<-gainAdjust.alignData(gainAdjust.fit,x, model="fix.both", ref.index=initsLM$ref.index)##<---for 4Pl, 
			input.aligned<-gainAdjust.prepareInput(ydata5PL,x.aligned)
			yinput.aligned<-input.aligned$y
			xinput.aligned<-input.aligned$x
			jpeg(filename=paste("fit1_dot_blk_",j, ".jpeg",sep=""))
			plot((xinput.aligned), (yinput.aligned), type="p", main="5-p logistic fitting of intensity vs. PMT gain",
				xlab="PMT gain (log)",ylab="intensity (log)",log="xy"
				)
			#plot(log(xinput.aligned), (yinput.aligned), type="p")
			xsim<-seq(4,7,by=0.01)
			
			lines(exp(xsim), (f5pl(c(a,d, gainAdjust.fit$par[1:3]),xsim)),lty=1,col=3,lwd=2)
			dev.off()
			}
		#done for first run
			
		#now do second fit, only shift data
			#do aggregation now
			#ydata_ag<-gainAdjust.aggregate(ydata, mode="geometric")
			#aggregated_3rd<-T
			#--input<-gainAdjust.prepareInput(ydata,x) #now including even the low expressing data points
			#--yinput<-input$y
			#--xinput<-input$x
			cat("\trefining fitting for shifts .....\n")	
		#get the initial values for fitting, by picking the data series within the middle range 
			initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			#parStart<-c(6,0.1,0.13,-0.2,-0.5,-0.55, -0.5,1)
				
			fitShift<-fitShifts5plSingle(ydata,x, c(a,d,gainAdjust.fit$par[1:3]), data.aggregated, debug);
		#}
			
			if(debug){
			gainAdjust.plot(fitShift,ydata, x,ylog=data.ylog , ref.index=initsLM$ref.index,
				model="fix.all",aggregated=data.aggregated,
				a=a, d=d, xmid=gainAdjust.fit$par[1], scal=gainAdjust.fit$par[2], g=gainAdjust.fit$par[3]
				#,filename=paste("fit2_blk_",j,".jpg", sep="")
				)#a=1.2440093153044
			
				
			x.aligned<-gainAdjust.alignData(fitShift,x, model="fix.all", ref.index=initsLM$ref.index)##<---for 4Pl, fixing highest end only

			input.aligned<-gainAdjust.prepareInput(ydata,x.aligned)
			yinput.aligned<-input.aligned$y
			xinput.aligned<-input.aligned$x		
			jpeg(filename=paste("fit2_dot_blk_",j,".jpg", sep=""))		
			plot((xinput.aligned), (yinput.aligned), type="p", main="5-p logistic fitting of intensity vs. PMT gain",
				xlab="PMT gain (log)",ylab="intensity (log)",log="xy"
				)
			#plot(log(xinput.aligned), (yinput.aligned), type="p")
			xsim<-seq(4,7,by=0.01)
			
			lines(exp(xsim), (f5pl(c(a,d, gainAdjust.fit$par[1:3]),xsim)),lty=1,col=3,lwd=2)
			dev.off();
			}			
		#done for second try
			
		#now do the last fitting for parameters, mainly for parameter a
		#-------->now fit 4 parameters on aligned data again to allow a to vary
			ydata_ag<-gainAdjust.aggregate(ydata, mode=aggregated.mode)
			aggregated_3rd<-T
			x.aligned<-gainAdjust.alignData(fitShift,x, aggregated=T, model="fix.all", ref.index=initsLM$ref.index)##<---for 4Pl, fixing highest end only

			input.aligned<-gainAdjust.prepareInput(ydata_ag,x.aligned)
			yinput.aligned<-input.aligned$y
			xinput.aligned<-input.aligned$x		
			parStartFpl<-c(a,
							gainAdjust.fit$par[1],gainAdjust.fit$par[2],gainAdjust.fit$par[3])
			cat("\tfinal fitting for parameters.....\n")
			gainAdjust.fitFinalP<-nls.lm(par=parStartFpl, 
					fn=gainAdjust.fnRes4pFpl,
					#fn=gainAdjust.fnRes3pFpl,
					#y=log(yinput.aligned), 
					y=yinput.aligned, 
					x=log(xinput.aligned),#xlen=xlen, aggregated=aggregated_3rd,
					ylog=ylog, model.weight="sqrt",order=2,d=d,#a=0.1,
					control = nls.lm.control(nprint=1, maxit=100)
			)
			if(gainAdjust.fitFinalP$par[1]<0)
			{
				gainAdjust.fitFinalP$par[1:4]<-c(a, gainAdjust.fit$par[1:3])
			}
			if(debug){
				gainAdjust.plot(fitShift,ydata_ag, x,ylog=data.ylog , ref.index=initsLM$ref.index,
						model="fix.all",
						a=gainAdjust.fitFinalP$par[1], d=d, 
						xmid=gainAdjust.fitFinalP$par[2], scal=gainAdjust.fitFinalP$par[3], g=gainAdjust.fitFinalP$par[4],
						aggregated=T
						#,filename=paste("fit3_blk_",j,".jpg")
				)#a=1.2440093153044
				
				jpeg(filename=paste("fit3_dot_blk_",j,".jpeg",sep=""))
				plot((xinput.aligned), (yinput.aligned), type="p", main="5-p logistic fitting of intensity vs. PMT gain",
					xlab="PMT gain (log)",ylab="intensity (log)",log="xy"
					)
				#plot(log(xinput.aligned), (yinput.aligned), type="p")
				xsim<-seq(4,7,by=0.01)
				
				lines(exp(xsim), (f5pl(c(gainAdjust.fitFinalP$par[1],d, gainAdjust.fitFinalP$par[2:4]),xsim)),lty=1,col=3,lwd=2)
				dev.off();
			}
			cat("\tDONE!\n");
			#shift<-
			#fitRes<-list("5PL"=c(gainAdjust.fitFinalP$par[1],d, gainAdjust.fitFinalP$par[2:4]), 
			#	"shift"=gainAdjust.insertAt(fitShift$par,0, initsLM$ref.index))
			#fr[[j]]<-fitRes;
			
			#from this point one, we start doing the "normalization" to pick the good ones for each array
			#first do control
			if(fit.mode=="all")
			{
				ydata<-y$C[idx.prC,]
				initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			}
			#get shift
			fitShift<-fitShifts5plSingle(ydata,x, c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]), data.aggregated);
			fitShift$par<-gainAdjust.insertAt(fitShift$par, 0, initsLM$ref.index)
			#generate output
			temp_gain<-NULL
			if(!missing(PMT.gain)&&!is.null(PMT.gain)){
				temp_gain<-PMT.gain$C[idx.prC]
			}
			#cat("1\n")
			ydataC.adj<-adjust.Matrix(ydata=ydata,x=x,par.5pl=c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]),fitShift=fitShift,PMT.gain=temp_gain, F);
			#cat("2\n")
			fr$C[idx.prC]<-ydataC.adj$E;
			fr$C.adj[idx.prC]<-ydataC.adj$E.adj;
			fr$C.gain[idx.prC]<-ydataC.adj$gain;
			
			#Now do fitting for shifts for target expression only
			idx.prE<-which(is.element(y$gene[,"Block"],c(bg:ed)))
			ydata<-(y$E[idx.prE,])
			initsLM<-gainAdjust.prepareInitsLM(ydata,x, aggregated=data.aggregated)
			temp_gain<-NULL
			fitShift<-fitShifts5plSingle(ydata,x, c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]), data.aggregated);
			fitShift$par<-gainAdjust.insertAt(fitShift$par, 0, initsLM$ref.index)
			if(!missing(PMT.gain)&&!is.null(PMT.gain)){
				temp_gain<-PMT.gain$E[idx.prE]
			}
			ydataE<-adjust.Matrix(ydata=ydata,x=x,par.5pl=c(gainAdjust.fitFinalP$par[1],d,gainAdjust.fitFinalP$par[2:4]),fitShift=fitShift,PMT.gain=temp_gain,F);
			fr$E[idx.prE]<-ydataE$E;
			fr$E.adj[idx.prE]<-ydataE$E.adj;
			fr$E.gain[idx.prE]<-ydataE$gain;
		} #end of for loop for analysis by block
	
	#prepare the output
		
		fr
	}#end of the function
	
	#this is the function determine the
	#fitShift including all the shifts for all data points for each individual Array data E or C
	#	This function includes the one with index of initsLM$ref.index.
	#
	#'@export
	adjust.Matrix<-function(ydata, x, par.5pl,  fitShift, PMT.gain, data.aggregated=F )
	{
		if(!data.aggregated)
		{
			fitShift$par<-rep(fitShift$par, rep(2, length(fitShift$par)))
		}
		##check for the data integrity
		if(missing(PMT.gain)||is.null(PMT.gain)){#this is the case we need to compare to get optimal pmt for each gene
			PMT.gain<-fitShift$par
			#now go through the list to pick the  best one
			xmid<-par.5pl[3]
			for(i in 1:length(fitShift$par))
			{
				PMT.gain[i]<-x[which.min(abs(log(x)+fitShift$par[i]-xmid))]
			}
		}
		#now for sure, we have PMT.gain reference setting, get 
		y.out<-list("E"=c(),"E.adj"=c(), "gain"=PMT.gain)
		y.out$E.adj<-f5pl(par.5pl,(log(PMT.gain)+fitShift$par))
		y.out$E<-fitShift$par
		for(j in 1:dim(ydata)[1])
		{
			y.out$E[j]<-ydata[j,which(x==PMT.gain[j])]
		}
		y.out
	}
	
	#function to insert value at "index" into an array
	#'@export
	gainAdjust.insertAt<-function(vec, value=0, index=1)
	{
		if(missing(vec))
		{
			stop("missing the input vector, please specify")
		}
		if(class(vec)!="numeric")
		{
			stop("the input vector is not in correct format. please specify")
		}
		
		arr<-c(0, vec)
		if(index==1)
		{
			arr[1]<-value
		}else if(index==length(vec)+1)
		{
			arr[1:length(vec)]<-vec
			arry[length(vec)+1]<-value
		}else  #in this case, index is not on th boundary
		{
			arr[1:(index-1)]<-vec[1:(index-1)]
			arr[index]<-value
			arr[(index+1):length(arr)]<-vec[index:length(vec)]
		}
		arr
	}
	
	###function to do 
	##assuming has done background correction
	##input, 
	###	x, is the PMT gain setting vector/matrix assuming all the arrays having the identical pmt vector
	##
	#'@export
	gainAdjust.fitArrayByBlock<-function(object, x, ylong=F, aggregated=F, debug=F )
	{
		if(missing(object))
		{
			stop("missing object, please specify....")
		}
		if(class(object)!="EListRaw")
		{
			stop("object is not of correct format. please specify......")
		}
	
		#now getting data ready to do the fitting......
		blkNum<-max(object$gene[,"Block"])
		arrys<-unique(object$target[,"Array"])
		arryNm<-length(arrys)
		if(class(x)=="numeric")
		{
			x<-matrix(rep(x, arryNm),nrow=arryNm, ncol=length(x),byrow=T)
			#x<-
		}
		#check for the dim
		if(dim(x)[1]<arryNm)
		{
			stop("x is not in correct form, please specify!")
		}
		cat("Start doing the fitting ......\n")
		for(i in 1:blnNum) #by blocks first
		{
			cat("block ", i,"/",blnNum, "...\n")
			for(j in 1:arrys) #by array 
			{
				#now take out the array data from the object
				aNm<-arrys[j]
				idx.aNm<-which(object$target[,"Array"]==aNm)
				#now we get the index of the one array j, now need to do fitting by block
				
				#for the array, we need to get the protein data by block
				idx.prE<-which(is.element(object$gene[,"Block"],c((i+4):(i+7)))) #which(object$gene[,"Block"]==i|object$gene[,"Block"]==i+1|object$gene[,"Block"]==i+2|object$gene[,"Block"]==i+3)
				idx.prC<-which(is.element(object$cgene[,"Block"],c((i+4):(i+7))))  #which(object$cgene[,"Block"]==i|object$cgene[,"Block"]==i+1|object$cgene[,"Block"]==i+2|object$cgene[,"Block"]==i+3)
				
				#now get idx, need to get the data into matrix
				y<-rbind(object$E[idx.prE, idx.aNm], object$C[idx.prC,idx.aNm])
				x.arry<-x[j,]
				
				#now doing fitting
				f<-gainAdjust.fit5PL(y, x.arry,F, aggregated=F, debug=T)
				
				plot(c(250,650),c(0.05, 65553),log="xy",type="n")
				for(k in 1:length(y[,1]))
				{
					lines(x.arry,(y[k,]), col=k)
				}
			}
		}
	
	}
	
	###function to get rid of small flot data set
	#'@export
	rmLows<-function(mtx, threshold=40, aggregated=F, index)
	{
		if(missing(mtx))
		{
			stop("please specify the input data matrix")
		}
		if(missing(index))
		{
			step("please specify the input index of data matrix to compare with")
		}
		
		idx<-which(mtx[,index]<threshold)
		if(!aggregated)  #need to figure out the indices
		{
			idx<-unique(ceiling(idx/2))
			idx<-rep(idx,rep(2,length(idx)))
			idx<-idx*2
			dif<-c(1,0)
			dif<-rep(dif, length(idx)/2)
			idx<-idx-dif
		}
		if(length(idx)>0)
		{
			mtx[-(idx),]
		}
		else{
			mtx
		}
	}
	
#this is the outermost driving function to call the functions to do the fitting and adjustment
#'@export
PMT.adjust<-function(data, #raw date after reading the GPR files
			x, #vector holding the PMT gains for reading the chips, the arrangement of data in PMT reading is identical to the element order in x.
			data.ylog=F, data.aggregated=F,#description of data, usually, don't need to change, for advance user only
			fit.mode="control", ##fit the 5pl with control data only, to speed up 
			fit.aggregated=F, #fit the 5pl with aggregated data, for advanced data
			aggregated.mode="geometric", #controlling the way to do aggregation, geometric is preferred
			block.size=12, #size of blocks to do fitting, the more the better??. but too slow, 12 is appropriate
			a=1.1, d=65535, #change only when you know what you do.
			debug=T  #showing the debugging information
			)
{
	#check data
	if(missing(data)||missing(x))
	{
		stop("please specify the data/x.")
	}
	if(class(data)!="EListRaw")
	{
		stop("the data object is not in correct format!!")
	}
	if(class(x)!="numeric")
	{
		stop("the x is not in correct format!!")
	}
	if(floor(dim(data$E)[2]/length(x))*length(x)!=dim(data$E)[2])
	{
		stop("data and x are not correctly set up, please check")
	}
	#now start picking the data 
	blkNum<-max(data$gene[,"Block"])
	arrys<-unique(data$target[,"Array"])
	arryNm<-length(arrys)

	#i<-j<-1
	#cat("block ", i,"/",blkNum, "...\n")
	y<-list("C"=NULL, "E"=NULL, "cgene"=NULL, "gene"=NULL)
	yAdj<-list("E"=matrix(0,ncol=arryNm, nrow=dim(data$E)[1],byrow=T), 
				"C"=matrix(0,ncol=arryNm, nrow=dim(data$C)[1],byrow=T), 
				"E.adj"=matrix(0,ncol=arryNm, nrow=dim(data$E)[1],byrow=T), 
				"C.adj"=matrix(0,ncol=arryNm, nrow=dim(data$C)[1],byrow=T),
				"E.PMT"=NULL, "C.PMT"=NULL
			)

	PMT.gain<-NULL;
	for( j in 1:arryNm)
	{
		#now take out the array data from the object
		aNm<-arrys[j]
	
		idx.aNm<-which(data$target[,"Array"]==aNm)
		#now we get the index of the one array j, now need to do fitting by block
		
		#now get idx, need to get the data into matrix
		
		y$C<-object$C[,idx.aNm] ; #rbind(object$E[idx.prE, idx.aNm], object$C[idx.prC,idx.aNm])
		y$E<-object$E[,idx.aNm];
		y$cgene<-object$cgene;
		y$gene<-object$gene;
		
		#if(j==1)
		#{
		#	
		#}
		rs<-gainAdjust.fit5plArray(y, x, data.ylog=data.ylog, data.aggregated=data.aggregated,
		fit.mode=fit.mode,
		fit.aggregated=fit.aggregated, aggregated.mode=aggregated.mode, order=1,
		block.size=block.size,a=a, d=d,threshold=60,
		PMT.gain=PMT.gain,
		debug=debug
		)
		#need to get PMT.gain using the first one as reference
		PMT.gain<-list("E"=rs$E.gain, "C"=rs$C.gain);
		
		#at end of each run, need to add the result to the output list
		yAdj$E[,j]<-rs$E
		yAdj$C[,j]<-rs$C
		yAdj$E.adj[,j]<-rs$E.adj
		yAdj$C.adj[,j]<-rs$C.adj
	}
	yAdj$E.PMT<-PMT.gain$E
	yAdj$C.PMT<-PMT.gain$C
	yAdj;
}