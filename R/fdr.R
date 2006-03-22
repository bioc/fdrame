###############################################################################################
#											      #
#			                  FDR MAIN FUNCTION				      #
#											      #
###############################################################################################

fdr.ma<- function(exp.arr=NA,design=NA,p.method="resampling",fdr.adj="BH-LSU",equal.var=TRUE,plot=c("pvlVSrank","adjVSstat"),perms.num=100)
{
	if (is.na(exp.arr[1])) stop("It's impossible to work without data. You must specify a valid \"exp.arr\" .")
	if (is.na(design[1])) stop("You must specify a valid \"design\" vector.")

	if (length(unique(design))==2)				# Decides whether to run f test or t test according to the number of groups in the design
	{
		if ((equal.var==FALSE)&&(p.method=="theoretic"))
			test<-"t.welch"
		else if ((equal.var==FALSE)&&(p.method=="resampling"))
			stop("Resampling p.method may be use only under equal variance assumption")
		else test<-"t.equal.var"		
	}
	else if (length(unique(design))>2)			
		test<-"f"
	else	
		stop("There is only one group in the design vector!. Please fix it and run the process again.") 
	 
	if (p.method=="resampling")				# Calling core functions that uses resampling
	{
		if (perms.num<=1)	
			stop("Number of pemutation must be >=2.")
		if (fdr.adj=="BH-LSU")
			fdr.output<-fdr.bh(exp.arr,design,perms.num=perms.num,ref.vector="HYBRID",test=test)	
		else if (fdr.adj=="point.est")
			fdr.output<-fdr.q(exp.arr,design,perms.num=perms.num,ref.vector="HYBRID",test=test)
		else if (fdr.adj=="upper.est")
			fdr.output<-fdr.qu(exp.arr,design,perms.num=perms.num,ref.vector="HYBRID",test=test)
		else if (fdr.adj=="adaptive")
			fdr.output<-fdr.adaptive.c.resampling(exp.arr,design,perms.num=perms.num,ref.vector="NULL",test=test)
		else 
		stop(paste("There is no such option!: p.method=",p.method," ; fdr.adj=", fdr.adj))

	}
	else if (p.method=="theoretic")				# calling core functions that doesn't uses resampling
	{
		perms.num<-1
		if (fdr.adj=="BH-LSU")
			fdr.output<-fdr.pt(exp.arr,design,ref.vector="HYBRID",test=test)
		else if (fdr.adj=="adaptive")
			fdr.output<-fdr.adaptive.c(exp.arr,design,ref.vector="NULL",test=test)
		else 
			stop(paste("ERROR: There is no such option!: p.method=",p.method," ; fdr.adj=", fdr.adj))

	}
	else 	stop(paste("ERROR: There is no such method!: p.method=",p.method))

	xy<-fdr.plot(info.list=fdr.output,plot=plot,test=test)	# calling plot to draw graphs and compute approximated values
	
	return(list(adj=xy$y,dif=fdr.output$dif))
}

###############################################################################################
#											      #
#			       FDR COMMON COMPUTATIONS FUNCTION		 	              #
#											      #
###############################################################################################

fdr.basic.comp<- function(exp.arr,design,test="t.welch",ref.vector="NULL",jpoint=2.5,perms.num=1) # Inits and builds rejections vector
{
	#Each one of the core functions call this function. Here the statistic is computed, refference vector is created, pvalues and resampled pvalues are computed.
	if (test=="t.welch")
		ref.vector<-"NULL"
	n.total<-length(design)
	n.groups<-length(unique(design))
	if (length(design)>dim(exp.arr)[2])
		stop(paste("ERROR: The number of columns in the file(=",dim(exp.arr)[2],") is smaller than the `design` size (=",length(design),")"))

	genes.num<-dim(exp.arr)[1]			#each row represents a gene 
	cs<-compute.statistic(exp.arr,design,test)
	statistic.vector<-cs$statistic.vector	#compute stat each gene (No Resampling)
	
	if ((test=="t.welch")|(test=="t.equal.var"))	jpoint<-2
	else if (test=="f") jpoint<- qf(0.95,n.groups-1,dim(exp.arr)[2]-n.groups)

	if (ref.vector[1]=="NULL")		 				#checks whether it has ref.vector
	{									# if not, it uses the real T values as refference
		sorted.statistic.vector<-sort(abs(statistic.vector), decreasing = TRUE)
		ref.vector<-sorted.statistic.vector
		ref.vector.real.values<-rep("TRUE",length(ref.vector))
	}
	else if (ref.vector[1]=="HYBRID")
	{
		sorted.statistic.vector<-sort(abs(statistic.vector), decreasing = TRUE)
		r1<-sorted.statistic.vector[sorted.statistic.vector>jpoint]
		r2<-seq(jpoint,0,length=200)
		ref.vector<-as.numeric(c(r1,r2))
		ref.vector.real.values<-c(rep("TRUE",length(r1)),rep("FALSE",length(r2)))
	}
	else if (is.numeric(ref.vector))
	{
		ref.vector.real.values<-rep("FALSE",length(ref.vector))
	}		
		
	
	ref.vector.size<-length(ref.vector)
	r.vector<-vector("numeric",ref.vector.size)

	for (i in 1:ref.vector.size)	
	{
		r.vector[i]<-sum(abs(statistic.vector)>ref.vector[i],na.rm=TRUE)	#r.vector[m] counts the number rejected values (bigger than ref.vector[m])
	}
	
	



	if ((test=="t.welch")|(test=="t.equal.var"))
		pvalues<-2*(1-pt(as.numeric(ref.vector),cs$df))	#compute pvalues from t-statistic
	#	
	#	pvalues<-2*(1-pt(as.numeric(ref.vector),(n.total-2)))	#compute pvalues from t-statistic
	
	else if (test=="f")
		pvalues<-(1-pf(as.numeric((ref.vector)),n.groups-1,(n.total-n.groups)))	#compute pvalues from f-statistic
	
	ud<-unique(design)
	groups.sizes<-vector("numeric",length(ud))
	for (i in 1:length(ud))
		groups.sizes[i]<-sum(design==ud[i])
	
	if (perms.num>1)
	{
		resamp.pvalues<-vector("numeric",ref.vector.size)
		res.stat.matrix<-get.resampling.statistic.array(exp.arr,design,perms.num,groups.sizes,test)
		for (i in 1:ref.vector.size)		
			resamp.pvalues[i]<-sum(abs(res.stat.matrix)> ref.vector[i])/(perms.num*genes.num)
		return(list(genes.num=genes.num,design=design,ref.vector=as.vector(ref.vector),pvalues=pvalues,r.vector=r.vector,statistic.vector=statistic.vector,ref.vector.size=ref.vector.size,ref.vector.real.values=ref.vector.real.values,dif=cs$dif,groups.sizes=groups.sizes,resamp.pvalues=resamp.pvalues,res.stat.matrix=res.stat.matrix))
	}
	else	return(list(genes.num=genes.num,design=design,ref.vector=as.vector(ref.vector),pvalues=pvalues,r.vector=r.vector,statistic.vector=statistic.vector,ref.vector.size=ref.vector.size,ref.vector.real.values=ref.vector.real.values,dif=cs$dif,groups.sizes=groups.sizes,test=test))
	
}



###############################################################################################

fdr.plot <-function(info.list,plot=c("pvlVSrank","adjVSrank"),test=test)
{
	px.values<-NA
	py.values<-NA
	
	if ((info.list$p.method=="theoretic")&&(info.list$fdr.adj=="BH-LSU"))
		y.values<-info.list$pt.adjust
	else if ((info.list$p.method=="theoretic")&&(info.list$fdr.adj=="adaptive"))
		y.values<-info.list$adapk
	else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="BH-LSU"))
		y.values<-info.list$bh.value
	else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="upper.est"))
		y.values<-info.list$qu.value
	else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="point.est"))
		y.values<-info.list$q.value

	

	x.values<-info.list$ref.vector
	

	x.leg<-paste("observed |",test,"|")
	
	
	
	fdr.pa<-approx(spline(x.values,y.values),xout=abs(info.list$statistic.vector))
	
	px.values<-fdr.pa$x
	py.values<-fdr.pa$y
	py.values<-ifelse(py.values>1,1,py.values)

	outx<-px.values
	outy<-py.values
	
	if (length(plot)>0)		#adjusted
	{
		
	
		if (sum(plot==rep("adjVSstat",length(plot)))>0)
		{
			if ((info.list$p.method=="theoretic")&&(info.list$fdr.adj=="BH-LSU"))
				leg<-"BH LSU - No Resampling"
			else if ((info.list$p.method=="theoretic")&&(info.list$fdr.adj=="adaptive"))
				leg<-"Adaptive - No Resampling"
			else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="BH-LSU"))
				leg<-"BH point-estimate"
			else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="upper.est"))
				leg<-"Res. upper limit"
			else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="point.est"))
				leg<-"Res. point-estimate"
			
			x11()
			par(mfcol=c(1,1))
				
			plot(px.values,py.values,col=1,"p",ylab=paste("Adjusted P-Values by ",leg),xlab=x.leg,main="FDR Adjusted P-Values VS Observed Statistic")
			abline(0.25,0,lwd=0.01,lty=2) 
			abline(0.2,0,lwd=0.01,lty=2)
			abline(0.15,0,lwd=0.01,lty=2)
			abline(0.1,0,lwd=0.01,lty=2)
			abline(0.05,0,lwd=0.01,lty=2)

		}
			
		if (sum(plot==rep("pvlVSrank",length(plot)))>0)
		{
			x.leg<-"rank of P-values"

			x11()
			par(mfcol=c(1,1))

			if ((info.list$p.method=="resampling"))
				#fdr.pa<-approx(spline(info.list$ref.vector,info.list$resamp.pvalues),xout=abs(info.list$statistic.vector))			
				fdr.pa$y<-info.list$resamp.pvalues
				else
				fdr.pa<-approx(spline(info.list$ref.vector,info.list$pvalues),xout=abs(info.list$statistic.vector))
			
			py.values<-fdr.pa$y
			py.values<-ifelse(py.values>1,1,py.values)
			
			leg<-paste("P-values based on ",test," test")

			px.values<-rank(py.values)
			plot(px.values,py.values,col=1,"p",ylab=leg,xlab=x.leg,main="P-Values VS Rank")
			abline(0.25,0,lwd=0.01,lty=2) 
			abline(0.2,0,lwd=0.01,lty=2)
			abline(0.15,0,lwd=0.01,lty=2)
			abline(0.1,0,lwd=0.01,lty=2)
			abline(0.05,0,lwd=0.01,lty=2)

		}

		if (sum(plot==rep("adjVSrank",length(plot)))>0)
		{
			if ((info.list$p.method=="theoretic")&&(info.list$fdr.adj=="BH-LSU"))
			{
				y.values<-info.list$pt.adjust
				leg<-"BH LSU - No Resampling"
			}
			else if ((info.list$p.method=="theoretic")&&(info.list$fdr.adj=="adaptive"))
			{
				y.values<-info.list$adapk
				leg<-"Adaptive - No Resampling"
			}
			else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="BH-LSU"))
			{
				y.values<-info.list$bh.value
				leg<-"BH point-estimate"
			}
			else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="upper.est"))
			{
				y.values<-info.list$qu.value
				leg<-"Res. upper limit"
			}
			else if ((info.list$p.method=="resampling")&&(info.list$fdr.adj=="point.est"))
			{
				y.values<-info.list$q.value
				leg<-"Res. point-estimate"
			}

			
			x.leg<-"rank of adjusted P-values"

			x11()
			par(mfcol=c(1,1))

			fdr.pa<-approx(spline(info.list$ref.vector,y.values),xout=abs(info.list$statistic.vector))
			py.values<-fdr.pa$y
			py.values<-ifelse(py.values>1,1,py.values)
			
			if ((info.list$p.method=="resampling"))
				fdr.pa$y<-info.list$resamp.pvalues
				else
				fdr.pa<-approx(spline(info.list$ref.vector,info.list$pvalues),xout=abs(info.list$statistic.vector))
			
			px.values<-fdr.pa$y
			px.values<-ifelse(px.values>1,1,px.values)


			
			#if ((info.list$p.method=="resampling"))
			#	y.values<-info.list$resamp.pvalues
			#else if  ((info.list$p.method=="theoretic"))
			#	y.values<-info.list$pvalues
			
		#	fdr.pa<-approx(spline(x.values,y.values),xout=abs(info.list$statistic.vector))			
		#	px.values<-rank(fdr.pa$y)		

			plot(rank(px.values),py.values,col=1,"p",ylab=paste("Adjusted P-Values by ",leg),xlab=x.leg,main="FDR Adjusted P-Values VS Rank")
			abline(0.25,0,lwd=0.01,lty=2) 
			abline(0.2,0,lwd=0.01,lty=2)
			abline(0.15,0,lwd=0.01,lty=2)
			abline(0.1,0,lwd=0.01,lty=2)
			abline(0.05,0,lwd=0.01,lty=2)

		}
	
	}	
	return(list(x=outx,y=outy))
}

read.matrix <-function(file)
{
	return (as.matrix(read.table(file,TRUE,",")))
}




###################################################################################################


gene.expression.normalization <- function(genes,genes.num,n.total)	#	Normalization
{		
	#	Definitions for Normalization
	
	a.exp.array<-matrix(NA,genes.num,n.total)
	m.exp.array<-matrix(NA,genes.num,n.total)
	exp.arr<-matrix(NA,genes.num,n.total)
	lowess.m<-vector("numeric",genes.num)
	for (j in 1:(n.total))	
	{	
		k<-2*j
		a.exp.array[,j]<-ifelse((genes[,k-1]==0 | genes[,k]==0),0,log(genes[,k],base=2)+log(genes[,k-1],base=2))
		m.exp.array[,j]<-ifelse((genes[,k-1]==0 | genes[,k]==0),0,log(genes[,k],base=2)-log(genes[,k-1],base=2))
							   
		lowess.m<-lowess(a.exp.array[,j],m.exp.array[,j],f=1/5)$y
		uns<-match(rank(a.exp.array[,j]), unique(sort(rank(a.exp.array[,j]))))
		exp.arr[,j]<-m.exp.array[,j]-lowess.m[uns]
	}
	return(exp.arr)
}

##############################################################################################
#											     #
#				Statistic computation functions				     #
#										             #
##############################################################################################
#	ComputeS t-statistic
#	(in this case, t for comparing two independent samples)

compute.t.statistic<-function(exp.arr,design,equal.var=TRUE)
{
	columns.b<-(1:length(design))[unique(design)[1]==design]
	columns.a<-(1:length(design))[unique(design)[2]==design]
	
	n.b<-length(columns.b)
	n.a<-length(columns.a)
	rep.one.n.b<-rep.int(1,n.b)
	rep.one.n.a<-rep.int(1,n.a)

	mean.b<-((exp.arr[,(columns.b)]%*%rep.one.n.b)/n.b)
	mean.a<-((exp.arr[,(columns.a)]%*%rep.one.n.a)/n.a)
	var.b<-(as.vector((((exp.arr[,columns.b]-
			as.vector((exp.arr[,columns.b]%*%rep.one.n.b)/n.b))^2)/(n.b-1))%*%rep.one.n.b)) 
	var.a<-(as.vector((((exp.arr[,columns.a]-
			as.vector((exp.arr[,columns.a]%*%rep.one.n.a)/n.a))^2)/(n.a-1))%*%
			rep.one.n.a))
#	mean.b<-apply(exp.arr[,(columns.b)],1,mean)
#	mean.a<-apply(exp.arr[,(columns.a)],1,mean)
#	var.b<-apply(exp.arr[,(columns.b)],1,var)
#	var.a<-apply(exp.arr[,(columns.a)],1,var)

	if (equal.var)
	{
		s.pooled<-(var.a*(n.a-1)+var.b*(n.b-1))/(n.a+n.b-2)
		statistic.vector<-as.vector(mean.b-mean.a)/sqrt(s.pooled*((1/n.b)+(1/n.a)))
	}
	else
		statistic.vector<-as.vector(mean.b-mean.a)/sqrt((var.b/n.b)+(var.a/n.a))

	exp.nas.rows.nums<-(1:dim(exp.arr)[1])[is.na(statistic.vector)]
	if (length(exp.nas.rows.nums)>0)
	{
		exp.nas<-exp.arr[exp.nas.rows.nums,]
		for (i in 1:length(exp.nas.rows.nums))
		{
			columns.b<-columns.b[is.na(exp.nas[i,columns.b])==FALSE]
			columns.a<-columns.a[is.na(exp.nas[i,columns.a])==FALSE]
		
			n.b<-length(columns.b)
			n.a<-length(columns.a)
			mean.b<-mean(exp.nas[i,(columns.b)])
			mean.a<-mean(exp.nas[i,(columns.a)])
			rep.one.n.b<-rep.int(1,n.b)
			rep.one.n.a<-rep.int(1,n.a)
			var.b<-(((((exp.nas[i,columns.b]-
				((exp.nas[i,columns.b]%*%rep.one.n.b)/n.b))^2)/(n.b-1))%*%rep.one.n.b)) 
			var.a<-(((((exp.nas[i,columns.a]-
				((exp.nas[i,columns.a]%*%rep.one.n.a)/n.a))^2)/(n.a-1))%*%rep.one.n.a))
			if (equal.var)
			{
				s.pooled<-(var.a*(n.a-1)+var.b*(n.b-1))/(n.a+n.b-2)
				statistic.vector[exp.nas.rows.nums[i]]<-(mean.b-mean.a)/sqrt(s.pooled*((1/n.b)+(1/n.a)))
			}
			else
				statistic.vector[exp.nas.rows.nums[i]]<-(mean.b-mean.a)/sqrt((var.b/n.b)+(var.a/n.a))

			
		}
	}
	
	if (!equal.var)
	{
		#print(var.b)
		#print(var.a)
		#print(n.a)
		c<-(var.a/n.a)/(var.a/n.a+var.b/n.b)
		
		df<-1/(c^2/(n.a-1)+(1-c)^2/(n.b-1))
		#print(c)
		#print((1-c^2)/(n.b-1))
		
	}
	else 
		df<-n.a+n.b-2
	return(list(statistic.vector=statistic.vector,dif=as.vector(mean.b-mean.a),df=df))
}


compute.f.statistic<-function(exp.arr,design)
{
	
	genes.num<-dim(exp.arr)[1]
	dif<-vector("numeric",genes.num)
	statistic.vector<-vector("numeric",genes.num)

	#for (i in 1:genes.num)
	#	statistic.vector[i]<-summary(aov(exp.arr[i,]~design))[[1]]$F[1]
	for (i in 1:genes.num)
	{
		summarylm<-summary(lm(exp.arr[i,]~design))
		statistic.vector[i]<-summarylm[10]$fstatistic[1]
		dif[i]<-summarylm[8]
	}
	
	return(list(statistic.vector=statistic.vector))
}

compute.statistic<-function(exp.arr,design,test)
{
	if (test=="t.welch") 
		return(compute.t.statistic(exp.arr,design,equal.var=FALSE))
	if (test=="t.equal.var")
		return(compute.t.statistic(exp.arr,design,equal.var=TRUE))
	else if (test=="f") return(compute.f.statistic(exp.arr,design))
	else 
	{
		stop("ERROR: It's impossible to run this test. Only t,f tests are implemented.")
	}
}


###############################################################################################
get.resampling.statistic.array<-function(exp.arr,design,perms.num,groups.sizes,test="t.equal.var")
{
	row.doesnt.contain.na<-apply(is.na(exp.arr),1,sum)==0
	genes.num<-sum(row.doesnt.contain.na)
	if (test=="t.equal.var")
	{
		statistic.array<-vector("numeric",genes.num*perms.num)
		statistic.array<-.C("compute_resampling_t_stat",as.double(as.vector(exp.arr[row.doesnt.contain.na])),as.integer(groups.sizes[1]),as.integer(groups.sizes[2]),as.integer(genes.num), as.integer(perms.num),as.double(statistic.array),PACKAGE="fdrame")[[6]]	
		statistic.matrix<-matrix(statistic.array,genes.num,perms.num)
	}
	else if (test=="f")
	{
		statistic.matrix<-compute.resampling.stat(exp.arr=exp.arr[(1:genes.num)[(row.doesnt.contain.na)],],design=design,genes.num=genes.num,perms.num=perms.num,test=test)
	}
	else
	{
		stop("test must be t.equal.var or f")
	}
		
	return(statistic.matrix)
}

compute.resampling.stat<- function(exp.arr,design,genes.num,perms.num,test="t.welch")
{
	statistic.matrix<-matrix(NA,genes.num,perms.num)
	
	
	n.total<-length(design)
	
	for (j in 1:perms.num)	#	Creates Sets of Permutations of the Original Data Set
	{
		
		shuffled.design<-design[sample(n.total,n.total,replace=FALSE)]	#shuffling the columns' numbers
		
		statistic.matrix[,j]<-compute.statistic(exp.arr=exp.arr,design=shuffled.design,test=test)$statistic.vector
		
	}
	return(statistic.matrix)
}


###############################################################################################
#											      #
#					CORE FUNCTIONS					      #
#											      #
###############################################################################################


fdr.pt<- function(exp.arr,design,ref.vector="NULL",test="t.welch")
{
	info.list<-fdr.basic.comp(exp.arr,design,test,ref.vector,perms.num=1) 
	pt.adjust<-info.list$pvalues*info.list$genes.num/info.list$r.vector	#BH Linear Step-Up adjusted p-values Vector - No Resampling
	pt.adjust<-ifelse(is.nan(pt.adjust),0,pt.adjust)
	pt.adjust<-rev(cummin(rev(pt.adjust)))	#Monotone Adjustment

	return(list(ref.vector=info.list$ref.vector,pt.adjust=pt.adjust,pvalues=info.list$pvalues,statistic.vector=info.list$statistic.vector,ref.vector.real.values=info.list$ref.vector.real.values,dif=info.list$dif,p.method="theoretic",fdr.adj="BH-LSU",test=test))
}



fdr.bh<- function(exp.arr,design,perms.num=100,ref.vector="NULL",test="t.equal.var")
{
	info.list<-fdr.basic.comp(exp.arr=exp.arr,design=design,test=test,ref.vector=ref.vector,perms.num=perms.num)
	
	m<-info.list$ref.vector.size
	
	pa<-approx(spline(info.list$ref.vector,info.list$resamp.pvalues),xout=abs(info.list$statistic.vector))	

	bh.value<-ifelse(info.list$r.vector>0,(info.list$resamp.pvalues*info.list$genes.num/(info.list$r.vector)),0) 
	bh.value<-rev(cummin(rev(bh.value)))	#Monotone Adjustment
	
	return(list(ref.vector=info.list$ref.vector,bh.value=bh.value,pvalues=info.list$pvalues,statistic.vector=info.list$statistic.vector,ref.vector.real.values=info.list$ref.vector.real.values,dif=info.list$dif,p.method="resampling",fdr.adj="BH-LSU",resamp.pvalues=pa$y,test=test))
}

fdr.qu<- function(exp.arr,design,perms.num=100,ref.vector="NULL",test="t.equal.var",alpha=0.05)  #upper.est
{
	info.list<-fdr.basic.comp(exp.arr=exp.arr,design=design,test=test,ref.vector=ref.vector,perms.num=perms.num)
	
	m<-info.list$ref.vector.size
	qu.value<-vector("numeric",info.list$ref.vector.size)	
	est.95qu.exc<-vector("numeric",info.list$ref.vector.size)	
	su.hat<-vector("numeric",info.list$ref.vector.size)	
	genes.num<-info.list$genes.num
	
	for (i in 1:info.list$ref.vector.size)
	{
		v<-apply(abs(info.list$res.stat.matrix) > info.list$ref.vector[i],2,sum)	#count the number of rejected genes for each permutation
		est.95qu.exc[i]<- quantile(v,prob=(1-alpha),,names = FALSE)	#compute qu for these values
		su.hat[i]<-info.list$r.vector[i]-est.95qu.exc[i]
		su.hat[i]<- ifelse(su.hat[i]<0,0,su.hat[i])			#insure values are not negative
		v<-ifelse((v+su.hat[i])!=0,v/(v+su.hat[i]),0)     # checks division by zero
		qu.value[i]<-mean(v)			#Resampling Upper Limit Estimate adjusted p-values Vector
	}
	qu.value<-cummax(qu.value)
	pa<-approx(spline(info.list$ref.vector,info.list$resamp.pvalues),xout=abs(info.list$statistic.vector))	
	
	return(list(ref.vector=info.list$ref.vector,qu.value=qu.value,pvalues=info.list$pvalues,statistic.vector=info.list$statistic.vector,ref.vector.real.values=info.list$ref.vector.real.values,dif=info.list$dif,p.method="resampling",fdr.adj="upper.est",resamp.pvalues=pa$y,test=test))	
}						                     

fdr.q<- function(exp.arr,design,perms.num=100,ref.vector="NULL",test="t.equal.var",alpha=0.05)
{
	info.list<-fdr.basic.comp(exp.arr=exp.arr,design=design,test=test,ref.vector=ref.vector,perms.num=perms.num)
	
	m<-info.list$ref.vector.size
	resamp.pvalues<-vector("numeric",m)
	est.95qu.exc<-vector("numeric",info.list$ref.vector.size)	
	s.hat<-vector("numeric",info.list$ref.vector.size)
	q.value<-vector("numeric",info.list$ref.vector.size)
	genes.num<-info.list$genes.num

	for (i in 1:info.list$ref.vector.size)	
	{
		v<-apply(abs(info.list$res.stat.matrix) > info.list$ref.vector[i],2,sum)
		est.95qu.exc[i]<- quantile(v,prob=(1-alpha),,names = FALSE)
		s.hat[i]<- info.list$r.vector[i]-mean(v)
		s.hat[i]<- ifelse(s.hat[i]<est.95qu.exc[i],0,s.hat[i])
		v<-ifelse((v+s.hat[i])!=0,v/(v+s.hat[i]),0)     # checks division by zero
	
		q.value[i]<-mean(v)		#Resampling Point Estimate adjusted p-values Vector
		
	}
	
	q.value<-rev(cummin(rev(q.value)))	#Monotone Adjustment
	pa<-approx(spline(info.list$ref.vector,info.list$resamp.pvalues),xout=abs(info.list$statistic.vector))	
		
	return(list(ref.vector=info.list$ref.vector,q.value=q.value,pvalues=info.list$pvalues,statistic.vector=info.list$statistic.vector,ref.vector.real.values=info.list$ref.vector.real.values,dif=info.list$dif,p.method="resampling",fdr.adj="point.est",resamp.pvalues=pa$y,test=test))	
} 


fdr.adaptive.c<- function(exp.arr,design,ref.vector="NULL",test="t.welch")
{	
	info.list<-fdr.basic.comp(exp.arr=exp.arr,design=design,test=test,ref.vector=ref.vector,perms.num=1)
	
	m<-info.list$ref.vector.size
	adapk<-vector("numeric",m)	
	adapk<-.C("adaptive",as.double(info.list$pvalues),as.integer(m),as.double(adapk),PACKAGE="fdrame")[[3]]		
	
	return(list(ref.vector=info.list$ref.vector,adapk=adapk,pvalues=info.list$pvalues,statistic.vector=info.list$statistic.vector,ref.vector.real.values=info.list$ref.vector.real.values,dif=info.list$dif,p.method="theoretic",fdr.adj="adaptive",test=test))		
}

fdr.adaptive.c.resampling<- function(exp.arr,design,perms.num=100,ref.vector="NULL",test="t.equal.var")
{
	info.list<-fdr.basic.comp(exp.arr=exp.arr,design=design,test=test,ref.vector=ref.vector,perms.num=perms.num)
	m<-info.list$ref.vector.size
	adapk<-vector("numeric",m)	
	statistic.array<-vector("numeric",info.list$genes.num*perms.num)
			
	pa<-approx(spline(info.list$ref.vector,info.list$resamp.pvalues),xout=abs(info.list$statistic.vector))	
	adapk<-.C("adaptive",as.double(info.list$resamp.pvalues),as.integer(m),as.double(adapk),PACKAGE="fdrame")[[3]]		
	
	return(list(ref.vector=info.list$ref.vector,adapk=adapk,pvalues=info.list$pvalues,statistic.vector=info.list$statistic.vector,ref.vector.real.values=info.list$ref.vector.real.values,dif=info.list$dif,p.method="theoretic",fdr.adj="adaptive",resamp.pvalues=pa$y,test=test))		
}
