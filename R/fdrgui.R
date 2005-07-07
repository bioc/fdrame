OnRbpm1 <-function()
{
	tkconfigure(entry.Perms, state="disabled")

}

OnRbpm2 <-function()
{
	tkconfigure(entry.Perms, state="normal")

}

OnRbfa1 <-function()
{
	#tkconfigure(rbpm1, state="normal")
	#tkconfigure(rbpm2, state="normal")
	#tkconfigure(entry.Perms, state="normal")
	
	tkconfigure(cbfa1, state="normal")
	tkconfigure(cbfa2, state="disabled")
}

OnRbfa2 <-function()
{
#	tkconfigure(rbpm1, state="normal")
#	tkconfigure(rbpm2, state="normal")
#	tkconfigure(entry.Perms, state="normal")
	tkconfigure(cbfa1, state="disabled")
	tkconfigure(cbfa2, state="normal")

}

OnRbfa3 <-function()
{
#	rbpmValue <<- tclVar("resampling")
#	tkconfigure(rbpm1,variable=rbpmValue,value="theoretic")
#	tkconfigure(rbpm2,variable=rbpmValue,value="resampling")
#	tkconfigure(rbpm2,variable=rbpmValue,value="resampling")
#	tkconfigure(rbpm1, state="disabled")
#	tkconfigure(rbpm2, state="normal")
	
#	tkconfigure(entry.Perms, state="normal")

	tkconfigure(cbfa1, state="disabled")
	tkconfigure(cbfa2, state="disabled")
	
	
}

OnRbfa4 <-function()
{
#	rbpmValue <<- tclVar("resampling")
#	tkconfigure(rbpm1,variable=rbpmValue,value="theoretic")
#	tkconfigure(rbpm2,variable=rbpmValue,value="resampling")
	
#	tkconfigure(rbpm1, state="disabled")
#	tkconfigure(rbpm2, state="normal")
#	tkconfigure(entry.Perms, state="normal")

	tkconfigure(cbfa1, state="disabled")
	tkconfigure(cbfa2, state="disabled")

}


OnBrowseInput <-function()
{
	filename.input.tk<<-tkgetOpenFile()
	if(tclvalue(filename.input.tk)!="")
	{	
		filename.input<<-tclVar(tclvalue(filename.input.tk))
		tkconfigure(entry.Filename.input,textvariable=filename.input)
	}
}


OnOK <- function()
{
	
	design.string <- vector.list[as.numeric(tkcurselection(tl))+1]
	tkdestroy(ttdesignlist)
	
	designVal<<-eval(parse(text=design.string))
	design.tcl<<-tclVar(design.string)
	tkconfigure(entry.Design,textvariable=design.tcl)
	tkfocus(fdrtt)
	designFlag<<-"FALSE"

}

OnCancel <- function()
{
	
	tkdestroy(ttdesignlist)
	tkfocus(fdrtt)
	designFlag<<-"FALSE"
}

OnBrowseDesign <-function()
{
	if (designFlag=="TRUE") return()
	else designFlag<<-"TRUE"
	object.list<-ls(envir=globalenv())
	vector.list<<-c();
	for (i in 1:length(object.list))
	{
		if ((is.vector(eval(parse(text=object.list[i])),"character"))||(is.vector(eval(parse(text=object.list[i])),"numeric")))
			vector.list<<-c(vector.list,object.list[i])
	}	
	if (length(vector.list)>0)
	{
		ttdesignlist<<-tktoplevel()
		tl<<-tklistbox(ttdesignlist,height=length(vector.list),selectmode="single",background="white")
		tkgrid(tklabel(ttdesignlist,text="Please choose a vector"))
		tkgrid(tl,columnspan=2)
		for (i in (1:length(vector.list)))
		{
    			tkinsert(tl,"end",vector.list[i])
		}
		OK.but <-tkbutton(ttdesignlist,text="   OK   ",command=OnOK)
		Cancel.but <-tkbutton(ttdesignlist,text=" Cancel ",command=OnCancel)
		tkgrid(OK.but,Cancel.but)
		tkfocus(ttdesignlist)
	}
	else print("There are no vectors in the memory")
	tkfocus(fdrtt)
}


OnBrowseOutput <-function()
{
	filename.output.tk<<-tkgetSaveFile()
	if(tclvalue(filename.output.tk)!="")
	{
		filename.output<<-tclVar(tclvalue(filename.output.tk))
		tkconfigure(entry.Filename.output,textvariable=filename.output)
	}
}

OnHelp <- function()
{

	
	tt  <- tktoplevel()
	scr <- tkscrollbar(tt, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
	txt <- tktext(tt,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
	tkgrid(txt,scr)
	tkgrid.configure(scr,sticky="ns")
	tkinsert(txt,"end","Help page for Fdr.gui() - FDR graphic User Interface for fdr.ma()\n")
	tkinsert(txt,"end","
This progran takes unnormalized or normalized expression data array, 
written to the input file specified by the user (without a header) experimental
design and computes adjusted p-values.
An output .Rda file is written. It contains:
* adj - FDR adjusted p-values
* dif - Diffrences between mean of groups (when there are two groups),
        Multiple R-Squared values (when there are more than 2 groups).
The user should specify which type of p-value adjustment method may be used.
The options are:
* Benjamini-Hochberg Linear Step Up porcedure (with or without resampling).
* Point estimation procedure (with or without resampling).
* upper estimation procedure (with resampling).
* two-stage adaptive procedure (with resampling).
When using resampling, the user may specify the number of permutations.
There are 3 types of plots that can be drawn:
* P-values VS Rank
* Adjusted P-values VS Rank
* Adjusted P-vlaues VS Observed statistic
The program starts it's proccessing stage when Go! button is pressed. 
It may take several minuetes.
Due to TCLTK focus problems on windows it's reccomended to use the 
RGui SDI mode.
References:
* Reiner A, Yekutieli D, Benjamini Y: Identifying differentially
     expressed genes using false discovery rate controlling procedures.
     Bioinformatics 19:368-375, 2003 
*Benjamini, Y., Krieger,A.M.,Yekutieli, D. (2001) 
     ?Two Staged Linear Step Up FDR Controlling Procedure?, Technical 
     Report Department of Statistics and O.R., Tel Aviv University.
")
	
	tkconfigure(txt, state="disabled")
	tkfocus(txt)


}


OnCompute <- function()
{
	
	tkconfigure(fdrtt,cursor="watch")
	norm<-as.character(tclvalue(rbnoValue))
	
	fdr.adj <- as.character(tclvalue(rbfaValue))
	
	cbpl1Val <- as.character(tclvalue(cbpl1Value))
	cbpl2Val <- as.character(tclvalue(cbpl2Value))
	cbpl3Val <- as.character(tclvalue(cbpl3Value))
	cbfa1Val <- as.character(tclvalue(cbfa1Value))
	cbfa2Val <- as.character(tclvalue(cbfa2Value))
	cbfa3Val <- as.character(tclvalue(cbfa3Value))
	cbfa4Val <- as.character(tclvalue(cbfa4Value))

	plot<-c()
	if (cbpl1Val=="1") plot<-c(plot,"pvlVSrank")
	if (cbpl2Val=="1") plot<-c(plot,"adjVSrank")
	if (cbpl3Val=="1") plot<-c(plot,"adjVSstat")

	
	filename.inputVal<- tclvalue(filename.input)
	filename.outputVal<- tclvalue(filename.output)
	design.string<- tclvalue(design.tcl)
	ch2<-FALSE
	try({designVal<-eval(parse(text=design.string));ch2<-TRUE})
	if (ch2==FALSE) fdr.error()
	permsNumVal<- tclvalue(permsNum)
	if (norm=="normalized")
	{	
		if (file.exists(filename.inputVal))
		{	
			ch2<-FALSE
			try({exp.arr<-as.matrix(read.table(filename.inputVal));ch2<-TRUE})
			if (ch2==FALSE) fdr.error()
		}
		else
		{
			print(paste("The file \"",filename.inputVal,"\" does not exist"))
			fdr.error()
		}

		if (filename.outputVal=="")
		{
			print(paste("You must specify an output file"))
			fdr.error()
		}	
		#stop(paste("Can't read the file \"",filename.inputVal,"\". Check whether each row has the same number of columns"))
	}	
	else if (norm=="unnormalized")
	{
		input<-read.table(filename.inputVal)
		exp.arr<-gene.expression.normalization(	input,dim(input)[1],dim(input)[2]/2)			#Normalizes samples for both groups
	}
	
	if (fdr.adj=="BH-LSU")
		p.method<-ifelse(cbfa1Val=="1","resampling","theoretic")
	if (fdr.adj=="adaptive")
		p.method<-ifelse(cbfa2Val=="1","resampling","theoretic")
	if (fdr.adj=="upper.est")
		p.method<-ifelse(cbfa3Val=="1","resampling","theoretic")
	if (fdr.adj=="point.est")
		p.method<-ifelse(cbfa4Val=="1","resampling","theoretic")



	time1<-proc.time()[3]
	#a<<-tryCatch(fdr.output<<-fdr.ma(exp.arr=exp.arr,design=designVal,p.method=p.method,fdr.adj=fdr.adj,plot=plot,perms.num=as.numeric(permsNumVal)),fdr.error())
	ch2<-FALSE
	a<<-try({fdr.output<-fdr.ma(exp.arr=exp.arr,design=designVal,p.method=p.method,fdr.adj=fdr.adj,plot=plot,perms.num=as.numeric(permsNumVal));ch2<-TRUE})
	if (ch2==FALSE) fdr.error()
	time2<-proc.time()[3]
	#print(paste("TIME:",(time2-time1)," s"))
	save(fdr.output,file=filename.outputVal)
	tkconfigure(fdrtt,cursor="arrow")
	tkfocus(fdrtt)
	
}

fdr.error<- function()
{
	tkmessageBox(message="An error has occurred!",icon="error",type="ok")
	tkconfigure(fdrtt,cursor="arrow")
	stop()
}	

fdr.gui <- function()
{
	require(tcltk)
	fdrtt<<-tktoplevel()
	designFlag<<-"FALSE"

	frameOverall <<- tkframe(fdrtt)
	frameCompute <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	frameInput <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	frameOutput <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	frameFileOutput <<- tkframe(frameOutput,relief="groove",borderwidth=2)
	framePmethod <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	frameNorm <<- tkframe(frameInput,relief="groove",borderwidth=2)
	
	frameAdj <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	framePerms <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	frameplot <<- tkframe(frameOverall,relief="groove",borderwidth=2)
	

	#filename.input <<- tclVar("exp.arr.normalized")
	filename.input <<- tclVar("")
	filename.output <<- tclVar("fdr.output")
	#design.tcl <<- tclVar("c(rep(0,8),rep(1,8))")
	design.tcl <<- tclVar("")
	permsNum<<- tclVar("100")

	entry.Filename.input <<-tkentry(frameInput,width="20",textvariable=filename.input)
	entry.Filename.output <<-tkentry(frameOutput,width="20",textvariable=filename.output)
	entry.Design <<-tkentry(frameInput,width="20",textvariable=design.tcl)
	entry.Perms <<-tkentry(frameAdj,width="6",textvariable=permsNum)


	rbpm1 <<- tkradiobutton(framePmethod,command=OnRbpm1)
	rbpm2 <<- tkradiobutton(framePmethod,command=OnRbpm2)

	rbno1 <<- tkradiobutton(frameInput)
	rbno2 <<- tkradiobutton(frameInput)
	rbfa1 <<- tkradiobutton(frameAdj,command=OnRbfa1)
	rbfa2 <<- tkradiobutton(frameAdj,command=OnRbfa2)
	rbfa3 <<- tkradiobutton(frameAdj,command=OnRbfa3)
	rbfa4 <<- tkradiobutton(frameAdj,command=OnRbfa4)
	
	cbpl1 <<- tkcheckbutton(frameplot)
	cbpl2 <<- tkcheckbutton(frameplot)
	cbpl3 <<- tkcheckbutton(frameplot)
	
	cbfa1 <<- tkcheckbutton(frameAdj)
	cbfa2 <<- tkcheckbutton(frameAdj)
	cbfa3 <<- tkcheckbutton(frameAdj)
	cbfa4 <<- tkcheckbutton(frameAdj)


	rbpmValue <<- tclVar("resampling")
	rbnoValue <<- tclVar("normalized")
	rbfaValue <<- tclVar("BH-LSU")
	cbpl1Value <<- tclVar("0")
	cbpl2Value <<- tclVar("0")
	cbpl3Value <<- tclVar("1")
	
	cbfa1Value <<- tclVar("1")
	cbfa2Value <<- tclVar("0")
	cbfa3Value <<- tclVar("1")
	cbfa4Value <<- tclVar("1")

	tkconfigure(rbpm1,variable=rbpmValue,value="theoretic")
	tkconfigure(rbpm2,variable=rbpmValue,value="resampling")
	tkconfigure(rbno1,variable=rbnoValue,value="normalized")
	tkconfigure(rbno2,variable=rbnoValue,value="unnormalized")
	tkconfigure(rbfa1,variable=rbfaValue,value="BH-LSU")
	tkconfigure(rbfa2,variable=rbfaValue,value="adaptive")
	tkconfigure(rbfa3,variable=rbfaValue,value="upper.est")
	tkconfigure(rbfa4,variable=rbfaValue,value="point.est")
	tkconfigure(cbpl1,variable=cbpl1Value)
	tkconfigure(cbpl2,variable=cbpl2Value)
	tkconfigure(cbpl3,variable=cbpl3Value)

	tkconfigure(cbfa1,variable=cbfa1Value)
	tkconfigure(cbfa2,variable=cbfa2Value)
	tkconfigure(cbfa3,variable=cbfa3Value)
	tkconfigure(cbfa4,variable=cbfa4Value)

	tkgrid(tklabel(framePmethod,text="P.method:"))
	tkgrid(tklabel(framePmethod,text="No resampling "),rbpm1)
	tkgrid(tklabel(framePmethod,text="Resampling "),rbpm2)



	tkgrid(tklabel(frameAdj,text="FDR adjustment:"),tklabel(frameAdj,text=""),(tklabel(frameAdj,text="Resampling?")))
	tkgrid(tklabel(frameAdj,text="BH-LSU "),rbfa1,cbfa1)
	tkgrid(tklabel(frameAdj,text="Adaptive-two stages"),rbfa2,cbfa2)
	tkgrid(tklabel(frameAdj,text="Upper estimation "),rbfa3,cbfa3)
	tkgrid(tklabel(frameAdj,text="Point estimation "),rbfa4,cbfa4)
	
	tkconfigure(cbfa2, state="disabled")
	tkconfigure(cbfa3, state="disabled")
	tkconfigure(cbfa4, state="disabled")
	
	tkgrid(tklabel(frameAdj,text="Number of permutations:       "),entry.Perms,columnspan=2)
	tkgrid(tklabel(frameplot,text="Plot:"))


	BrowseInput.but <<-tkbutton(frameInput,text="    Browse    ",command=OnBrowseInput)
	BrowseDesign.but <<-tkbutton(frameInput,text="  Select a Vector  ",command=OnBrowseDesign)
	BrowseOutput.but <<-tkbutton(frameOutput,text="    Browse    ",command=OnBrowseOutput)

	Compute.but <<-tkbutton(frameCompute,text="   Go!   ",command=OnCompute)
	Help.but <<-tkbutton(frameCompute,text=" Help ",command=OnHelp)	

	tkwm.title (fdrtt,"FDR")
	tkwm.resizable(fdrtt, 0, 0)

	tkgrid(tklabel(frameInput,text="Input Filename: "),columnspan=2)
	tkgrid(entry.Filename.input,columnspan=2)
	
	tkgrid(BrowseInput.but,columnspan=2)
	
	tkgrid(tklabel(frameInput,text="Design:"),columnspan=2)
	tkgrid(entry.Design,columnspan=2)	
	
	tkgrid(BrowseDesign.but,columnspan=2)
	

	tkgrid(tklabel(frameInput,text=" "))
	tkgrid(tklabel(frameInput,text="Data is:"))
	tkgrid(tklabel(frameInput,text="Normalized "),rbno1)
	tkgrid(tklabel(frameInput,text="Unnormalized "),rbno2)


	tkgrid(tklabel(frameOutput,text="Output .Rda Filename: "))
	tkgrid(entry.Filename.output)
	tkgrid(BrowseOutput.but)
	


	tkgrid(tklabel(frameOutput,text=" "))
	tkgrid(tklabel(frameOutput,text=" "))
	tkgrid(tklabel(frameOutput,text=" "))
	tkgrid(tklabel(frameOutput,text=" "))
	tkgrid(tklabel(frameOutput,text=" "))
	tkgrid(tklabel(frameOutput,text=" "))	
	tkgrid(tklabel(frameOutput,text=" "))
	

	tkbind(fdrtt, "<Return>",OnBrowseInput)
	tkbind(fdrtt, "<Return>",OnBrowseDesign)
	tkbind(fdrtt, "<Return>",OnBrowseOutput)

	tkbind(fdrtt, "<Return>",OnCompute)
	tkbind(fdrtt, "<Return>",OnHelp)

	tkgrid(tklabel(frameplot,text="P-values VS Rank "),cbpl1)
	tkgrid(tklabel(frameplot,text="Adjusted P-values VS Rank "),cbpl2)
	tkgrid(tklabel(frameplot,text="Adjusted P-values VS Observed statistic "),cbpl3)
	
	tkgrid(Help.but,Compute.but)
	
	#tkgrid(frameFileOutput)
	tkgrid(frameInput,frameOutput,sticky="n",ipadx=5,ipady=5)
	tkgrid(frameOutput,ipady=10)
	
	tkgrid(frameAdj,columnspan=2,ipadx=33,ipady=5)
	tkgrid(frameplot,sticky="s",columnspan=2,ipadx=24)
	tkgrid(frameCompute,columnspan=2)
	
	
	tkgrid(frameCompute)
	tkgrid(frameOverall)

	tkfocus(fdrtt)


}
