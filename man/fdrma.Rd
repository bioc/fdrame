\name{fdr.ma}
\alias{fdr.ma}
\title{ fdr.ma}
\description{
This function takes normalized expression data array, experimental design and 
computes adjusted p-values. It returns the fdr adjusted p-values and plots, according to the methods
described in (Reiner, Yekutieli and Benjamini 2002)
}
\usage{fdr.ma(exp.arr,design,p.method="resampling",fdr.adj="BH-LSU",plot=c("pvlVSrank","adjVSstat"),perms.num=100)}
\arguments{
\item{exp.arr}{A data frame or matrix, with m rows corresponding to the gene expression level of an observation and n columns to the observations. In the case of gene expression data, rows correspond to genes and columns to mRNA samples. The data can be read using 'read.table'. The matrix is as the X matrix in multtest package.}
\item{design}{A vector of charachters corresponding to observation (column) class labels. For k classes, the labels must be k different charachters. The vector is similar to classlabel vector in multtest package.}
\item{p.method}{A character string specifying whether resampling mehtod may be used. If 'p.method="resampling"', resampling method is used. If 'p.method="t.test"', t.test is computed without resampling.}
\item{fdr.adj}{A character string specifying which type of p-value adjustment method may be used. If 'fdr.adj="BH-LSU"', Benjamini-Hochberg Linear Step Up porcedure is used. If 'fdr.adj="point.est"', point estimation procedure is used.(p.method value must be "resampling") If 'fdr.adj="upper.est"', upper estimation procedure is used.(p.method value must be "resampling") If 'fdr.adj="adaptive"', two-stage adaptive procedure is used.}
\item{plot}{A vector of character strings specifying the plots that may be drawn. 
If the vector contains "pvlVSrank", a plot of p-values vs rank is drawn. 
If the vector contains "adjVSrank", a plot of adjusted p-values vs rank is drawn.
If the vector contains "adjVSstat", a plot of adjusted p-values vs computed statistic (f or t, depends depends on the number of groups on the design) is drawn.  }
\item{perms.num}{An integer specifying the number of permutations that would be used.}
}
\value{A data frame with components
\item{adj}{adjusted p-values.}
\item{dif}{When there are only two groups dif contains the differences between the means of the two groups. when there are more than two groups dif contains the multiple R-Squared values.}
\item{res.pvalues}{resampling p-values. (in case resampling was requested)}
}
\references{Reiner A, Yekutieli D, Benjamini Y: Identifying differentially expressed genes using false discovery rate controlling procedures. Bioinformatics 19:368-375, 2003
Benjamini, Y., Krieger, A.M.,Yekutieli, D. (2001) "Two Staged Linear Step Up FDR Controlling Procedure", Technical Report Department of Statistics and O.R., Tel Aviv University.
}
\author{Yoav Benjamini, Ephraim Kenigsberg, Anat Reiner, Daniel Yekutieli}
\examples{
data(exp.arr,design)
output<-fdr.ma(exp.arr,design,p.method="theoretic",fdr.adj="BH-LSU",plot=c("pvlVSrank","adjVSrank","adjVSstat"))
}
\keyword{htest}
