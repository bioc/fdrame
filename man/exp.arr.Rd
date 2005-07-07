\name{exp.arr}
\alias{exp.arr}
\docType{data}
\title{ Normalized Expression Array }
\description{The source of this Data is a lipid metablolism study described and  analyzed in Dudoit,S., Yang,Y.H, Callow,M.J. and Speed,T.P. (2002) Statistical Methods for Identifying Differentially Expressed Genes in Replicated cDNA Microarray Experiments. Stat Sinica, 12, 111-139.
}
\usage{data(exp.arr)}
\format{
  The format is:
 num [1:6384, 1:16] -0.2500  0.0329 -0.2065 -0.2240 -0.8542 ...
 - attr(*, "dimnames")=List of 2
  Each row represents a gene: chr [1:6384] "1" "2" "3" "4" ...
  Each column represents a sample : chr [1:16] "X1" "X2" "X3" "X4" ...
}
\details{
	The data is normalized.
}
\source{
 http://www.stat.berkeley.edu/users/terry/zarray/Html/matt.html	
}
\references{
  Reiner A, Yekutieli D, Benjamini Y: Identifying differentially expressed genes using false discovery rate controlling procedures. Bioinformatics 19:368-375, 2003
}
\examples{
data(exp.arr)
}
\keyword{datasets}
