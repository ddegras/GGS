\name{GGSCrossVal}
\alias{GGSCrossVal}

\title{
Greedy Gaussian Segmentation: Cross-Validation
}

\description{
10-fold cross-validation for selecting the number of breakpoints 
and the regularization parameter in the model fit
}

\usage{
GGSCrossVal(data, Kmax = 25, lambdaList = c(0.1, 1, 10), features = NULL, 
	verbose = FALSE, parallel = FALSE)
}

\arguments{
  \item{data}{
data matrix (rows=time points, columns=variables)
}
  \item{Kmax}{
maximum number of (inner) breakpoints allowed
}
  \item{lambdaList}{
regularization parameter values (numeric vector)
}
  \item{features}{
Variables to include in study (integer vector). 
Defaut = \code{NULL} meaning all variables)
}
  \item{verbose}{
logical; if TRUE, print out segmentations and log-likelihood. 
Default = \code{FALSE} 
}
\item{parallel}{logical; if \code{TRUE}, run calculations in parallel. 
	Default = \code{FALSE}  }
}
\details{
To perform the cross-validation in parallel, 
a parallel backend must be registered first. 
See example below. 
}

\value{
A list with components
\item{LLtrain}{matrix of average log-likelihood per sample point for the training data 
	(rows=\# inner break points, columns=regularization parameters)}
\item{LLtest}{matrix of average log-likelihood per sample point for the testing data 
	(rows=\# inner break points, columns=regularization parameters)}
\item{KBest}{number $K$ of (inner) breakpoints yielding the maximum 
	(average) log-likelihood with the testing data} 
\item{lambdaBest}{number $\lambda$ yielding the maximum 
	(average) log-likelihood with the testing data}
\item{Kmax}{input argument \code{Kmax}}
\item{lambdaList}{input argument \code{lambdaList}}
}
\references{
Hallac, Nystrup, and Boyd (2019). Greedy Gaussian Segmentation of Multivariate Time Series. 
}

\seealso{
\code{\link{calculateLikelihood}}, \code{\link{GGS}}, \code{\link{GGSMeanCov}} 

}

\examples{

## For calculations in parallel
\dontrun{
library(foreach)
library(doParallel) # other parallel backends like 'doSNOW' or 'doMC' work too
registerDoParallel()
fpath <- system.file("extdata","Returns.txt",package="GGS")
x <- as.matrix(read.table(fpath))
cv <- GGSCrossVal(x, Kmax=50, lambdaList=c(.1,1,2,5,10,20,50), parallel=TRUE)
}
	
}

