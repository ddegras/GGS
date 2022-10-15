\name{GGS}
\alias{GGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Greedy Gaussian Segmentation of Multivariate Time Series}
\description{
Fit piecewise Gaussian model to multivariate time series by the GGS approach of Hallac et al. (2019)
}
\usage{
GGS(data, Kmax, lambda, features = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data matrix: rows = time points, columns = variables
}
  \item{Kmax}{
maximum number of (inner) breakpoints (nonnegative integer)
}
  \item{lambda}{
regularization parameter for maximum likelihood estimation
}
  \item{features}{
  of variables to include in study (integer vector; default = \code{NULL} meaning all variables)}
  \item{verbose}{
Logical: if \code{TRUE}, prints out segmentation results
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with components 
\item{breaks}{list of breakpoints by increasing order of \code{K}}
\item{logLik}{associated log-likelihood values}
%% ...
}
\references{
Hallac, Nystrup, and Boyd (2019). Greedy Gaussian Segmentation of Multivariate Time Series. 
}



\seealso{
\code{\link{calculateLikelihood}}, \code{\link{GGSCrossval}}, \code{\link{GGSMeanCov}} 
}

% \examples{
% }

