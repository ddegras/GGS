\name{calculateLikelihood}
\alias{calculateLikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Greedy Gaussian Segmentation: Calculate Lok-Likelihood for Model fit.
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateLikelihood(data, breaks, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, breaks, lambda) 
{
    K <- length(breaks) - 2
    n <- ncol(data)
    LL <- 0
    for (i in 1:(K + 1)) {
        m <- breaks[i + 1] - breaks[i]
        sig <- if (m > 1) {
            cov(data[breaks[i]:(breaks[i + 1] - 1), ]) * ((m - 
                1)/m)
        }
        else {
            matrix(0, n, n)
        }
        if (lambda > 0) 
            diag(sig) <- diag(sig) + (2 * lambda/m)
        R <- tryCatch(chol(sig), error = function(e) NULL)
        if (is.null(R)) 
            return(-Inf)
        LL <- LL - m * sum(log(diag(R)))
    }
    return(LL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
