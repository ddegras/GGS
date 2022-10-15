GGSMeanCov <-
function(data, breakpoints, lambda, features=NULL)
{
    ## Select the desired features
    if (!is.null(features))
        data <- data[,features]
	
	K <- length(breakpoints) - 2 # number of inner breakpoints
    n <- ncol(data)
    mu <- matrix(,n,K+1)
    sig <- array(dim=c(n,n,K+1))

    for (i in 1:(K+1)) {
        ## Get mean and regularized covariance of current segment
        idx <- breakpoints[i]:(breakpoints[i+1]-1)
        m <- breakpoints[i+1] - breakpoints[i]
        if (m > 1) {
	        mu[,i] <- colMeans(data[idx,])
	        sigi <- cov(data[idx,]) * ((m-1)/m)
		} else {
			mu[,i] <- data[idx,]
			sigi <- matrix(0,n,n)
		}
        diag(sigi) <- diag(sigi) + (2*lambda/m) 
        sig[,,i] <- sigi
	}
    return(list(mean=mu, cov=sig))
}
