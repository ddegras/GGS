addBreak <-
function(data, lambda)
{

	## Data dimensions
	m <- nrow(data) # number of time points
	n <- ncol(data) # number of variables

	## Trivial case: 1 time point
	if (is.null(m) || m == 1) return(c(0,0))
	
	## Log-likelihook for entire segment
	origLL <- calculateLikelihood(data,c(1,m+1),lambda)
	
	## Initialize best performance
	maxLL <- origLL
	maxInd <- 0
	
	## Data sums and cross-products
	totSum <- colSums(data)
	runSum <- numeric(n)	
	totCP <- crossprod(data)
	runCP <- matrix(0,n,n)

	## Loop over candidate breakpoints
	lambdapositive <- (lambda > 0)
	for (i in 1:(m-1)) {
		runSum <- runSum + data[i,]
		runCP <- runCP + tcrossprod(data[i,])
		## Contribution of left segment
		mu <- runSum/i
		sig <- runCP/i - tcrossprod(mu)
		if (lambdapositive)
			diag(sig) <- diag(sig) + (2*lambda/i)
		R <- tryCatch(chol(sig), error = function(e) NULL)
		if (is.null(R) || kappa(sig) > 1e8) next
		LL <- -i * sum(log(diag(R))) 
		## Contribution of right segment
		mu <- (totSum-runSum)/(m-i)
		sig <- (totCP-runCP)/(m-i) - tcrossprod(mu)
		if (lambdapositive)
			diag(sig) <- diag(sig) + (2*lambda/(m-i))
		R <- tryCatch(chol(sig), error = function(e) NULL)
		if (is.null(R) || kappa(sig) > 1e8) next
		LL <- LL - (m-i) * sum(log(diag(R)))
		if (LL > maxLL) {
			maxLL <- LL
			maxInd <- i
		}				
	}
	
	return(c(maxInd,maxLL-origLL))
}
