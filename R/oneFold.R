oneFold <-
function(data, train, test, Kmax, lambda)
{
	## Data dimensions
	m <- nrow(data)
	n <- ncol(data)
	
	## Fit piecewise Gaussian model on training data
	fit <- GGS(data[train,], Kmax, lambda)
	Keff <- length(fit$logLik) - 1
	
	## Calculate standard log-likelihoods (w/o regularization)
	LLtrain <- LLtest <- rep(NA,Kmax+1)
	for (k in 1:(Keff+1)) {
		if (is.infinite(fit$logLik[k]))
			{LLtrain[k] <- LLtest[k] <- -Inf; next}
			
		## Calculate mean and covariance for each segment
		breaks <- fit$breaks[[k]] # breakpoints in training set
		meancovs <- GGSMeanCov(data[train,], breaks, lambda)					
		
		## Calculate log-likelihood for all sample points
		LL <- numeric(m)	
		## Breakpoints in entire set
		breaks <- c(1,train[breaks[-c(1,k+1)]],m+1) 
		for (j in 1:k) {
			mu <- meancovs$mean[,j]
			sig <- meancovs$cov[,,j]
			idx <- breaks[j]:(breaks[j+1]-1)
			R <- chol(sig)
			e <- forwardsolve(t(R),t(data[idx,,drop=FALSE])-mu)
			LL[idx] <- -0.5 * colSums(e^2) - sum(log(diag(R)))
		}
		## Average log-likelihood
		LLtrain[k] <- mean(LL[train])
		LLtest[k] <- mean(LL[test])
	}
	## Put back additive constant
	LLtrain <- LLtrain - 0.5 * n * log(2*pi)
	LLtest <- LLtest - 0.5 * n * log(2*pi)
	return(c(LLtrain,LLtest))
}
