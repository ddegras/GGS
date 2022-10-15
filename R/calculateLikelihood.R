calculateLikelihood <-
function(data, breaks, lambda)
{
	K <- length(breaks) - 2
	n <- ncol(data)
	LL <- 0
	for (i in 1:(K+1)) {
		m <- breaks[i+1] - breaks[i] 
		sig <- if (m > 1) {
			cov(data[breaks[i]:(breaks[i+1]-1),]) * ((m-1)/m)
		} else { matrix(0,n,n) }
		if (lambda > 0) 	
			diag(sig) <- diag(sig) + (2*lambda/m)
		## Attempt Cholesky decomposition
		R <- tryCatch(chol(sig), error = function(e) NULL)
		if (is.null(R)) return(-Inf)
		LL <- LL - m * sum(log(diag(R)))
	}
	# Full log-likelihood
	# LL <- LL - 0.5 * length(data) - 0.5 * length(data) * log(2*pi)
	return(LL)
}
