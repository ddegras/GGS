adjustBreaks <-
function(data, breaks, lambda, 
	verbose = FALSE, maxShuffles = 250)
{
	K <- length(breaks) - 2 # number of inner breakpoints
	if (K <= 1) return(breaks)
	## Initialize "active" breakpoints (i.e. that are changed in current pass) 
	active <- c(FALSE,rep(TRUE,K),FALSE) 
	set.seed(0)
	for (z in 1:maxShuffles) {
		ordering <- sample(K)
		## Breakpoints active at previous pass
		active.prev <- active
		for (i in ordering) {
			## If breakpoints b(i) and b(i+2) 
			## were not changed in previous pass, skip
			if (!(active.prev[i] || active.prev[i+2])) 
				{ active[i+1] <- FALSE; next }
			idx <- breaks[i]:(breaks[i+2]-1)
			out <- addBreak(data[idx,],lambda)
			ind <- out[1]
			val <- out[2]
			t <- breaks[i] + ind # adjusted breakpoint
			if (t != breaks[i+1]) {
				active[i+1] <- TRUE
				bp <- breaks[i+1]
				breaks[i+1] <- t				
				if (verbose) cat("Breakpoint at sample number",
					bp,"moved to",t,"LL = ",
					calculateLikelihood(data,breaks,lambda),"\n")
			} else {
				active[i+1] <- FALSE
			}
		}
		if (all(!active)) break
	}
	return(breaks)
}
