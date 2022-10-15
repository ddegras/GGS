GGS <-
function(data, Kmax, lambda, features = NULL, verbose = FALSE)
{

	if (!is.null(features))
		data <- data[,features]
    m <- nrow(data) # number of time points
    n <- ncol(data) # number of variables
    
    ## Initialize breakpoints
	breaks = c(1,m+1)
	breakPoints <- vector("list",Kmax+1)
    breakPoints[[1]] <- breaks
    plotPoints <- numeric(Kmax+1)
    ll <- calculateLikelihood(data,breaks,lambda)
    plotPoints[1] <- ll
	if (Kmax == 0) 
		return(list(breaks=breakPoints, logLik=plotPoints))
		
    ## Start GGS Algorithm
    for (K in 1:Kmax) {  	
        ## For each segment, find breakpoint and increase in LL
        newVal <- 0
        for (i in 1:K) {
        		idx <- breaks[i]:(breaks[i+1]-1)
            out <- addBreak(data[idx,],lambda)
            ind <- out[1]
            val <- out[2] 
            if (val > newVal) {
                newInd <- ind + breaks[i]
                newVal <- val
			}
		}
			
        ## Check algorithm termination
        if (newVal == 0) {
        		if (verbose) 
        			cat("We are done adding breakpoints!\n",
	            		"Breakpoints:",breaks,"\n")
			breakPoints <- breakPoints[1:K]
			plotPoints <- plotPoints[1:K]
            break
		}
		
        ## Add new breakpoint
		breaks <- sort(c(breaks,newInd))
		if (verbose) {
			cat("Breakpoint occurs at sample number:", 
				newInd,", LL =",ll+newVal,"\n")
			cat("Breakpoints:",breaks,"\n")
		}
		
		## Adjust current locations of the breakpoints
		# browser()
		breaks <- adjustBreaks(data,breaks,lambda,verbose)

		## Calculate likelihood
		ll <- calculateLikelihood(data,breaks,lambda)
		breakPoints[[K+1]] <- breaks
		plotPoints[K+1] <- ll
	}

	return(list(breaks=breakPoints, logLik=plotPoints))
}
