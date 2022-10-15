GGSCrossVal <-
function(data, Kmax = 25, lambdaList = c(0.1, 1, 10), 
	 features = NULL, parallel = FALSE, verbose = FALSE)
{
	
	## Extract features if required
	if (!is.null(features))
		data <- data[,features]
		
	## Data dimensions	
	m <- nrow(data)
	n <- ncol(data)
	
	## Set up (10) folds
	set.seed(0)
	ordering <- sample(m)
	foldBreaks <- round(seq.int(0,m,length=11))
	
	## Check parallel computing status
	test.parallel <- parallel && require(foreach) && getDoParRegistered()
	if (parallel && !test.parallel)
		warning(paste("Library 'foreach' not installed and/or",
		"no parallel backend registered.", 
		"Running calculations sequentially."))
	if (test.parallel && verbose) verbose <- FALSE
	
	nlam <- length(lambdaList)
	LL <- matrix(,2*(Kmax+1),nlam)

	## Run 10-fold cross-validation for each lambda
	nfold <- 10
	for (l in 1:nlam) {
		lambda <- lambdaList[l]
		if (test.parallel) { # PARALLEL EXECUTION
			result <- foreach(i=1:nfold, .combine=cbind) %dopar% {
				test <- ordering[(foldBreaks[i]+1):foldBreaks[i+1]] # testing samples
				train <- (1:m)[-test] # training samples			
				oneFold(data, train, test, Kmax, lambda)		
			} 			
		} else { # SEQUENTIAL EXECUTION
			result <- matrix(,2*(Kmax+1),nfold)
			for (i in 1:nfold) {
				test <- ordering[(foldBreaks[i]+1):foldBreaks[i+1]]
				train <- (1:m)[-test]
				result[,i] <- oneFold(data, train, test, Kmax, lambda)			
			}  
		} 
		LL[,l] <- rowMeans(result)						
	}
	
	## Training and testing average log-likelihoods
	LLtrain <- LL[1:(Kmax+1),]
	LLtest <- LL[(Kmax+2):(2*Kmax+2),]
	
	## Best K and lambda (for testing)
	# browser()
	idx <- arrayInd(which.max(LLtest),c(Kmax+1,nlam))
	KBest <- idx[1] - 1
	lambdaBest <- lambdaList[idx[2]]
	
	return(list(LLtrain=LLtrain, LLtest=LLtest, 
		KBest=KBest, lambdaBest=lambdaBest, Kmax=Kmax, 
		lambdaList=lambdaList))
	
}
