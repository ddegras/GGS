# GGS
Greedy Gaussian Segmentation

`R` Implementation of the Greedy Gaussian Segmentation algorithm of Hallac, Nystrup, and Boyd (2019) for multivariate time series. 
Original Python code: https://github.com/cvxgrp/GGS

**To install the `R` package:**
```
library(devtools)
install_github("https://github.com/ddegras/GGS")
```

**Package functions:**
- `GGS`: Finds breakpoints in the data for a given regularization parameter.
- `GGSMeanCov`: Finds the means and regularized covariances of each segment, given a set of breakpoints.
- `GGSCrossVal`: Runs 10-fold cross validation, and returns the train and test set likelihood for all specified numbers of breakpoints and regularization parameters.

#### Reference:
Hallac, D., Nystrup, P. & Boyd, S. Greedy Gaussian segmentation of multivariate time series. *Advances in Data Analysis and Classification* 13, 727–751 (2019). https://doi.org/10.1007/s11634-018-0335-0
