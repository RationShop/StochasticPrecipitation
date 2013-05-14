# library of functions for the modified wilks method

#-------------------------------------------------------------------------------
# function to compute markov transition probabilities for each station, Eq 1 & 2
FnComputeMarkovProbs <- function(tsData) {
  
  occ1 <- tsData[1:(length(tsData)-1)]
  occ2 <- tsData[2:length(tsData)]
  
  n01 <- ifelse(occ1==0 & occ2==1, 1, 0)
  n00 <- ifelse(occ1==0 & occ2==0, 1, 0)
  n11 <- ifelse(occ1==1 & occ2==1, 1, 0)
  n10 <- ifelse(occ1==1 & occ2==0, 1, 0)
  
  return (c(sum(n01,na.rm=TRUE), sum(n00,na.rm=TRUE), sum(n11,na.rm=TRUE), sum(n10,na.rm=TRUE)))
}

#-------------------------------------------------------------------------------
# function to compute joint station probabilities
FnComputeJointProbs <- function(ts1, ts2) {
  
  ind1 <- which(!is.na(ts1))
  ind2 <- which(!is.na(ts2))
  
  inds <- intersect(ind1, ind2)
  ts1 <- ts1[inds]
  ts2 <- ts2[inds]
  
  x00 <- sum(ifelse(ts1==0 & ts2==0, 1, 0), na.rm=TRUE)
  x01 <- sum(ifelse(ts1==0 & ts2==1, 1, 0), na.rm=TRUE)
  x11 <- sum(ifelse(ts1==1 & ts2==1, 1, 0), na.rm=TRUE)
  x10 <- sum(ifelse(ts1==1 & ts2==0, 1, 0), na.rm=TRUE)
  
  pi00 <- x00/length(ts1)
  pi01 <- x01/length(ts1)
  pi11 <- x11/length(ts1)
  pi10 <- x10/length(ts1)
  
  return (list("pi00"=pi00, "pi01"=pi01, "pi11"=pi11, "pi10"=pi10))
}

#-------------------------------------------------------------------------------
# function to compute continuity ratio for amounts
FnComputeAmtContinuityRatio <- function(ts1, ts2) {
  
  ind1 <- which(!is.na(ts1))
  ind2 <- which(!is.na(ts2))
  
  inds <- intersect(ind1, ind2)
  ts1 <- ts1[inds]
  ts2 <- ts2[inds]
  
  #numerator, E(Xi|Xi>0 & Xj==0)
  numer <- sum(ifelse(ts1>prcpTHRESH & ts2<=prcpTHRESH, ts1, NA), na.rm=TRUE)
  #denominator, E(Xi|Xi>0 & Xj>0)
  denom <- sum(ifelse(ts1>prcpTHRESH & ts2>prcpTHRESH, ts1, NA), na.rm=TRUE)
  
  return (numer/denom)
}

#-------------------------------------------------------------------------------
# function to compute continuity ratio for occurrence
FnComputeOccContinuityRatio <- function(ts1, ts2) {
  
  ind1 <- which(!is.na(ts1))
  ind2 <- which(!is.na(ts2))
  
  inds <- intersect(ind1, ind2)
  ts1 <- ts1[inds]
  ts2 <- ts2[inds]
  
  #numerator, E(Xi|Xi>0 & Xj==0)
  numer <- sum(ifelse(ts1==1 & ts2==0, ts1, NA), na.rm=TRUE)
  #denominator, E(Xi|Xi>0 & Xj>0)
  denom <- sum(ifelse(ts1==1 & ts2==1, ts1, NA), na.rm=TRUE)
  
  return (numer/denom)
}

#-------------------------------------------------------------------------------
# function to compute the "A" parameter in Eqn 16
FnParamA <- function(x) {
  return (log(mean(x, na.rm=TRUE)) - log(exp(mean(log(x), na.rm=TRUE))))
}

#-------------------------------------------------------------------------------
# function for the 1st order markov Chain precip occurence process; Eqn
FnNewMarkovState <- function(randNum, prob) {
  return (ifelse(randNum <= prob, 1, 0))
}

#-------------------------------------------------------------------------------
# computes monthly totals (number of rainy days, monthly rainfall total, etc)
FnCountMonthlySum <- function(tsData, daysInMonth) {
  newData <- matrix(tsData, nrow=daysInMonth)
  #   sums <- colSums(newData, na.rm=TRUE)
  #   return (mean(sums, na.rm=TRUE))
  return (colSums(newData, na.rm=TRUE))
  
}

# ------------------------------------------------------------------------------
# function to make a correlation matrix positive definite
eigenSmall <- 1e-6
FnFixCorrMatrix <- function(corrMatrix) {
  
  cholStatus <- try(outChol <- chol(corrMatrix), silent=TRUE) 
  cholError  <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  
  # fix the correlation matrix if chol fails
  while(cholError) {
    
    #replace -ve eigen values with a small +ve number
    neweig  <- eigen(corrMatrix)
    neweig$values[neweig$values < 0] <- eigenSmall
    
    #inv = transp for eig vectors
    corrMatrix <- neweig$vectors %*% diag(neweig$values) %*% t(neweig$vectors) 
    
    #try chol again
    cholStatus <- try(outChol <- chol(corrMatrix), silent=TRUE) 
    cholError  <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  }
  
  return (outChol)
}

#-------------------------------------------------------------------------------
# calculate skew
FnSkew <- function(x) {
  x <- x[!is.na(x)] #remove missing values
  return (sum((x - mean(x))^3)/(length(x) * sd(x)^3))     
}
