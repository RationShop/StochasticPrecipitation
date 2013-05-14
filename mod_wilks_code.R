# modified Wilks approach, based on Mhanna & Bauwens (2012)

# parameters
numStn        <- 15      #number of stations
prcpTHRESH    <- 0.01    #precipitation threshold, in inches
daysInMonth   <- 31      #days in month
numYears      <- 1000    #number of years of simulation
numDays       <- daysInMonth * numYears
plotSummary   <- FALSE   # TRUE will plot the summary results, takes a while
plotFileName  <- "results_summary.pdf"

# library file
source("mod_wilks_lib.R")

# input station precipitation data; "janData.txt" has values in inches
# should be consistent with units of prcpTHRESH
prcpData <- read.csv("janData.txt", header=FALSE, as.is=TRUE)

# subsets for analysis
occDataObs <- prcpData[, 3:(2 + numStn)]
occDataObs <- ifelse(occDataObs > prcpTHRESH, 1, 0)
amtDataObs <- prcpData[, 3:(2 + numStn)]


#-------------------------------------------------------------------------------
# estimate occurrence parameters

#markov transition probabilities for each station, Eq 1 & 2
occProbsMrkObs <- NULL
for (eachStn in 1:numStn) {
  occProbsMrkObs <- rbind(occProbsMrkObs, FnComputeMarkovProbs(occDataObs[,eachStn]))
}
occProbsMrkObs <- as.data.frame(occProbsMrkObs)
colnames(occProbsMrkObs) <- c("n01", "n00", "n11", "n10")

occProbsMrkObs$p01 <- occProbsMrkObs$n01/(occProbsMrkObs$n01 + occProbsMrkObs$n00)
occProbsMrkObs$p11 <- occProbsMrkObs$n11/(occProbsMrkObs$n11 + occProbsMrkObs$n10)


#joint probability of occurence
occProbsJntObs00 <- array(0, dim=c(numStn,numStn))
occProbsJntObs01 <- array(0, dim=c(numStn,numStn))
occProbsJntObs11 <- array(0, dim=c(numStn,numStn))
occProbsJntObs10 <- array(0, dim=c(numStn,numStn))
for (stn1 in 1:numStn) {
  for (stn2 in 1:numStn) {
    
    jointProbs <- FnComputeJointProbs(occDataObs[,stn1], occDataObs[,stn2])
    
    occProbsJntObs00[stn1,stn2] <- jointProbs$pi00
    occProbsJntObs01[stn1,stn2] <- jointProbs$pi01
    occProbsJntObs11[stn1,stn2] <- jointProbs$pi11
    occProbsJntObs10[stn1,stn2] <- jointProbs$pi10    
  }
}


# Eqn (21)
# oddsRatio <- (occProbsJntObs00 + occProbsJntObs11)/(occProbsJntObs10 + occProbsJntObs01)
oddsRatio <- (occProbsJntObs00 * occProbsJntObs11)/(occProbsJntObs10 * occProbsJntObs01)
# infinite values along the diagonal of Jnt10 and Jnt01 need to be replaced 
# for cholesky factorization
oddsRatio <- ifelse(is.infinite(oddsRatio), 1, oddsRatio)

# Eqn (20)
gammaCoeff <- (oddsRatio - 1)/(oddsRatio + 1)
diag(gammaCoeff) <- 1

#cholesky stuff
occChol <- FnFixCorrMatrix(gammaCoeff)

#-------------------------------------------------------------------------------
# estimate amount parameters
wetDataObs <- amtDataObs
wetDataObs[wetDataObs <= prcpTHRESH] <- NA #only for the wet days

meanDlyWetObs <- colMeans(wetDataObs, na.rm=TRUE)

aParam     <- apply(wetDataObs, 2, FUN=FnParamA)
alphaParam <- (1 + sqrt(1 + (4*aParam/3)))/(4*aParam)
betaParam  <- meanDlyWetObs/alphaParam

# correl
# amtCor <- cor(amtDataObs, use="pairwise.complete.obs", method="spearman")
amtCor   <- cor(wetDataObs, use="pairwise.complete.obs")
binorCor <- 2*sin(pi*amtCor/6)

# cholesky stuff
amtChol <- FnFixCorrMatrix(binorCor)

#-------------------------------------------------------------------------------
# generate synthetic data

# generate occurrences
occDataSyn <- array(0,dim=c(numStn, numDays))
#for computational ease (markov chain), p01 and p11 are cast as matrices and not vectors
p01 <- matrix(rep(occProbsMrkObs$p01, numDays), ncol=numDays)
p11 <- matrix(rep(occProbsMrkObs$p11, numDays), ncol=numDays)

# random numbers for occurrences
randNormOcc <- array(rnorm(numStn * numDays), dim=c(numStn, numDays)) #dim c(numStn, numDays)
inputrandOcc <- t(occChol) %*% randNormOcc #dim c(numStn, numDays)
inputrandOcc <- pnorm(inputrandOcc)

# synthetic occurrence
state01 <- FnNewMarkovState(inputrandOcc, p01)
state11 <- FnNewMarkovState(inputrandOcc, p11)
for (eachDay in 2:numDays) {
  occDataSyn[,eachDay] <- ifelse(occDataSyn[,eachDay-1] == 0, 
                                 state01[,eachDay], 
                                 state11[,eachDay])
}

# generate amounts
amtDataSyn <- array(0,dim=c(numStn, numDays))

# random numbers for amounts
randNormAmt <- array(rnorm(numStn * numDays), dim=c(numStn, numDays)) #dim c(numStn, numDays)
inputrandAmt <- t(amtChol) %*% randNormAmt #dim c(numStn, numDays)
inputrandAmt <- pnorm(inputrandAmt)

# synthetic amounts
for (eachStn in 1:numStn) {
  amtDataSyn[eachStn,] <- qgamma(inputrandAmt[eachStn,], 
                                 shape=alphaParam[eachStn], 
                                 scale=betaParam[eachStn]) * occDataSyn[eachStn,]
}
amtDataSyn <- ifelse(amtDataSyn > 0 & amtDataSyn < prcpTHRESH, prcpTHRESH, amtDataSyn)

