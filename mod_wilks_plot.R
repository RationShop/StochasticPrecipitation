# plot summary stats

if(plotSummary) {
  
  pdf(plotFileName)
  
  # compare number of rainy days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rainyDaysObs <- apply(occDataObs, 2, FUN=FnCountMonthlySum, daysInMonth)
  rainyDaysSyn <- apply(occDataSyn, 1, FUN=FnCountMonthlySum, daysInMonth)
  boxplot(rainyDaysObs, 
          range=0, 
          at = 1:numStn*3 + 1, 
          xlim=c(3, 3 + numStn*3), 
          ylim = c(0,30), 
          xaxt = "n", 
          main="Monthly Rainy Days")
  boxplot(rainyDaysSyn, 
          range=0, 
          at = 1:numStn*3 + 2, 
          xaxt = "n", 
          border="blue", 
          add = TRUE)
  axis(1, at = 1:numStn*3 + 1.5, labels = 1:numStn, tick = TRUE)
  
  # compare monthly total ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  monTotObs <- apply(amtDataObs, 2, FUN=FnCountMonthlySum, daysInMonth)
  monTotSyn <- apply(amtDataSyn, 1, FUN=FnCountMonthlySum, daysInMonth)
  boxplot(monTotObs, 
          range=0, 
          at = 1:numStn*3 + 1, 
          xlim=c(3, 3 + numStn*3), 
          ylim = c(0,50), 
          xaxt = "n", 
          main="Monthly Total (inches)")
  boxplot(monTotSyn, 
          range=0, 
          at = 1:numStn*3 + 2, 
          xaxt = "n", 
          border="blue", 
          add = TRUE)
  axis(1, at = 1:numStn*3 + 1.5, labels = 1:numStn, tick = TRUE)
  
  # compare means, number of wet days and monthly totals ~~~~~~~~~~~~~~~~~~~~~~~
  meanWetDaysObs <- apply(rainyDaysObs, 2, FUN=mean)
  meanWetDaysSyn <- apply(rainyDaysSyn, 2, FUN=mean)
  meanMonTotObs  <- apply(monTotObs, 2, FUN=mean)
  meanMonTotSyn  <- apply(monTotSyn, 2, FUN=mean)
  
  # means, daily all and daily wet days
  meanDlyObs <- apply(amtDataObs, 2, FUN=mean, na.rm=TRUE)
  meanDlySyn <- apply(amtDataSyn, 1, FUN=mean, na.rm=TRUE)
  wetDataSyn <- amtDataSyn
  wetDataSyn[wetDataSyn <= prcpTHRESH] <- NA #only for the wet days
  meanDlyWetSyn <- rowMeans(wetDataSyn, na.rm=TRUE)
  
  # sd and skew
  sdDlyObs <- apply(amtDataObs, 2, FUN=sd, na.rm=TRUE)
  sdDlySyn <- apply(amtDataSyn, 1, FUN=sd, na.rm=TRUE)
  skewDlyObs <- apply(amtDataObs, 2, FUN=FnSkew)
  skewDlySyn <- apply(amtDataSyn, 1, FUN=FnSkew)
  
  par(mfrow=c(3,2))
  # panel 1
  plot(meanWetDaysObs, 
       meanWetDaysSyn, 
       xlab="OBS", 
       ylab="SIM", 
       main="Mean No of Wet Days")
  abline(0,1)
  # panel 2
  plot(meanMonTotObs, 
       meanMonTotSyn, 
       xlab="OBS", 
       ylab="SIM", 
       main="Mean Monthly Total")
  abline(0,1)
  # panel 3
  plot(meanDlyObs, 
       meanDlySyn, 
       xlab="OBS", 
       ylab="SIM", 
       main="Mean Daily")
  abline(0,1)
  # panel 4
  plot(meanDlyWetObs, 
       meanDlyWetSyn, 
       xlab="OBS", 
       ylab="SIM", 
       main="Mean Daily, Wet Days Only")
  abline(0,1)
  # panel 5
  plot(sdDlyObs, 
       sdDlySyn, 
       xlab="OBS", 
       ylab="SIM", 
       main="Std. Dev. Daily")
  abline(0,1)
  # panel 1
  plot(skewDlyObs, 
       skewDlySyn, 
       xlab="OBS", 
       ylab="SIM", 
       main="Skew Daily")
  abline(0,1)
  par(mfrow=c(1,1))
  
  # compare markov probs and joint probs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # markov probs
  occProbsMrkSyn <- NULL
  for (eachStn in 1:numStn) {
    occProbsMrkSyn <- rbind(occProbsMrkSyn, FnComputeMarkovProbs(occDataSyn[eachStn,]))
  }
  occProbsMrkSyn <- as.data.frame(occProbsMrkSyn)
  colnames(occProbsMrkSyn) <- c("n01", "n00", "n11", "n10")
  occProbsMrkSyn$p01 <- occProbsMrkSyn$n01/(occProbsMrkSyn$n01 + occProbsMrkSyn$n00)
  occProbsMrkSyn$p11 <- occProbsMrkSyn$n11/(occProbsMrkSyn$n11 + occProbsMrkSyn$n10)
  
  # joint probs
  occProbsJntSyn00 <- array(0, dim=c(numStn,numStn))
  occProbsJntSyn01 <- array(0, dim=c(numStn,numStn))
  occProbsJntSyn11 <- array(0, dim=c(numStn,numStn))
  occProbsJntSyn10 <- array(0, dim=c(numStn,numStn))
  for (stn1 in 1:(numStn-1)) {
    for (stn2 in stn1:numStn) {
      
      jointProbs <- FnComputeJointProbs(occDataSyn[stn1,], occDataSyn[stn2,])
      
      occProbsJntSyn00[stn1,stn2] <- jointProbs$pi00
      occProbsJntSyn01[stn1,stn2] <- jointProbs$pi01
      occProbsJntSyn11[stn1,stn2] <- jointProbs$pi11
      occProbsJntSyn10[stn1,stn2] <- jointProbs$pi10    
    }
  }
  occProbsJntSyn00 <- occProbsJntSyn00 + t(occProbsJntSyn00)
  diag(occProbsJntSyn00) <- 1
  occProbsJntSyn01 <- occProbsJntSyn01 + t(occProbsJntSyn01)
  diag(occProbsJntSyn01) <- 1
  occProbsJntSyn11 <- occProbsJntSyn11 + t(occProbsJntSyn11)
  diag(occProbsJntSyn11) <- 1
  occProbsJntSyn10 <- occProbsJntSyn10 + t(occProbsJntSyn10)
  diag(occProbsJntSyn10) <- 1
  
  par(mfrow=c(2,2))
  plot(occProbsMrkObs$p01, 
       occProbsMrkSyn$p01, 
       xlab="OBS", 
       ylab="SIM", 
       main="p01")
  abline(0,1)
  plot(occProbsMrkObs$p11, 
       occProbsMrkSyn$p11, 
       xlab="OBS", 
       ylab="SIM", 
       main="p11")
  abline(0,1)
  
  plot(occProbsJntObs11[lower.tri(occProbsJntObs11)], 
       occProbsJntSyn11[lower.tri(occProbsJntSyn11)], 
       xlab="OBS", 
       ylab="SIM", 
       main="Joint pi_11")
  abline(0,1)
  plot(occProbsJntObs00[lower.tri(occProbsJntObs00)], 
       occProbsJntSyn00[lower.tri(occProbsJntSyn00)], 
       xlab="OBS", 
       ylab="SIM", 
       main="Joint pi_00")
  abline(0,1)
  par(mfrow=c(1,1))
  
  # cross correlation, occurrence and amount ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  occCrossCorObs <- cor(occDataObs, use="pairwise.complete.obs")
  occCrossCorSyn <- cor(t(occDataSyn), use="pairwise.complete.obs")
  
  amtCrossCorObs <- cor(wetDataObs, use="pairwise.complete.obs", method="spearman")
  amtCrossCorSyn <- cor(t(wetDataSyn), use="pairwise.complete.obs", method="spearman")
  
  amtCrossCorObs2 <- cor(wetDataObs, use="pairwise.complete.obs")
  amtCrossCorSyn2 <- cor(t(wetDataSyn), use="pairwise.complete.obs")
  
  par(mfrow=c(2,2))
  plot(occCrossCorObs[lower.tri(occCrossCorObs)], 
       occCrossCorSyn[lower.tri(occCrossCorSyn)], 
       xlab="OBS", 
       ylab="SIM", 
       main="Occurrence Cross-Correlation")
  abline(0,1)
  plot(amtCrossCorObs[lower.tri(amtCrossCorObs)], 
       amtCrossCorSyn[lower.tri(amtCrossCorSyn)], 
       xlab="OBS", 
       ylab="SIM", 
       main="Amount Cross-Correlation (Rank)")
  abline(0,1)
  plot(amtCrossCorObs2[lower.tri(amtCrossCorObs2)], 
       amtCrossCorSyn2[lower.tri(amtCrossCorSyn2)], 
       xlab="OBS", 
       ylab="SIM", 
       main="Amount Cross-Correlation (Product-Moment)")
  abline(0,1)
  par(mfrow=c(1,1))
  
  # continuity ratio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # as defined by Wilks (1998) and implmented by Brissette et al (2007)
  # joint probs
  occContRatioObs <- array(0, dim=c(numStn,numStn))
  occContRatioSyn <- array(0, dim=c(numStn,numStn))
  amtContRatioObs <- array(0, dim=c(numStn,numStn))
  amtContRatioSyn <- array(0, dim=c(numStn,numStn))
  
  for (stn1 in 1:numStn) {
    for (stn2 in 1:numStn) {    
      occContRatioObs[stn1,stn2] <- FnComputeOccContinuityRatio(occDataObs[,stn1], occDataObs[,stn2])
      occContRatioSyn[stn1,stn2] <- FnComputeOccContinuityRatio(occDataSyn[stn1,], occDataSyn[stn2,])
      amtContRatioObs[stn1,stn2] <- FnComputeAmtContinuityRatio(amtDataObs[,stn1], amtDataObs[,stn2])
      amtContRatioSyn[stn1,stn2] <- FnComputeAmtContinuityRatio(amtDataSyn[stn1,], amtDataSyn[stn2,])    
      
    }
  }
  diag(occContRatioObs) <- NA
  diag(occContRatioSyn) <- NA
  diag(amtContRatioObs) <- NA
  diag(amtContRatioSyn) <- NA
  
  par(mfrow=c(1,2))
  plot(as.vector(occContRatioObs), 
       as.vector(occContRatioSyn), 
       xlab="OBS", 
       ylab="SIM", 
       main="Occurrence Continuity Ratio")
  abline(0,1)
  plot(as.vector(amtContRatioObs), 
       as.vector(amtContRatioSyn), 
       xlab="OBS", 
       ylab="SIM", 
       main="Amount Continuity Ratio")
  abline(0,1)
  par(mfrow=c(1,1))
  
  garbage <- dev.off()
}
