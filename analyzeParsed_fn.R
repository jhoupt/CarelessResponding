## Model Based Analysis ##
fitWeibullStan <- function(dat) { 
   require(rstan)
   dat <- subset(dat, !is.na(Response))
   dat$Subject <- factor(dat$Subject)
   
   nSubjects <- length(levels(dat$Subject))
   nConditions <- length(levels(dat$Condition))
   nTotalTrials <- length(dat$Response)
   
   subject <- as.numeric(dat$Subject)
   trialOrder <- dat$Order #* 500
   correct <- 1-dat$Response # Convert to 0: Correct; 1: Incorrect
   #correct <- dat$Response 
   chanceCorrect <- 2/length(unique(inData$Response[!is.na(inData$Response)]));
   
   
   conditionX <- with(dat, as.numeric(Condition))
   condition <- rep(NA, nSubjects)
   for ( sj in 1:nSubjects) {
      condition[sj] <- with(dat, conditionX[Subject==levels(Subject)[sj]][1])
   }
   
   
   datalist <- list(nSubjects=nSubjects, nConditions=nConditions, nTotalTrials=nTotalTrials, subject=subject, condition=condition, correct=correct, trialOrder=trialOrder, chanceCorrect=chanceCorrect)
   init1<-list(midpoint_group=rep(1,nConditions), 
               slope_group=rep(.5,nConditions), 
               midpoint_subject=rep(1,nSubjects), 
               slope_subject=rep(.5,nSubjects)) 
   init2<-list(midpoint_group=rep(.5, nConditions),
               slope_group=rep(.1,nConditions), 
               midpoint_subject=rep(.5,nSubjects), 
               slope_subject=rep(.1,nSubjects))
   init3<-list(midpoint_group=rep(.5,nConditions), 
               slope_group=rep(.1,nConditions), 
               midpoint_subject=rep(.5,nSubjects), 
               slope_subject=rep(.1,nSubjects))
   inits <- list(init1, init2, init3)
   output.stan <- stan(file="depletion.stan", iter=4000, data=datalist, init=inits,chains=3) 
   return(output.stan)
}


depletionPct <- function(trial, midpoint, spread) {
   dep <- 1- exp(-exp(2 * spread * midpoint / log(2) * ( log(trial) - log(midpoint) ) + log(log(2)) ) )
   return(dep)
}

pIER <- function(trial, midpoint, spread, initialDepletion, maxDepletion) { 
   dep <- depletionPct(trial, midpoint, spread)
   return(initialDepletion + (maxDepletion-initialDepletion)*dep)
}


pcorrect <- function(trial, midpoint, spread, initialDepletion, maxDepletion, chanceCorrect) {
   ier <- pIER(trial, midpoint, spread, initialDepletion, maxDepletion)
   pCorrect <- 1-ier  + ier * chanceCorrect
   return(pCorrect)
}

##  Posterior Group Plot ##
plotWeibullPosterior <- function(output.stan, factorName="Fatigue", chanceCorrect = 2/7) {
   filename <- paste(factorName, ".png", sep="")

   output.post <- extract(output.stan, c("midpoint_group", "slope_group"), 
                           permute=TRUE)

   samps <- sample(1:3000, 100)
   colvec1 <- rgb(matrix(c(.5, .5, .7), 3000, 3, byrow=TRUE)  + 
                  matrix( rnorm(3*3000, 0, .05), 3000, 3))
   colvec2 <- rgb(matrix(c(.7, .5, .5), 3000, 3, byrow=TRUE)  + 
                  matrix( rnorm(3*3000, 0, .05), 3000, 3))
   tr <- (1:500)/500
   
   
   png(filename, 1600, 800, pointsize=24)
   #png(filename, 1200, 400)
   #setEPS()
   #postscript(filename, 
   par(mfrow=c(1,2))
   #plot(c(0,500), c(0,1), type='n', main=paste("Posterior Group Mean", factorName, sep="\n"), xlab="Item", ylab=factorName)
   plot(c(0,500), c(0,1), type='n', main=paste("Posterior Group Mean", "Latent Depletion", sep="\n"), xlab="Item", ylab=factorName)
   for ( s in samps ) { 
      lines(500*tr, depletionPct(tr, midpoint=output.post$midpoint_group[s,1], spread=output.post$slope_group[s,1]), col=colvec1[s])
      lines(500*tr, depletionPct(tr, midpoint=output.post$midpoint_group[s,2], spread=output.post$slope_group[s,2]), col=colvec2[s])
   }
   lines(500*tr, depletionPct(tr, midpoint=mean(output.post$midpoint_group[,1]), spread=mean(output.post$slope_group[,1])), col='blue', lwd=2)
   lines(500*tr, depletionPct(tr, midpoint=mean(output.post$midpoint_group[,2]), spread=mean(output.post$slope_group[,2])), col='red', lwd=2)
   legend('topleft', c("Control", "Warning"), lty=1, col=c('blue','red'))
   
   
   #plot(c(0,500), c(0,1), type='n', main=paste("Posterior Group Mean\nProbability of", factorName, "Response", sep=" "), xlab="Item", ylab="P(Response)")
   #for ( s in samps ) { 
   #   lines(500*tr, pIER(tr, midpoint=output.post$midpoint_group[s,1], spread=output.post$slope_group[s,1], initialDepletion=0, maxDepletion=1), col=colvec1[s])
   #   lines(500*tr, pIER(tr, midpoint=output.post$midpoint_group[s,2], spread=output.post$slope_group[s,2], initialDepletion=0, maxDepletion=1), col=colvec2[s])
   #}
   #lines(500*tr, pIER(tr, midpoint=mean(output.post$midpoint_group[,1]), spread=mean(output.post$slope_group[,1]), initialDepletion=0, maxDepletion=1), col='blue', lwd=2)
   #lines(500*tr, pIER(tr, midpoint=mean(output.post$midpoint_group[,2]), spread=mean(output.post$slope_group[,2]), initialDepletion=0, maxDepletion=1), col='red', lwd=2)
   #
   
   plot(c(0,500), c(0,1), type='n', main="Posterior Group Mean\nProbability of Endorsing Response", xlab="Item", ylab="P(Endorse)")
   for ( s in samps ) { 
      lines(500*tr, 1-pcorrect(tr, midpoint=output.post$midpoint_group[s,1], spread=output.post$slope_group[s,1], initialDepletion=0, maxDepletion=1, chanceCorrect=chanceCorrect), col=colvec1[s])
      lines(500*tr, 1-pcorrect(tr, midpoint=output.post$midpoint_group[s,2], spread=output.post$slope_group[s,2], initialDepletion=0, maxDepletion=1, chanceCorrect=chanceCorrect), col=colvec2[s])
   }
   lines(500*tr, 1-pcorrect(tr, midpoint=mean(output.post$midpoint_group[,1]), spread=mean(output.post$slope_group[,1]), initialDepletion=0, maxDepletion=1, chanceCorrect=chanceCorrect), col='blue', lwd=2)
   lines(500*tr, 1-pcorrect(tr, midpoint=mean(output.post$midpoint_group[,2]), spread=mean(output.post$slope_group[,2]), initialDepletion=0, maxDepletion=1, chanceCorrect=chanceCorrect), col='red', lwd=2)
   
   dev.off()
}




calcLongString <- function(dat) { 
   subjects <- levels(dat$Subject)

   lsmax <- rep(NA, length(subjects))
   lsavg <- rep(NA, length(subjects))
   cond.summary <- rep(NA, length(subjects))
   cond.full <- c()
   lsall <- c()
   for (sn in 1:length(subjects)) { 
      sj <- subjects[sn]
      ix <- sort(dat$Order[dat$Subject==sj], index=TRUE)$ix
      xx <- dat$Response[dat$Subject==sj]
      xx <- xx[ix]

      lsvec <- rep(NA, 25) 

      for (pg in 1:25 ) { 
         startx <- (pg-1)*20 + 1
         endx <- pg*20
         lsvec[pg] <- max( rle(xx[startx:endx])$lengths) 
      }
      cond <- unique(dat$Condition[dat$Subject==sj])
      cond.summary[sn] <- cond
      cond.full <- c(cond.full, rep(cond, 25))
      lsall <- c(lsall, lsvec)
      lsmax[sn] <- max(lsvec)
      lsavg[sn] <- mean(lsvec)
   }
   cond.summary <- factor(cond.summary)
   cond.full <- factor(cond.full)
   rval1 <- data.frame(Subject=subjects, Condition=cond.summary, 
                        LS_max=lsmax, LS_avg=lsavg)
   rval2 <- data.frame(Subject=rep(subjects, each=25), Condition=cond.full,
                        Page=rep(1:25, length(subjects)), LongString=lsall)
   return(list(Summary=rval1, Full=rval2))
}

calcSumInfrequency <- function(infreq) { 
   subjects <- levels(infreq$Subject)
   suminfreq <- rep(NA, length(subjects))
   for (sn in 1:length(subjects)) { 
      suminfreq[sn] <- with(infreq, sum(Response[Subject==subjects[sn]]))
   }
   rval <- data.frame(Subject=subjects, Sum_infreq=suminfreq)
   return(rval)
}

calcMeanPageTime <- function(datrt) { 
   subjects <- levels(datrt$Subject)
   meanrt <- rep(NA, length(subjects))
   for (sn in 1:length(subjects)) { 
      meanrt[sn] <- with(datrt, mean(RT[Subject==subjects[sn]]))
   }
   rval <- data.frame(Subject=subjects, MeanPageTime=meanrt)
   return(rval)
}



findMaxTrials<-function (accuracy, pctsubjects, midpoint, slope, sdmid, sdslope) { 
   lowermid <- qnorm(1-pctsubjects, midpoint, sdmid)
   upperslp <- qnorm(pctsubjects, slope, sdslope)
   erf <- function(tr) {
      (accuracy - pcorrect(tr, lowermid, upperslp, 0, 1, 1/7))^2
   }

   tr <- optimize(erf, c(0,3))
   print(tr)
   return(tr$minimum * 500)
}
