rm(list=ls())
require(rstan)
require(lme4)
require(lmerTest)
source("analyzeParsed_fn.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Indicator for whether to re-fit Bayesian latent carelessness model
fitNew <- FALSE # TRUE # 

Exp <- 1;
inData <- read.csv(paste("ParsedData_", Exp, ".csv", sep=""))

if (Exp == 2) { 
  inData$Proctor <- NA
  inDataProctor <- read.csv("../Experiment 2/IP Study 2 Data.csv")
  for (subject in levels(inData$Subject)) { 
    inData$Proctor[inData$Subject == subject] <- with(inDataProctor, 
      levels(proctor)[proctor[levels(id)[id] == subject]])
  }
  inData$Proctor <- factor(inData$Proctor)
}


# Remove subjects with missing data
badSubjects <- unique(subset(inData, is.na(Response) | 
                                     Response==-99)$Subject)
for (sj in badSubjects) { 
  inData <- subset(inData, Subject!= sj)
}

inData$Subject <- factor(inData$Subject)
inData$Condition <- factor(inData$Condition)
subjects <- levels(inData$Subject)


#####################################
### Calculate long string indices ###
#####################################
longString <- calcLongString(inData)


######################
### Reverse Coding ###
######################
if (Exp == 1) { 
  reverseItems <- c("Q296", "Q11", "Q16", "Q110", "Q162", "Q527", "Q263", "Q325", "Q519", "Q505", "Q510", "Q330", "Q333")
} else if (Exp == 2) { 
  reverseItems <- c("Q292", "Q7", "Q12", "Q106", "Q158", "Q523", "Q259", "Q321", "Q515", "Q501", "Q503", "Q506", "Q326", "Q329")
}
for (item in reverseItems) { 
   inData$Response[inData$Question==item] <- 8-inData$Response[inData$Question==item]
}


#################################
### Infrequency Item Analysis ###
#################################
# 1: Incorrect; 0: Correct
infrequency <- subset(inData, Scale=="Infrequency")
infrequency$Response[infrequency$Response <= 2] <- 0
infrequency$Response[infrequency$Response > 2] <- 1
infrequency$Order <- infrequency$Order  / 500
if (Exp ==1 ) { 
  infrequency <- infrequency[,1:5]
} else {
  infrequency <- infrequency[,c(1:5, 7)]
}


sumInfrequency <- calcSumInfrequency(infrequency)

bins <- seq(0,400, 10)
prx <- matrix(NA, 381, length(bins))
for (sn in 1:length(subjects)) { 
   for ( i in 1:length(bins) ) {
      responses <- subset(infrequency, Subject==subjects[sn] & 
                                       (Order*500) > bins[i] & 
                                        Order <= (bins[i]+100) )
      prx[sn,i] <- mean(responses$Response, na.rm=TRUE)
   }
}

png(paste("InfrequencyRaw_", Exp, ".png", sep=""), 400, 400)
   matplot(bins, t(prx), type='l', main="Infrequency Items", 
            ylab="P(Incorrect)", xlab="Trial")
   lines(bins, apply(prx, 2, mean, na.rm=TRUE), lwd=3)
dev.off()
#xx <- apply(prx, 2, mean, na.rm=TRUE)


## Bayesian LGC anlaysis 
if (fitNew) { 
   fit.infrequency <- fitWeibullStan(infrequency)
   save(fit.infrequency, file=paste("fitinfrequency_", Exp, ".Rdata", sep=""))
} else { 
   load(paste("fitinfrequency_", Exp, ".Rdata", sep=""))
}


subject.mid <- extract(fit.infrequency, c("midpoint_subject"), permute=TRUE)$midpoint_subject
subject.mid.mn <- apply(subject.mid, 2, mean)

subject.slp <- extract(fit.infrequency, c("slope_subject"), permute=TRUE)$slope_subject
subject.slp.mn <- apply(subject.slp, 2, mean)

latentInfrequency <- data.frame(Subject=subjects, Midpoint=subject.mid.mn, Slope=subject.slp.mn)

latentDepletion <- rep(NA, dim(inData)[1])
attach(inData)
for (i in 1:dim(inData)[1]) { 
  latentDepletion[i] <- depletionPct(Order[i]/500, 
                                     midpoint=subject.mid.mn[Subject[i]],
                                     spread=subject.slp.mn[Subject[i]])
}
detach(inData)
inData$latentDepletion <- latentDepletion


cleanData <- inData
subjects <- levels(factor(infrequency$Subject))
for (sn in 1:length(subjects)) { 
   if(subject.mid.mn[sn] < 1) { 
      cleanData <- subset(cleanData, Subject != subjects[sn])
   }
}
cleanData$Subject <-factor(cleanData$Subject)


## NHST Analysis
infrequency <- subuset(ind
x1.infrequency <- glmer(Response ~ Order*Condition + (1|Subject), 
                        family=binomial, data=infrequency)
x2.infrequency <- glmer(Response ~ Order*Condition + (Order|Subject), 
                        family=binomial, data=infrequency)
anova(x1.infrequency, x2.infrequency)
exp(summary(x1.infrequency)$coefficients[2:4])


infrequency <- subset(inData, Scale=="Infrequency")
infrequency$Response[infrequency$Response <= 2] <- 0
infrequency$Response[infrequency$Response > 2] <- 1
infrequency$Order <- infrequency$Order  / 500
if (Exp ==1 ) { 
  infrequency <- infrequency[,1:5]
} else {
  infrequency <- infrequency[,c(1:5, 7:8)]
}
cl1 <- glmer(Response ~ latentDepletion + (1|Subject), family=binomial, 
               data=infrequency)
clc1 <- glmer(Response ~ latentDepletion*Condition + (1|Subject), family=binomial, 
               data=infrequency)
cc1 <- glmer(Response ~ Condition + (1|Subject), family=binomial, 
               data=infrequency)


## Check for Proctor effect on careless responding
if (Exp == 2) { 
  proctor.effect <- glmer(Response ~ Proctor + (1|Subject), 
                          family=binomial, data=infrequency)
  summary(proctor.effect)
}

#################################
#################################


###########################################
### State Positive Affectivity Analysis ###
###########################################
positive <- subset(inData, Scale=="Positive")
positive$Order <- positive$Order  / 500

poc1 <- lmer(Response ~ Order*Condition + (1|Subject), 
              data=positive)
po1 <- lmer(Response ~ Order + (1|Subject), 
              data=positive)
pl1 <- lmer(Response ~ latentDepletion + (1|Subject), 
              data=positive)
plc1 <- lmer(Response ~ latentDepletion*Condition + (1|Subject), 
              data=positive)
plo1 <- lmer(Response ~ Order + latentDepletion+ (1|Subject), 
              data=positive)
ploc1 <- lmer(Response ~ Order*Condition + latentDepletion*Condition + 
              (1|Subject), 
              data=positive)
p2 <- lmer(Response ~ Order*Condition + (Order|Subject), 
              data=positive)
summary(ploc1)
anova(plc1, ploc1)
anova(poc1, ploc1)
###########################################
###########################################



###########################################
### State Negative Affectivity Analysis ###
###########################################
negative <- subset(inData, Scale=="Negative")
negative$Order <- negative$Order  / 500

noc1 <- lmer(Response ~ Order*Condition + (1|Subject), 
              data=negative)
no1 <- lmer(Response ~ Order + (1|Subject), 
              data=negative)
nl1 <- lmer(Response ~ latentDepletion + (1|Subject), 
              data=negative)
nlc1 <- lmer(Response ~ latentDepletion*Condition + (1|Subject), 
              data=negative)
nlo1 <- lmer(Response ~ Order + latentDepletion+ (1|Subject), 
              data=negative)
nloc1 <- lmer(Response ~ Order*Condition + latentDepletion*Condition+ (1|Subject), 
              data=negative)
n2 <- lmer(Response ~ Order*Condition + (Order|Subject), 
              data=negative)
summary(nloc1)
anova(nlc1, nloc1)
anova(noc1, nloc1)
###########################################
###########################################



###########################################
### State Frustration Analysis ###
###########################################
frustration <- subset(inData, Scale=="Frustration")
frustration$Order <- frustration$Order  / 500

foc1 <- lmer(Response ~ Order*Condition + (1|Subject), 
              data=frustration)
fo1 <- lmer(Response ~ Order + (1|Subject), 
              data=frustration)
fl1 <- lmer(Response ~ latentDepletion + (1|Subject), 
              data=frustration)
flc1 <- lmer(Response ~ latentDepletion*Condition + (1|Subject), 
              data=frustration)
flo1 <- lmer(Response ~ Order + latentDepletion+ (1|Subject), 
              data=frustration)
floc1 <- lmer(Response ~ Order*Condition + latentDepletion*Condition+ 
              (1|Subject), 
              data=frustration)
f2 <- lmer(Response ~ Order*Condition + (Order|Subject), 
              data=frustration)

summary(floc1)
anova(flc1, floc1)
anova(foc1, floc1)
###########################################
###########################################


###########################################
### State Self Efficacy Analysis ###
###########################################
self.efficacy <- subset(inData, Scale=="Self_Efficacy")
#self.efficacy$Response[self.efficacy$Response <= 2] <- 0
#self.efficacy$Response[self.efficacy$Response > 0] <- 1
self.efficacy$Order <- self.efficacy$Order  / 500

eoc1 <- lmer(Response ~ Order*Condition + (1|Subject), 
              data=self.efficacy)
eo1 <- lmer(Response ~ Order + (1|Subject), 
              data=self.efficacy)
el1 <- lmer(Response ~ latentDepletion + (1|Subject), 
              data=self.efficacy)
elc1 <- lmer(Response ~ latentDepletion*Condition + (1|Subject), 
              data=self.efficacy)
elo1 <- lmer(Response ~ Order + latentDepletion+ (1|Subject), 
              data=self.efficacy)
eloc1 <- lmer(Response ~ Order*Condition + latentDepletion*Condition + 
              (1|Subject), 
              data=self.efficacy)
e2 <- lmer(Response ~ Order*Condition + (Order|Subject), 
              data=self.efficacy)
summary(eloc1)
anova(elc1, eloc1)
anova(eoc1, eloc1)
###########################################
###########################################



################################
### State Depletion Analysis ###
################################
depletionData <- subset(inData, Scale=="Depletion")
depletionData$Order <- depletionData$Order  / 500
pois <- FALSE

# Latent depletion analysis
if (pois) { 
x <- glmer(Response ~ (1|Subject), family=poisson, data=depletionData)
yq <- glmer(Response ~ Order + (1|Subject), family=poisson, 
               data=depletionData)
yl <- glmer(Response ~ latentDepletion + (1|Subject), family=poisson, 
               data=depletionData)
yql <- glmer(Response ~ Order + latentDepletion + (1|Subject), 
               family=poisson, data=depletionData)
} else { 

x1.statedepletion <- lmer(Response ~ Order*Condition + (1|Subject), 
                          data=depletionData)
x2.statedepletion <- lmer(Response ~ Order*Condition + (Order|Subject), 
                          data=depletionData)

yq <- lmer(Response ~ Order + (1|Subject), data=depletionData)
yl <- lmer(Response ~ latentDepletion + (1|Subject), data=depletionData)
yql <- lmer(Response ~ Order + latentDepletion + (1|Subject), 
            data=depletionData)
}

anova(x,yq)
anova(x,yl)
anova(yl,yql)
anova(yq,yql)



#  Latent Depletion * Condition Analysis
if (pois) { 
x1 <- glmer(Response ~ Condition + (1|Subject), family=poisson, 
            data=depletionData)
yq1 <- glmer(Response ~ Condition*Order + (1|Subject), family=poisson, 
               data=depletionData)
yl1 <- glmer(Response ~ Condition*latentDepletion + (1|Subject), 
               family=poisson, data=depletionData)
yql1 <- glmer(Response ~ Condition * Order + Condition * latentDepletion 
            + (1|Subject), family=poisson, data=depletionData)
} else { 
x1 <- lmer(Response ~ Condition + (1|Subject), 
            data=depletionData)
yq1 <- lmer(Response ~ Condition*Order + (1|Subject),
               data=depletionData)
yl1 <- lmer(Response ~ Condition*latentDepletion + (1|Subject), 
               data=depletionData)
yql1 <- lmer(Response ~ Condition * Order + Condition * latentDepletion 
            + (1|Subject), data=depletionData)
}

anova(x1,yq1)
anova(x1,yl1)
anova(yl1,yql1)
anova(yq1,yql1)


exp(summary(yq1)$coefficients[2:4])


##  With cleaned data ##
depletionClean <- subset(cleanData, Scale=="Depletion")
depletionClean$Subject <- factor(depletionClean$Subject)
depletionClean$Order <- depletionClean$Order  / 500

# Latent depletion analysis
#  Latent Depletion * Condition Analysis
if (pois) { 
x2 <- glmer(Response ~(1|Subject), family=poisson, 
            data=depletionClean)
yq2 <- glmer(Response ~ Condition*Order + (1|Subject), family=poisson, 
               data=depletionClean)
yl2 <- glmer(Response ~ Condition*latentDepletion + (1|Subject), 
               family=poisson, data=depletionClean)
yql2 <- glmer(Response ~ Condition * Order + Condition * latentDepletion 
            + (1|Subject), family=poisson, data=depletionClean)
} else { 
   x1.depletionClean <- lmer(Response ~ Order * Condition + (1|Subject), 
                             data=depletionClean)
   x2.depletionClean <- lmer(Response ~ Order * Condition + (Order|Subject),
                             data=depletionClean)

   yq2 <- lmer(Response ~ Condition*Order + (1|Subject), 
                  data=depletionClean)
   yl2 <- lmer(Response ~ Condition*latentDepletion + (1|Subject), 
                  data=depletionClean)
   yql2 <- lmer(Response ~ Condition * Order + Condition * latentDepletion 
               + (1|Subject), data=depletionClean)
}

anova(x2,yq2)
anova(x2,yl2)
anova(yl2,yql2)
anova(yq2,yql2)

################################
################################



########################
### Fatigue Analysis ###
########################
fatigue <- subset(inData, Scale=="Fatigue")
fatigue$Order <- fatigue$Order  / 500

fgoc1 <- lmer(Response ~ Order*Condition + (1|Subject), 
              data=fatigue)
fgo1 <- lmer(Response ~ Order + (1|Subject), 
              data=fatigue)
fgl1 <- lmer(Response ~ latentDepletion + (1|Subject), 
              data=fatigue)
fglo1 <- lmer(Response ~ Order + latentDepletion+ (1|Subject), 
              data=fatigue)
fg2 <- lmer(Response ~ Order*Condition + (Order|Subject), 
              data=fatigue)

########################
########################



##########################
### Page Time Analysis ###
##########################
inDataRT <- read.csv(paste("ParsedDataRT_", Exp, ".csv", sep=""))
for (sj in badSubjects) { 
 inDataRT <- subset(inDataRT, Subject!= sj)
}

inDataRT$Subject <- factor(inDataRT$Subject)
if (Exp == 1) { 
inDataRT$Question <- as.numeric(gsub("Page", "", 
                              levels(inDataRT$Question))[inDataRT$Question])
} else if (Exp == 2) { 
inDataRT$Question <- as.numeric(gsub("Time.", "", 
                              levels(inDataRT$Question))[inDataRT$Question])
}
inDataRT$Order <- inDataRT$Question / 25
meanPageTime <- calcMeanPageTime(inDataRT)

#badsubj <- with(inDataRT, unique(Subject[is.na(RT)]))
#for (sj in badsubj) { 
#   inDataRT <- subset(inDataRT, Subject != sj)
#}
#inDataRT$Subject <- factor(inDataRT$Subject)

## 20 Items per page
## Trial scale is out of 500
## Take midpoint of every page for estimated depletion
tr <- seq(10,490, by=20) / 500

latentDepletion <- rep(NA, dim(inDataRT)[1])
attach(inDataRT)
for (i in 1:dim(inDataRT)[1]) { 
  latentDepletion[i] <- depletionPct(tr[Question[i]], 
                                     midpoint=subject.mid.mn[Subject[i]],
                                     spread=subject.slp.mn[Subject[i]])
}
detach(inDataRT)
inDataRT$latentDepletion <- latentDepletion
rm(latentDepletion)


x <- lmer(log(RT) ~ (1|Subject), data=inDataRT)
yq <- lmer(log(RT) ~ Order + (1|Subject), data=inDataRT)
yl <- lmer(log(RT) ~ latentDepletion + (1|Subject), data=inDataRT)
yql <- lmer(log(RT) ~ Order + latentDepletion + (1|Subject), data=inDataRT)
yql.b <- lmer(log(RT) ~ Order * latentDepletion + (1|Subject), data=inDataRT)

anova(x,yq)
anova(x,yl)
anova(yl,yql)
anova(yq,yql)

x1.pagetime <- lmer(log(RT) ~ Condition * Order + (1|Subject), 
                    data=inDataRT)
x2.pagetime <- lmer(log(RT) ~ Condition * Order + (Order|Subject),
                    data=inDataRT)
anova(x1.pagetime, x2.pagetime)


yl1 <- lmer(log(RT) ~ Condition * latentDepletion + (1|Subject), data=inDataRT)
yql1 <- lmer(log(RT) ~ Order + latentDepletion + (1|Subject), data=inDataRT)

anova(x,yq1)
anova(x,yl1)
anova(yl,yql1)
anova(yq,yql1)
##########################
########################## 



############################
### Long String Analysis ###
############################
longString$Full$latentDepletion <- NA
pgseq <- round(seq(0.04,0.96, by=0.04), 3)
longString$Full$Order <- longString$Full$Page/25 - 1/50
for (sj in subjects) {
   for (pg in  pgseq) {
   #for (pg in 1:25) {
      if (Exp ==1  & (sj == "102532" | sj=="103870" | sj=="103189")) { 
        next 
      }
#      if (Exp ==2  & (sj== "R_0AKq2s7rVVgZwmp")) { 
#        next 
#      }
      longString$Full$latentDepletion[longString$Full$Subject==sj & 
               longString$Full$Page==pg*25] <- 
               with(inDataRT, latentDepletion[Subject==sj & Order==pg])
               
            
   }
}
longStringLatent <- subset(longString$Full, !is.na(latentDepletion))

longStringLatent$PageN <- (longStringLatent$Page - 
                           mean(longStringLatent$Page)) / 
                           sd(longStringLatent$Page)

xl <- glmer(LongString ~ (1|Subject), family=poisson, longStringLatent)
yl <- glmer(LongString ~ latentDepletion + (1|Subject), family=poisson, 
            longStringLatent)
yq  <- glmer(LongString ~ PageN + (1|Subject), family=poisson, 
            longStringLatent)
yql <- glmer(LongString ~ PageN + latentDepletion + (1|Subject), family=poisson, 
            longStringLatent)

yq1 <- glmer(LongString ~ Page*Condition + (1|Subject), family=poisson, longString$Full)


x1.longstring <- glmer(LongString ~ Order*Condition + (1|Subject), 
                       family=poisson, longString$Full)
x2.longstring <- glmer(LongString ~ Order*Condition + (Page|Subject), 
                       family=poisson, longString$Full)
anova(x1.longString, x2.longString)



exp(summary(yq)$coefficients[2])/ exp(summary(yq)$coefficients[1])

############################
############################




withinTable <- cbind(sumInfrequency, meanPageTime[,2], 
                     longString$Summary[,c("LS_max", "LS_avg")],
                     latentInfrequency[,2:3])
names(withinTable)[3] <- "MeanPageTime"
round(cor(withinTable[,2:7], use="complete.obs"), 2)
attach(withinTable)
cor.test(Sum_infreq, MeanPageTime)
cor.test(Sum_infreq, LS_max)
cor.test(Sum_infreq, LS_avg)
cor.test(Sum_infreq, Midpoint)
cor.test(Sum_infreq, Slope)
cor.test(MeanPageTime, LS_max)
cor.test(MeanPageTime, LS_avg)
cor.test(LS_max, LS_avg)
detach(withinTable)



####################################################
########  Between Experiment Comparison ############
####################################################
infrequency <- subset(inData_all, Scale=="Infrequency")
infrequency$Response[infrequency$Response <= 2] <- 0
infrequency$Response[infrequency$Response > 2] <- 1
infrequency$Order <- infrequency$Order  / 500
infrequency <- infrequency[,c(1:5,8)]

x1.infrequency <- glmer(Response ~ Order*Experiment + (1|Subject),
                        family=binomial, data=infrequency)
exp(summary(x1.infrequency)$coefficients[2:4])


load(paste("fitinfrequency_1.Rdata", sep=""))
samples.exp1 <- extract(fit.infrequency, c("midpoint_group", "slope_group"), permute=TRUE)
load(paste("fitinfrequency_2.Rdata", sep=""))
samples.exp2 <- extract(fit.infrequency, c("midpoint_group", "slope_group"), permute=TRUE)

png("BetweenComparison.png", 1200, 400)
par(mfrow=c(1,2))
vioplot(500 * samples.exp1$midpoint_group[,1], 
        500 * samples.exp1$midpoint_group[,2], 
        500 * samples.exp2$midpoint_group[,1], 
        500 * samples.exp2$midpoint_group[,2], 
        col=grey(.8), names=c("Online\u2013Control", "Online\u2013Warning", 
                              "Lab\u2013Control", "Lab\u2013Warning"))
        mtext("Number of Items", 2, line=2)
title("Number of items to reach 50% latent depletion")
abline(h=1)

vioplot(1/500 * samples.exp1$slope_group[,1], 
        1/500 * samples.exp1$slope_group[,2], 
        1/500 * samples.exp2$slope_group[,1], 
        1/500 * samples.exp2$slope_group[,2], 
        col=grey(.8), names=c("Online\u2013Control", "Online\u2013Warning", 
                              "Lab\u2013Control", "Lab\u2013Warning"))
        mtext("Increase in Probability per Item", 2, line=2)
title("Rate of increase in latent depletion")
dev.off()



infrequency$Location <- 0
infrequency$Location[infrequency$Condition == "2" 
                     | infrequency$Condition == "3"] <- 1
infrequency$Location <- factor(infrequency$Location, 
                               labels=c("Online", "Lab"))

infrequency$Warning <- 0
infrequency$Warning[infrequency$Condition =="2" 
                     | infrequency$Condition =="4"] <- 1
infrequency$Warning <- factor(infrequency$Warning, 
                               labels=c("Control", "Warning"))


x1.infrequency <- glmer(Response ~ Location + (1|Subject), 
                        family=binomial, data=infrequency)
x2.infrequency <- glmer(Response ~ Order*Condition + (Order|Subject), 
                        family=binomial, data=infrequency)
anova(x1.infrequency, x2.infrequency)
exp(summary(x1.infrequency)$coefficients[2:4])



### Between RT
inDataRT_1 <- read.csv(paste("ParsedDataRT_1.csv", sep=""))
inDataRT_2 <- read.csv(paste("ParsedDataRT_2.csv", sep=""))

inDataRT_1$Question <- as.numeric(gsub("Page", "", 
                         levels(inDataRT_1$Question))[inDataRT_1$Question])
inDataRT_2$Question <- as.numeric(gsub("Time.", "", 
                         levels(inDataRT_2$Question))[inDataRT_2$Question])
inDataRT_1$Location <- 1
inDataRT_2$Location <- 0
inDataRT_all <- rbind(inDataRT_1, inDataRT_2)

inDataRT_all$Order <- inDataRT_all$Question / 25
#meanPageTime <- calcMeanPageTime(inDataRT_all)

#badsubj <- with(inDataRT, unique(Subject[is.na(RT)]))
#for (sj in badsubj) { 
#   inDataRT <- subset(inDataRT, Subject != sj)
#}
inDataRT_all$Subject <- factor(inDataRT_all$Subject)
inDataRT_all$Location <- factor(inDataRT_all$Location, 
                                 labels=c("Online", "Lab"))


rt1.location <- lmer(log(RT) ~ Location*Order + (1|Subject), 
                     data=inDataRT_all)
rt2.location <- lmer(log(RT) ~ Order*Location + (Order|Subject), 
                     data=inDataRT_all)




## LongString

longString_1$Full$Location <- 0
longString_2$Full$Location <- 1
longString_all <- rbind(longString_1$Full, longString_2$Full)
pgseq <- round(seq(0.04,0.96, by=0.04), 3)
longString_all$Order <- longString_all$Page/25 - 1/50


longString_all$PageN <- (longString_all$Page - 
                           mean(longString_all$Page)) / 
                           sd(longString_all$Page)


longString_all$Location <- factor(longString_all$Location, 
                                 labels=c("Online", "Lab"))

x1.longstring <- glmer(LongString ~ Order*Location + (1|Subject), 
                       family=poisson, longString_all)
x2.longstring <- glmer(LongString ~ Order*Location + (Page|Subject), 
                       family=poisson, longString_all)
anova(x1.longString, x2.longString)



### Number of trials...

# Experiment 1
#findMaxTrials(.8, .9, 1.09, .72, .1, .1)
#findMaxTrials(.8, .9, 1.48, .65, .1, .1)


#findMaxTrials(.8, .9, 1.28, .59, .1, .1)
#findMaxTrials(.8, .9, 1.58, .51, .1, .1)
