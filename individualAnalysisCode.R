#############################################################
#
# individualAnalysisCode.r => analysis code for paper 
# English et al. 2022 "Investigating the unaccounted ones: insights on age-dependent reproductive loss in a viviparous fly"
#
# last updated 29/9/2022
# author: Sinead English (sinead.english@bristol.ac.uk)
#
#############################################################

library(binom)
library(lme4)
library(interactions)
library(car)
library(ggplot2)
library(gridExtra)

#setwd("[ADD FILE LOCATION HERE]")
indDat <- read.csv("indLarvipositionDat.csv")

############### SECTION 3 analysis ######################
# QUANTIFYING PRESUMED ABORTIONS IN INDIVIDUALLY TRACKED FEMALES

# First, use Hampel filter to establish thresholds

thresholdDat <- indDat[indDat$abortion==0,c("parity","name","timeBetLarv")]

firstLarv <- na.omit(thresholdDat$timeBetLarv[thresholdDat$parity=="first"])
lower_bound <- median(firstLarv) - 3 * mad(firstLarv, constant = 1);lower_bound  # 16
upper_bound <- median(firstLarv) + 3 * mad(firstLarv, constant = 1);upper_bound # 28

subLarv <- na.omit(thresholdDat$timeBetLarv[thresholdDat$parity=="subsequent"])
lower_bound <- median(subLarv) - 3 * mad(subLarv, constant = 1);lower_bound # 8
upper_bound <- median(subLarv) + 3 * mad(subLarv, constant = 1);upper_bound # 14


# Second, apply thresholds to create 'Previous birth outcome' code (column name = abortCode)

indDat$abortCode <- NA

# 0 = in 'normal range' (no abortion)
indDat$abortCode[indDat$parity=="first" & indDat$timeBetLarv >=16 & indDat$timeBetLarv <= 28] <- 0
indDat$abortCode[indDat$parity=="subsequent" & indDat$timeBetLarv >=8 & indDat$timeBetLarv <= 14] <- 0

# 1 = normal birth outcome but too-short gestation (data error potentially, not used in analysis)
indDat$abortCode[indDat$parity=="first" & indDat$timeBetLarv <16] <- 1
indDat$abortCode[indDat$parity=="subsequent" & indDat$timeBetLarv <8] <- 1

# 2 = presumed missed abortion (took too long to produce next offspring)
indDat$abortCode[indDat$parity=="first" & indDat$timeBetLarv>28] <- 2
indDat$abortCode[indDat$parity=="subsequent" & indDat$timeBetLarv>14] <- 2

# 3 = confirmed abortion (stage given)
indDat$abortCode[indDat$abortion==1] <- 3 # identified abortion

table(indDat$abortCode)

### for analysis, will want for each viable offspring what was previous outcome
### but add missed abortions for 'previous status' (actually current status)

# add column for previous code (only when ID is matched), make sure first sorted to order of ID then larvNum

indDat <- indDat[
  order(indDat[,"adults_id"], indDat[,"larvNum"]),
]

indDat$prevAdults_id <- c(NA,indDat$adults_id[-length(indDat$adults_id)])
indDat$prevAbortCode <- c(NA,indDat$abortCode[-length(indDat$abortCode)])

indDat$prevAbortCode[indDat$adults_id != indDat$prevAdults_id] <- NA

indDat$prevAbortCode[indDat$abortCode==2] <- 4

table(indDat$abortCode, indDat$prevAbortCode)

indDat$prevAbortCodeF <- NA
indDat$prevAbortCodeF[indDat$prevAbortCode==0] <- "1.PrevNormal"
indDat$prevAbortCodeF[indDat$prevAbortCode==3] <- "2.PrevConfirmedAbortion"
indDat$prevAbortCodeF[indDat$prevAbortCode==2] <- "3.PrevMissedAbortion"
indDat$prevAbortCodeF[indDat$prevAbortCode==4] <- "4.CurrMissedAbortion" # missed abortion in previous larviposition (i.e. potentially 2 births ago)

table(indDat$prevAbortCodeF, indDat$prevAbortCode)

# Third, look at pattern of assumed abortion with age

# only include those with known viable birth (0) or presumed abortion (2)
dat1 <- indDat[indDat$abortCode %in% c(0,2),]
dat1$CurrAbort <- ifelse(dat1$abortCode==0,0,1)

table(dat1$CurrAbort,dat1$abortCode)

m1 <- glm(CurrAbort ~ mAgeDays + I(mAgeDays^2) + name, data = dat1, family = "binomial")
summary(m1)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -3.3871476  0.7946162  -4.263 2.02e-05 ***
#   mAgeDays      -0.0005214  0.0294107  -0.018    0.986    
# I(mAgeDays^2)  0.0001618  0.0002535   0.638    0.523    
# namenutrition  0.9711233  0.2460567   3.947 7.92e-05 ***

# plot predicted results
newDat <- expand.grid(mAgeDays = 17:98, name = c("control","nutrition"))

newDat$fit <- predict(m1, newdata = newDat, se.fit = TRUE)$fit
newDat$se.fit <- predict(m1, newdata = newDat, se.fit = TRUE)$se.fit

# for points: bin average and proportion actual abortions
#***********calculate mean weight and confidence intervals for females per 10 days***********

dat1$mAgeBins <- cut(as.numeric(dat1$mAgeDays)
                     ,c(17,27,37,47,57,67,77,87,98)
                     ,labels=c("22","32","42","52","62","72","82","92"))
ctrl <- dat1[dat1$name %in% "control",]
nuts <- dat1[dat1$name %in% "nutrition",]

#************************************************************************************
ctrl.means <- lapply(unique(ctrl$mAgeBins),function(x){
  temp <- ctrl[ctrl$mAgeBins %in% x,]
  m <- binom.confint(sum(temp$CurrAbort),length(temp$CurrAbort),method="exact")
  m$mAgeBin <- x
  return(m)
})
ctrl.means <- do.call(rbind.data.frame,ctrl.means)

nuts.means <- lapply(unique(nuts$mAgeBins),function(x){
  temp <- nuts[nuts$mAgeBins %in% x,]
  m <- binom.confint(sum(temp$CurrAbort),length(temp$CurrAbort),method="exact")
  m$mAgeBin <- x
  return(m)
})
nuts.means <- do.call(rbind.data.frame,nuts.means)

nuts.means$name <- "nutrition"
ctrl.means$name <- "control"
abort.means <- rbind.data.frame(nuts.means,ctrl.means)

abort.means$mAgeBin <- as.numeric(as.character(abort.means$mAgeBin))
abort.means$name <- as.factor(abort.means$name)
levels(abort.means$name) <- c("control","nutrition")

names(abort.means)[4] <- "CurrAbort"
names(abort.means)[7] <- "mAgeDays"

cbPalette2 <- c("#E69F00", "#56B4E9")

#tiff("Fig3_presAbortGLM.tiff", height = 5, width = 5, units = 'in', compression="lzw", res=400)
predM1 <- ggplot()  + 
  scale_fill_manual(values=cbPalette2) + 
  scale_colour_manual(values=cbPalette2) + 
  geom_line(data=newDat, aes(x = mAgeDays, y = exp(fit), col=name)) + 
  geom_ribbon(data=newDat, aes(x = mAgeDays, ymax = exp(fit + 1.96 * se.fit),
                               ymin = exp(fit - 1.96 * se.fit),
                               group = name,
                               fill = name
  ), alpha=0.3,linetype=0) +
  geom_point(data=abort.means,aes(mAgeDays,CurrAbort,col=name)) +
  geom_errorbar(data=abort.means,aes(x=mAgeDays,ymin=lower,ymax=upper,col=name),width=.1) +
  labs(  x="Mother age (days)"
         ,y="Probability of presumed abortion") + 
  theme_classic()
predM1
#dev.off()


############### SECTION 3 analysis ######################

# select only those which produce viable offspring (abortCode==0) and where previous birth was normal, confirmed or presumed abortion
dat2 <- indDat[indDat$prevAbortCodeF %in% c("1.PrevNormal","2.PrevConfirmedAbortion","3.PrevMissedAbortion") & indDat$abortCode==0,]
dat2$mAgeDays_z <- (dat2$mAgeDays - mean(dat2$mAgeDays))/sd(dat2$mAgeDays)

# 1. How does time to next larva depend on previous birth outcome? 

durm1 <- lmer(timeBetLarv ~ prevAbortCodeF*name + mAgeDays_z + I(mAgeDays_z^2) + (1|adults_id), data=dat2)
summary(durm1)
anova(durm1)

# Look at means to show what is the difference in days 
tapply(dat2$timeBetLarv, list(dat2$name, dat2$prevAbortCodeF), mean)
# 
# 1.PrevNormal 2.PrevConfirmedAbortion 3.PrevMissedAbortion
# control       10.79187                11.81250             10.94118
# nutrition     10.95050                10.27273             11.16667

# 2. Effect of previous birth outcome on subsequent weight of offspring

pupm1 <- lmer(wet_weight ~ prevAbortCodeF*name + mAgeDays_z + I(mAgeDays_z^2) + (1|adults_id), data=dat2)
summary(pupm1)
anova(pupm1)

# Fig 4 based on lmer output
#tiff("Fig4A_larDur.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
p1 <- cat_plot(durm1, pred = name, modx = prevAbortCodeF, line.thickness=0.8,errorbar.width = 0.5,x.label="Experimental treatment",y.label = "Larviposition duration",
               plot.points=T, legend.main="Previous birth outcome")
p1
#dev.off()

#tiff("Fig4B_pupWeight.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
p2 <- cat_plot(pupm1, pred = name, modx = prevAbortCodeF, line.thickness=0.8,errorbar.width = 0.5,x.label="Experimental treatment",y.label = "Pupal weight (mg)",
               plot.points=T, legend.main="Previous birth outcome")
p2
#dev.off()

grid.arrange(p1,p2,ncol=2)


# 3. Effect of previous birth outcome on sex of offspring produced

sexm1 <- glmer(sexN ~ prevAbortCodeF*name + mAgeDays_z + I(mAgeDays_z^2) + (1|adults_id), family=binomial, data=dat2)
summary(sexm1)
anova(sexm1)
car::Anova(sexm1)