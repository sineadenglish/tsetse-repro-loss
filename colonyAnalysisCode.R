#############################################################
#
# colonyAnalysisCode.r => analysis code for paper 
# English et al. 2022 "Investigating the unaccounted ones: insights on age-dependent reproductive loss in a viviparous fly"
#
# last updated 29/9/2022
# author: Sinead English (sinead.english@bristol.ac.uk)
#
#############################################################

library(mgcv)
library(AICcmodavg)
library(ggplot2)

#setwd("[ADD FILE LOCATION HERE]")
dat <- read.csv("colonyAbortions.csv")

# data processing
dat$femAged <- dat$femAge*7                         # female age in days
dat$totAbort <- dat$Eggs + dat$L1 + dat$L2 + dat$L3 # total abortions recorded per tray
dat$totRepro <- dat$totAbort + dat$pupNo            # total number of offspring produced (pupae + abortions) per tray

# for plots:
dat$Egg.propPerTot <- dat$Eggs/dat$totRepro
dat$L1.propPerTot <- dat$L1/dat$totRepro
dat$L2.propPerTot <- dat$L2/dat$totRepro
dat$L3.propPerTot <- dat$L3/dat$totRepro


# GAM analysis - 

# EGGS
#**********************choose knots***************
eggGAMfits <- lapply(3:10,function(x){
  g1Gam <- gam(cbind(dat$Eggs,dat$totRepro-dat$Eggs) ~ s(dat$femAged,bs="cr",k=x),family=binomial,data=dat,se=T,method="ML")
  return(c(x,AICc(g1Gam,nobs=length(dat[,1]))))
})

eggGAMfits <- do.call(rbind.data.frame,eggGAMfits)
names(eggGAMfits) <- c("Number of knots","AICc"); eggGAMfits
# k = 7 
#************************************************

g1Gam <- gam(cbind(dat$Eggs,dat$totRepro-dat$Eggs) ~ s(dat$femAged,k=7),family=binomial,data=dat,se=T,method="ML")
summary(g1Gam)
gam.check(g1Gam)

g1M2 <- predict(g1Gam
                ,level=1
                ,type="response",se=T,unconditional=T)

dat$Egg.pred <- g1M2$fit
dat$Egg.se <- g1M2$se.fit

# L1
#**********************choose knots***************
l1GAMfits <- lapply(3:10,function(x){
  g2Gam <- gam(cbind(dat$L1,dat$totRepro-dat$L1) ~ s(dat$femAged,bs="cr",k=x),family=binomial,data=dat,se=T,method="ML")
  return(c(x,AICc(g2Gam,nobs=length(dat[,1]))))
})

l1GAMfits <- do.call(rbind.data.frame,l1GAMfits)
names(l1GAMfits) <- c("Number of knots","AICc"); l1GAMfits
# k = 6
#************************************************

g2Gam <- gam(cbind(dat$L1,dat$totRepro-dat$L1) ~ s(dat$femAged,bs="cr",k=6),family=binomial,data=dat,se=T,method="ML")
summary(g2Gam)
gam.check(g2Gam)

g2M2 <- predict(g2Gam
                ,level=1
                ,type="response",se=T,unconditional=T)

dat$L1.pred <- g2M2$fit
dat$L1.se <- g2M2$se.fit

# L2
#**********************choose knots***************
l2GAMfits <- lapply(3:10,function(x){
  g3Gam <- gam(cbind(dat$L2,dat$totRepro-dat$L2) ~ s(dat$femAged,bs="cr",k=x),family=binomial,data=dat,se=T,method="ML")
  return(c(x,AICc(g3Gam,nobs=length(dat[,1]))))
})

l2GAMfits <- do.call(rbind.data.frame,l2GAMfits)
names(l2GAMfits) <- c("Number of knots","AICc"); l2GAMfits
# k = 4
#************************************************

g3Gam <- gam(cbind(dat$L2,dat$totRepro-dat$L2) ~ s(dat$femAged,bs="cr",k=4),family=binomial,data=dat,se=T,method="ML")
summary(g3Gam)

g3M2 <- predict(g3Gam
                ,level=1
                ,type="response",se=T,unconditional=T)

dat$L2.pred <- g3M2$fit
dat$L2.se <- g3M2$se.fit


# L3
#**********************choose knots***************
l3GAMfits <- lapply(3:10,function(x){
  g4Gam <- gam(cbind(dat$L3,dat$totRepro-dat$L3) ~ s(dat$femAged,bs="cr",k=x),family=binomial,data=dat,se=T,method="ML")
  return(c(x,AICc(g4Gam,nobs=length(dat[,1]))))
})

l3GAMfits <- do.call(rbind.data.frame,l3GAMfits)
names(l3GAMfits) <- c("Number of knots","AICc"); l3GAMfits
# k = 6
#************************************************

g4Gam <- gam(cbind(dat$L3,dat$totRepro-dat$L3) ~ s(dat$femAged,bs="cr",k=6),family=binomial,data=dat,se=T,method="ML")
summary(g4Gam)
gam.check(g4Gam)

g4M2 <- predict(g4Gam
                ,level=1
                ,type="response",se=T,unconditional=T)

dat$L3.pred <- g4M2$fit
dat$L3.se <- g4M2$se.fit

############

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
scale_fill_manual(values=cbPalette)
scale_colour_manual(values=cbPalette)

dat <- dat[order(dat$femAged),]
as.matrix(names(dat))

plotDat1 <- reshape(dat[,c('rowID','femAged','Egg.pred','Egg.se','L1.pred','L1.se','L2.pred','L2.se','L3.pred','L3.se')],
                    direction='long',
                    varying=c('Egg.pred','Egg.se','L1.pred','L1.se','L2.pred','L2.se','L3.pred','L3.se'),
                    timevar='stage',
                    times=c('Egg','L1','L2','L3'),
                    v.names=c('pred','se'),
                    idvar='rowID')

plotDat2 <- reshape(dat[,c('rowID','femAged','Egg.propPerTot','L1.propPerTot','L2.propPerTot','L3.propPerTot')],
                    direction='long',
                    varying=c('Egg.propPerTot','L1.propPerTot','L2.propPerTot','L3.propPerTot'),
                    timevar='stage',
                    times=c('Egg','L1','L2','L3'),
                    v.names=c('propPerTot'),
                    idvar='rowID')


plotDat <- merge(plotDat1, plotDat2, by=c("rowID","stage","femAged"))


#tiff("Fig2_abortionGAM.tiff", height = 5, width = 5, units = 'in', compression="lzw", res=400)
plotGAM <- ggplot(data=plotDat, aes(y=pred, x=femAged, group=stage, colour=stage,
                                    fill=stage)) +
  scale_fill_manual(values=cbPalette) + 
  scale_colour_manual(values=cbPalette) + 
  geom_line() +
  geom_ribbon(aes(ymin=pred-2*se, ymax=pred+2*se), alpha=.3, linetype=0) +
  geom_point(aes(x=femAged,y=propPerTot,col=stage),size=1) + 
  xlim(20,85) +
  ylim(0,0.5) +
  labs(  x="Mother age (days)"
         ,y="Proportion of offspring aborted") + 
  theme_classic()
plotGAM
#dev.off()
