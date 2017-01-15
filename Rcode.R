################################################################################
# Updated version of the code for the analysis in:
#
#   "Water supply interruptions and suspected cholera incidence:
#     a time-series regression"
#   Aurelie Jeandron and colleagues
#   http://www.ag-myresearch.com/2015_jeandron_plosmed.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2015_jeandron_PlosMed_Rcodedata
################################################################################

# LOAD THE PACKAGES REQUIRED FOR THE ANALYSIS
library(splines)
library(tsModel)
library(dlnm)
library(gplots)

# SET THE OPTIONS
options(na.action="na.exclude")

# LOAD THE DATA
data <- read.csv("Dataset.csv")[,c("cholera","volume","date","rain","chlorine",
  "chol_highcov", "chol_lowcov")]

# PREPARE VARIABLES AND CUBIC SPLINE FUNCTION OF DATE
data$vol1000 <- (data$volume)/1000
data$date2<-as.Date(data$date, format="%d/%m/%Y", origin="30/12/1899")
data$dow <- weekdays(data$date2)
data <- data[order(data$date2),]
volcen <- quantile(data$vol1000, probs=c(0.95), na.rm=T)
spl <- bs(data$date2, degree=3, df=64)

################################################################################
# PLOTTING SERIES OF CHOLERA, CHLORINE, VOLUME and RAIN (Figure 1)
tiff(filename="figure1.tif", width=3, height=5, units="in",pointsize=8, res=300)
layout(matrix(c(1,2,3,4),4,1, byrow=TRUE), heights=c(1.5,1.5,1,1))
layout.show(4)
par(mar=c(2,6,2,2))
plot(data$date2, data$cholera, xlab="", ylab="Suspected cholera cases \n admitted to CTC",type="l")
plot(data$date2, data$vol1000, xlab="",ylab="Volume of water supplied \n(in 1,000 m3)", type="l")
plot(data$date2, data$chlorine, xlab="",ylab="Residual chlorine at \nplant output (in mg/l)", type="l")
plot(data$date2, data$rain, xlab="", ylab="Daily precipitation rate \n (in mm/h)", type="l")
par(mfrow=c(1,1))
dev.off()
par(mar=c(5,4,4,1)+0.1)

################################################################################
# PLOTTING MEAN INCIDENCE OF SUSPECTED CHOLERA AT LAG 6 BY QUINTILE OF VOLUME SUPPLIED AT LAG 0(Figure 2)
data$chol_lag6<-Lag(data$cholera, 6)
cutoffs<-round(quantile(data$vol1000, probs=0:5/5, na.rm=TRUE), digits=2)
data$volq<-cut(data$vol1000, breaks=cutoffs, include.lowest=TRUE)
tiff(filename="figure2.tif", width=5, height=3, units="in",pointsize=8, res=300)
par(mar=c(6,6,2,2))
plotmeans(chol_lag6~volq, data, ylim=c(2, 4),
  xlab=expression(paste("Quintiles of volume of water supplied at lag 0 (in 1,000 ",m^3,")")),
  ylab="Mean and 95% CI of the number of \nsuspected cholera cases at lag 6",
  barcol="black", legends=c("0 to 2.87","2.88 to 3.53","3.54 to 3.96","3.97 to 4.36","4.37 to 5.75"))
dev.off()
################################################################################
# MAIN MODEL WITH CENTERING OF VOLUME AT 95th PERCENTILE (volcen)

# CREATE THE CROSS-BASIS
model.basis<- crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=2, int=F))
# FIT THE MODEL WITH AUTOREGRESSIVE TERMS
model0 <- glm(cholera ~ model.basis+spl+dow+rain, data, family=quasipoisson)
dres0 <- residuals(model0, "deviance")
model0_ac <- update(model0,.~.+Lag(dres0,1:2))
# PREDICT
model0_ac.pred <- crosspred(model.basis, model0_ac, at=0:60/10, bylag=0.1,
  cen=volcen)
# PLOTS
# 12-DAY CUMULATIVE ASSOCIATION OF VOLUME OF TAP WATER SUPPLIED WITH SUSPECTED CHOLERA INCIDENCE IN UVIRA (Figure 3) 
tiff(filename="figure3.tif", width=4, height=3, units="in",pointsize=8, res=300)
plot(model0_ac.pred,"overall", xlab=expression(paste("Volume supplied (in 1,000 ",m^3," / day)")),ylab="RR suspected cholera incidence",ci="area", ylim=c(0,5),col="black")
optstr1<-expression("Optimal volume supply (4,790" ~ m^{3}~"/ day)")
axis(side=1, at=4.79, tick=TRUE, labels=c(optstr1), cex.axis=0.8, las=1, padj=-12, tck=0.5, lty=3, font=3)
dev.off()
# ASSOCIATION OF ABSENCE OF WATER SUPPLY WITH SUSPECTED CHOLERA INCIDENCE ON THE 12 FOLLOWING DAYS (Figure 4)
tiff(filename="figure4.tif", width=4, height=3, units="in",pointsize=8, res=300)
plot(model0_ac.pred, "slices", var=0,xlab="Lag (days)",ylab="RR cholera incidence",col="black",ci="area", ylim=c(0.8,1.3) )
dev.off()
# OVERALL CUMULATIVE RR AT VOLUME 0, WITH 95%CI
model0_ac.pred$allRRfit["0"]
model0_ac.pred$allRRlow["0"];model0_ac.pred$allRRhigh["0"]

################################################################################

# ATTRIBUTABLE FRACTION FOR LOW VOLUME - ALL TOWN

# LOAD THE FUNCTION
source('attrdl.R')

attrdl(data$vol1000, model.basis, data$cholera, model0_ac, type="af", 
  cen=volcen, range=c(0,volcen))

# EMPIRICAL CONFIDENCE INTERVALS
sim <- attrdl(data$vol1000, model.basis, data$cholera, model0_ac, type="af",
  cen=volcen, range=c(0,volcen),sim=T)
quantile(sim,c(2.5,97.5)/100)

################################################################################
# MODEL STRATIFIED BY NEIGHBOURHOOD WITH HIGH/LOW ACCESS TO TAP WATER

# MODEL CHOLERA ADMISSIONS FROM NEIGHBOURHOODS WITH LOWER TAP WATER CONSUMPTION
model0_low<- update(model0, chol_lowcov~.+chol_highcov)
dres0_low<-residuals(model0_low, "deviance")
acterm_low<-Lag(dres0_low, 1:2)
model0_ac_low<-update(model0_low, .~.+acterm_low)

# PREDICT VOLUME EFFECT IN NEIGHBOURHOODS WITH LOWER TAP WATER CONSUMPTION
model0_ac_low.pred<-crosspred(model.basis, model0_ac_low, at=0:60/10, bylag=0.1,
  cen=volcen)

# MODEL CHOLERA ADMISSIONS FROM NEIGHBOURHOODS WITH HIGHER TAP WATER CONSUMPTION
model0_high<- update(model0, chol_highcov~.+chol_lowcov)
dres0_high<-residuals(model0_high, "deviance")
acterm_high<-Lag(dres0_high, 1:2)
model0_ac_high<-update(model0_high, .~.+acterm_high)

# PREDICT VOLUME EFFECT IN NEIGHBOURHOODS WITH HIGHER TAP WATER CONSUMPTION
model0_ac_high.pred<-crosspred(model.basis, model0_ac_high, at=0:60/10, 
  bylag=0.1, cen=volcen)


# PLOTS
# 12-DAY CUMULATIVE ASSOCIATION OF VOLUME OF TAP WATER SUPPLIED WITH CHOLERA INCIDENCE STRATIFIED BY NEIGHBOURHOODS WITH LOWER AND HIGHER ACCESS TO TAP WATER IN UVIRA (FIGURE NOT SHOWN)
plot(model0_ac_low.pred,"overall",xlab=expression(paste("Volume supplied (in 1,000 ",m^3," / day)")),ylab="RR cholera incidence",ci="area", ci.arg=list(density=30, angle=-45, col=1), ylim=c(0,8),lty=2, lwd="2", col="black")
lines(model0_ac_high.pred,"overall",ci="area", ci.arg=list(density=10, angle=45, col="black"), col="black", lty=3, lwd=2)
legend("topright", seg.len=c(3,3), c("Neighbourhoods with lower tap water consumption","Neighbourhoods with higher tap water consumption"), lty=c(2,3), lwd=c(2,2), bty="n", col=c("black", "black"))
axis(side=1, at=4.79, tick=TRUE, labels=c(optstr1), cex.axis=0.8, las=1, padj=-10, tck=0.5, lty=3, font=3)

# PREDICTED INCIDENCE RATES FOR NEIGHBOURHOODS WITH HIGHER AND LOWER TAP WATER CONSUMPTION AS A FUNCTION OF THE VOLUME OF WATER SUPPLIED (FIGURE 5)
# GENERATE CROSSBASIS FOR VOLUMES TO BE PREDICTED
volpred <- matrix(0:60/10,nrow=length(0:60/10),ncol=12+1)
modelbasispred <- do.call(crossbasis,list(x=volpred,lag=12,argvar=attr(model.basis,"argvar"),arglag=attr(model.basis,"arglag")))

# GENERATE THE PREDICTION TERMS FOR COVARIATES
# Spline function of date at mid-point 
splpred <- matrix(spl[974,],nrow=nrow(modelbasispred),ncol=ncol(spl), byrow=T)
# Median rainfall
rainpred <- rep(median(data$rain,na.rm=T),nrow(modelbasispred))
# Median cholera cases admitted from the other neighbourhood
chol_highcovpred <- rep(median(data$chol_highcov,na.rm=T),nrow(modelbasispred))
chol_lowcovpred <- rep(median(data$chol_lowcov,na.rm=T),nrow(modelbasispred))
# Set day of the week
dowpred <- rep("Thursday",nrow(modelbasispred))
# Median deviance residuals at lag 1 and 2 for autocorrelation terms
medac_high<-matrix(ncol=2, nrow=1)
medac_high[1,1]<-median(acterm_high[,1], na.rm=T)
medac_high[1,2]<-median(acterm_high[,2], na.rm=T)
medac_low<-matrix(ncol=2, nrow=1)
medac_low[1,1]<-median(acterm_low[,1], na.rm=T)
medac_low[1,2]<-median(acterm_low[,2], na.rm=T)
lagreshighpred<-matrix(medac_high[1,], nrow=nrow(modelbasispred), ncol=ncol(acterm_high), byrow=T)
lagreslowpred<-matrix(medac_low[1,], nrow=nrow(modelbasispred), ncol=ncol(acterm_low), byrow=T)

# PREDICTIONS
predhigh <- predict(model0_ac_high,newdata=list(model.basis=modelbasispred,spl=splpred,rain=rainpred,dow=dowpred,chol_lowcov=chol_lowcovpred, acterm_high=lagreshighpred), type="response")*10^5/98000
predlow <- predict(model0_ac_low,newdata=list(model.basis=modelbasispred,spl=splpred,rain=rainpred,dow=dowpred,chol_highcov=chol_highcovpred, acterm_low=lagreslowpred), type="response")*10^5/106000

# PREDICTED INCIDENCE RATE OF SUSPECTED CHOLERA FOR 10'000 AS A FUNCTION OF WATER VOLUME SUPPLIED, STRATIFIED BY NEIGHBOURHOODS WITH HIGHER OR LOWER TAP WATER CONSUMPTION (FIGURE 5)
tiff(filename="figure5.tif", width=4, height=3, units="in",pointsize=8, res=300)
par(mar=c(5,5,4,1)+0.1)
plot(0:60/10,predlow,type="l",ylab="Suspected cholera incidence rate \nper 100,000 per day",xlab=expression(paste("Volume supplied (in 1,000 ",m^3," / day)")),lty=2,lwd="1",ylim=c(0,3))
lines(0:60/10,predhigh,lty=3,lwd=1)
legend("topright", c("Areas with lower tap water consumption","Areas with higher tap water consumption"), lty=c(2,3), lwd=c(1,1),cex=0.8,col=c("black", "black"))
dev.off()
par(mar=c(5,4,4,1)+0.1)
# OVERALL CUMULATIVE RR AT VOLUME 0, WITH 95%CI
# IN NEIGHBOURHOODS WITH LOWER TAP WATER CONSUMPTION
model0_ac_low.pred$allRRfit["0"]
model0_ac_low.pred$allRRlow["0"];model0_ac_low.pred$allRRhigh["0"]
# IN NEIGHBOURHOODS WITH HIGHER TAP WATER CONSUMPTION
model0_ac_high.pred$allRRfit["0"]
model0_ac_high.pred$allRRlow["0"];model0_ac_high.pred$allRRhigh["0"]


# ATTRIBUTABLE NUMBER AND FRACTION FOR LOW VOLUME IN NEIGHBOURHOODS WITH HIGHER TAP WATER CONSUMPTION
attrdl(data$vol1000, model.basis, data$chol_highcov, model0_ac_high, type="an",
  cen=volcen, range=c(0,volcen))
attrdl(data$vol1000, model.basis, data$chol_highcov, model0_ac_high, type="af",
  cen=volcen, range=c(0,volcen))

# EMPIRICAL CONFIDENCE INTERVALS
sim <- attrdl(data$vol1000, model.basis, data$chol_highcov, model0_ac_high,
  type="af",cen=volcen,range=c(0,volcen),sim=T)
quantile(sim,c(2.5,97.5)/100)

################################################################################
# SENSITIVITY ANALYSES FOR THE MODEL
#### VARIATION OF #SPLINES ####
splA <- bs(data$date2, degree=3, df=17)
modelA.basis<- crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=2, int=F))
modelA <- glm(cholera ~ modelA.basis+spl+dow+rain, data, family=quasipoisson)
dresA <- residuals(modelA, "deviance")
modelA_ac <- update(modelA, .~.+Lag(dresA, 1:2))
modelA_ac.pred <- crosspred(modelA.basis, modelA_ac, at=0:6, cen=volcen)
print(c( modelA_ac.pred$allRRfit["0"],modelA_ac.pred$allRRlow["0"],modelA_ac.pred$allRRhigh["0"]), digits=3)


splB <- bs(data$date2, degree=3, df=68)
modelB.basis<- crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=2, int=F))
modelB <- glm(cholera~modelB.basis+splB+rain+dow, data, family=quasipoisson)
dresB <- residuals(modelB, "deviance")
modelB_ac <- update(modelB, .~.+Lag(dresB, 1:2))
modelB_ac.pred <- crosspred(modelB.basis, modelB_ac, at=0:6, cen=volcen)
print(c( modelB_ac.pred$allRRfit["0"],modelB_ac.pred$allRRlow["0"],modelB_ac.pred$allRRhigh["0"]), digits=3)

##### VARIATIONS in LAG-RESPONSE FUNCTION #######
modelC.basis<-crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=2, int=T))
modelD.basis<-crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=3, int=F))
modelE.basis<-crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=3, int=T))
modelF.basis<-crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="strata", breaks=c(3, 6)))
modelG.basis<-crossbasis(data$vol1000, lag=12, argvar=list(df=1),
  arglag=list(fun="strata", breaks=c(2, 4, 6, 8)))

modelC<-glm(cholera ~ modelC.basis+spl+dow+rain, data, family=quasipoisson)
modelD<-glm(cholera ~ modelD.basis+spl+dow+rain, data, family=quasipoisson)
modelE<-glm(cholera ~ modelE.basis+spl+dow+rain, data, family=quasipoisson)
modelF<-glm(cholera ~ modelF.basis+spl+dow+rain, data, family=quasipoisson)
modelG<-glm(cholera ~ modelG.basis+spl+dow+rain, data, family=quasipoisson)

dresC <- residuals(modelC, "deviance")
dresD <- residuals(modelD, "deviance")
dresE <- residuals(modelE, "deviance")
dresF <- residuals(modelF, "deviance")
dresG <- residuals(modelG, "deviance")

modelC_ac<-update(modelC, .~.+Lag(dresC, 1:2))
modelD_ac<-update(modelD, .~.+Lag(dresD, 1:2))
modelE_ac<-update(modelE, .~.+Lag(dresD, 1:2))
modelF_ac<-update(modelF, .~.+Lag(dresF, 1:2))
modelG_ac<-update(modelG, .~.+Lag(dresG, 1:2))

modelC_ac.pred <- crosspred(modelC.basis, modelC_ac, at=0:6, cen=volcen)
modelD_ac.pred <- crosspred(modelD.basis, modelD_ac, at=0:6, cen=volcen)
modelE_ac.pred <- crosspred(modelE.basis, modelE_ac, at=0:6, cen=volcen)
modelF_ac.pred <- crosspred(modelF.basis, modelF_ac, at=0:6, cen=volcen)
modelG_ac.pred <- crosspred(modelG.basis, modelG_ac, at=0:6, cen=volcen)

print(c( modelC_ac.pred$allRRfit["0"],modelC_ac.pred$allRRlow["0"],modelC_ac.pred$allRRhigh["0"]), digits=3)
print(c( modelD_ac.pred$allRRfit["0"],modelD_ac.pred$allRRlow["0"],modelD_ac.pred$allRRhigh["0"]), digits=3)
print(c( modelE_ac.pred$allRRfit["0"],modelE_ac.pred$allRRlow["0"],modelE_ac.pred$allRRhigh["0"]), digits=3)
print(c( modelF_ac.pred$allRRfit["0"],modelF_ac.pred$allRRlow["0"],modelF_ac.pred$allRRhigh["0"]), digits=3)
print(c( modelG_ac.pred$allRRfit["0"],modelG_ac.pred$allRRlow["0"],modelG_ac.pred$allRRhigh["0"]), digits=3)


######### INCLUSION OF RESIDUAL CHLORINE LEVELS WITH DLNM STRUCTURE ####
# GENERATE 95TH QUANTILE OF RESIDUAL CHLORINE LEVELS AND CHLORINE CROSS-BASIS CENTERED AT 95TH QUANTILE
chlcen <- quantile(data$chlorine, probs=c(0.95), na.rm=T)
data$chlorine[data$vol1000==0] <- mean(data$chlorine,na.rm=T)
chlorine.basis<- crossbasis(data$chlorine, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=2, int=F))

#Model H
modelH<- glm(cholera~model.basis+chlorine.basis+spl+dow+rain, data, family=quasipoisson)
dresH<-residuals(modelH, "deviance")
modelH_ac<-update(modelH, .~.+Lag(dresH, 1:2))
modelH_ac.pred <- crosspred(model.basis, modelH_ac, at=0:6, cen=volcen)
print(c( modelH_ac.pred$allRRfit["0"],modelH_ac.pred$allRRlow["0"],modelH_ac.pred$allRRhigh["0"]), digits=3)

###### ABSENCE OF EFFECT OF LOW CHLORINE LEVEL ####
# GENERATE CHLORINE CROSS-BASIS CENTERED AT 0.8 MG/L
chlorinehigh.basis<- crossbasis(data$chlorine, lag=12, argvar=list(df=1),
  arglag=list(fun="poly", degree=2, int=F))
# RUNNING MODEL I 
modelI<- glm(cholera~model.basis+chlorinehigh.basis+spl+dow+rain, data, family=quasipoisson)
dresI<-residuals(modelI, "deviance")
modelI_ac<-update(modelI, .~.+Lag(dresI, 1:2))
# PREDICTING CHLORINE EFFECT
modelI_ac.pred <- crosspred(chlorinehigh.basis, modelI_ac, at=4:8/10, cen=0.8)
# CALCULATE RR OF CHLORINE EFFECT AT 0.4MG/L
print(c( modelI_ac.pred$allRRfit["0.4"],modelI_ac.pred$allRRlow["0.4"],modelI_ac.pred$allRRhigh["0.4"]), digits=3)


