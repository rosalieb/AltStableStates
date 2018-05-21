########################################## General informations ##########################################
### This code was written to test for regime shift hypothesis
### If you plan on using this code, please refer to our paper:
### 
### 'Seeking alternative stable states in a deep lake'
### Freshwater Biology, 2018
###
### Rosalie Bruel, Aldo Marchetto, Anaëlle Bernard, Andrea Lami, Pierre Sabatier, Victor Frossard, and Marie-Elodie Perga
### 
### DOI: 10.1111/fwb.13093
###
### contact: rosaliebruel@gmail.com AND/OR marie-elodie.perga@unil.ch
###
### If you see any error in the code, please contact us.
###
########################################## Beginning of the code ##########################################

# Add and ID for your sequence/data e.g. VAR10_10_clado
ID="JOU16_clado"
ID="example"

# Load packages
source(paste(getwd(),"/Script/packages.R",sep=""))

# Load datasets
# At the bare minimum, you should have one driver and one response to test for alternative stable state hypothesis
# Each dataframe ('response' and 'driver') should be organise with 'Years' in the first column and then 1+ columns of driver(s)/response(s)
source(paste(getwd(),"/Script/load_datasets.R",sep=""))

# 

#Check whether the data are in the correct input format
if(!exists('response') || !is.data.frame(get('response')))
  cat("\nWarning! You need to define your response(s) variable. You can prepare these data\nin the load_datasets.R script\n\n")

if(!exists('driver')|| !is.data.frame(get('driver')))
  cat("\nWarning! You need to define your driver(s) variable. You can prepare these data\nin the load_datasets.R script\n\n")

#### Change points detection ####
# Different methods: AMOC for At Most One Change, and PELT, BinSeg and SegNeigh for multiple changes
# The BinSeg is quick but approximate, PELT is exact and quick but cannot be used in all distributions, SegNeigh is exact but slow.
# Penalty: Lower penalty value results in more changepoints identified. It is up to the user to decide what value is appropriate. 
# Use an "elbow" plot to decide the penalty

# Change here the univariate data vector where you wish to detect the changepoint
y=response$Axis1
time=response$Year

# AMOC method
cptmean.AMOC <- cpt.mean(y, method="AMOC") 
(cptmean.AMOC.date<-as.numeric(time[cpts(cptmean.AMOC)]))

# PELT method with AIC penalty
cptmean.AIC <- cpt.mean(y, method="PELT",penalty = "AIC") 
(cptmean.AIC.date<-as.numeric(time[cpts(cptmean.AIC)]))

# PELT method with manual penalty (elbow method)
plot(0,0,xlim=c(0,10), ylim=c(0,15),col ="white", ylab="Number of changes in mean", xlab="Penalty", font=2,main="Elbow plot - choice of the penalty for changepoints analysis \nPruned Exact Linear Time (PELT) method")
abline(h = c(seq(from=0, to=10, by=1),15), lty=3, col=adjustcolor("black",alpha.f = 0.5))
abline(h=length(cptmean.AIC.date), lty=3, col="black", lwd=2)
for (i in seq(from = 0, to = 10, by = 0.05)) {
  cptmean=cpt.mean(y, method="PELT", penalty = "Manual", pen.value = i) 
  nb <- length(cpts(cptmean))
  points(nb~i, pch=20)
}
legend("topright",legend = c("AIC","Manual - fluctuating penalty"), title = "Methods to chose the number of changes in mean",bty="n", pch = c(NA,20), lty=c(3,NA), lwd=c(2,NA), cex = 0.8)
# change the pen.value according to the elbow plot
cptmean.Manual <- cpt.mean(y, method="PELT",penalty = "Manual", pen.value = 1) 
(cptmean.Manual.date<-as.numeric(time[cpts(cptmean.Manual)]))

# Select here the final changepoint(s) you want to keep
my_cpt <- (cptmean.Manual.date<-as.numeric(time[cpts(cptmean.Manual)]))
if(is.null(my_cpt))
  cat("\nWarning! No changepoints were detected. One of the condition for alternative stable state\nis an abrupt transition (see Andersen et al 2009, and our discussion in Bruel et al 2O18).\n\n")
if(!is.null(my_cpt)) {
  my_cpt <- as.data.frame(matrix(c(my_cpt,rep(NA,length(my_cpt))),ncol=2,byrow=F))
  for (i in 1:nrow(my_cpt)) my_cpt[i,2]=time[which(time==my_cpt[i,1])+1]
  # Re-arrange you my_cpt matrix
  my_cpt <- my_cpt[order(my_cpt[,1],decreasing = T),order(my_cpt[1,],decreasing = T)]
}

# Remove object created
rm(y,time,cptmean.AIC,cptmean.AIC.date,cptmean.AMOC,cptmean.AMOC.date,cptmean.Manual,cptmean.Manual.date)

#### Identification of dominant forcing ####
# Sometimes the package mgcv needs to be reloaded (at least on my R version)
# If you get an error message running the following code, it might be just that...
# detach("package:mgcv", unload=TRUE);library("mgcv")

# You need to create the dataframe that will be used for your analysis here. By default, there are 2 responses and 2 drivers. 
YourData <- as.data.frame(cbind(response$Year,response$Axis1,response$Axis2,driver$TP,driver$SAT))
colnames(YourData) <- c("Year","Resp1","Resp2","TP","SAT")

# Here you need to select the principal drivers
# To do that, you compare the performance of the models
modV1<-gam(Resp1~s(TP),data=YourData,method="REML")
summary(modV1)
modV2<-gam(Resp1~s(TP)+s(SAT),data=YourData,method="REML")
summary(modV2)  
modV3<-gam(Resp1~s(SAT),data=YourData,method="REML")
summary(modV3) 

# Lowest AIC indicates best model
AIC(modV1,modV2,modV3)
# Run log-likelihood ratio tests too: chose the max likelihood, and the lowest p-value
lrtest(modV1,modV2,modV3)

# 'Save' here which model you selected
selected.gam <- modV2

#### Autocorrelation test GAM ####
gam.ac<-gamm(Resp1~s(TP)+s(SAT),data=YourData,family=gaussian(link="identity"),method="ML")
summary(gam.ac$gam)
par(mfrow=c(1,2))
acf(residuals(gam.ac$lme,type="normalized"))
pacf(residuals(gam.ac$lme,type="normalized"))
acf(residuals(gam.ac.axis2$lme,type="normalized"))
pacf(residuals(gam.ac.axis2$lme,type="normalized"))

# Random tests sur les residus 
gam.res<-residuals(selected.gam,type="deviance")
plot(YourData$Year,gam.res,xlab="dates",ylab="residuals TP",main="DCA1")
gam.res.smooth <- smooth.spline(YourData$Year,gam.res,spar=0.8)
lines(predict(gam.res.smooth),col="blue",lwd=2)
bartels.rank.test(gam.res)        #meme resultats entre les 2 tests
difference.sign.test(gam.res)

plot(YourData$Resp1,gam.res,xlab="observed",ylab="residuals driver",main="")
gam.res.smooth <- smooth.spline(YourData$Resp1,gam.res,spar=0.8)
lines(predict(gam.res.smooth),col="blue",lwd=2)
gam.fit<-lm(gam.res~YourData$Resp1)
summary(gam.fit)  

# If AR(1) and non random, use GAM with Correlation structure AR(1) or CAR(1). CAR(1) because uneven sampling in time
# The need to include the correlation structure can likewise be determined using a likelihood
# ratio test (LRT) or via AIC (or BIC) by comparing a model fitted using the sturcture with a
# model with independance of observations. The later model is a simpler model with fewer
# parameters, and in the case of AR(1) and CAR(1), we are testing whether phi-estimate is
# significantly different from 0. REML should be used to fit the models when performing these
# tests, and furthermore, as 0 is on the boundary of allowed values for phi, the p-value may be
# anticonservative (source: Simpson and Anderson, 2009, Limnology & Oceanography).
# Note: penalized quasilikelihood (PQL) and restricted maximum likelihood (REML)
gam.no.c<-gam(Resp1~s(TP)+s(SAT),data=YourData,method="REML")
gam.no.c.te<-gam(Resp1~te(TP)+te(SAT),data=YourData,method="REML")
gam.c <- gamm(Resp1~s(TP)+s(SAT),data=YourData,family=gaussian(link="identity"),method="REML") 
lrtest(gam.no.c,gam.no.c.te)
lrtest(gam.c$gam)
lrtest(gam.no.c)
gam.c.lme <- gam.c$lme
gam.c.lme$logLik
# More info on log-likelihood ratio test:
# http://stats.stackexchange.com/questions/6505/likelihood-ratio-test-in-r


#### OPTIONAL - GAM selection - looking at your data ####
# Look at your data: plot the model you selected
# You can skip this step and go directly to computing the synthesis plot
plot.gam(gam.no.c)
# legend:  upper and lower lines are added to the 1-d plots at 2 standard errors above and below the estimate of the smooth
# sewithmean: if TRUE the component smooths are shown with confidence intervals that include the uncertainty about the overall mean. If FALSE then the uncertainty relates purely to the centred smooth itself. Marra and Wood (2012) suggests that TRUE results in better coverage performance, and this is also suggested by simulation.
# Note that, if seWithMean=TRUE, the confidence bands include the uncertainty about the overall mean. In other words although each smooth is shown centred, the confidence bands are obtained as if every other term in the model was constrained to have average 0, (average taken over the covariate values), except for the smooth concerned. This seems to correspond more closely to how most users interpret componentwise intervals in practice, and also results in intervals with close to nominal (frequentist) coverage probabilities by an extension of Nychka's (1988) results presented in Marra and Wood (2012).
# Here the plot will be saved in your file
description = "PC1_example" # Add here something to better describe your data
if(is.null(description)) description=Sys.Date()
# Change the labels in the code (xlab and ylab terms)
pdf(paste(getwd(),"/Output/Figures/",ID,"_",description,"_GAM.pdf",sep=""), width=10, height=5, family="Helvetica")
par(mfrow=c(1,2))
plot(gam.no.c, residuals=F, rug=T, se=T, scale=-1, select = 1, xlab="TP", ylab="Smooth function (TP)", shade=T, seWithMean=T, lwd=2)
abline(h = 0, col="black", lty = "dotted")
plot(gam.no.c, residuals=F, rug=T, se=T, scale=-1, select = 2, xlab="SAT anomalies (°C)", ylab="Smooth function (SAT)", shade=T,seWithMean=T, lwd=2)
abline(h = 0, col="black", lty = "dotted")
dev.off()

# Contributions plots
env.gam<-YourData[,c(1,4,5)] #Select only the driver variable
predictgam<-predict.gam(selected.gam,env.gam,type="terms",se.fit=TRUE)
predictgam.selected.model<-predictgam
low<-predictgam$fit-1.96*predictgam$se.fit
high<-predictgam$fit+1.96*predictgam$se.fit
# s(TP) - It was our first driver, you can change it for whichever driver you found significant
plot(YourData$Year,predictgam$fit[,1],type="l",ylim=c(min(low[,1]), max(high[,1])),xlab="Year",ylab="TP",col="black",lwd=2,family="Helvetica")
polygon.x <- c(YourData$Year,rev(YourData$Year))
polygon.y <- c(low[,1], rev(high[,1]))
polygon(x=polygon.x, y=polygon.y, col=adjustcolor("black", alpha.f=0.4), border=NA)
abline(h = 0, col="black", lty = "dotted")
for (i in 1:nrow(my_cpt)) rect(my_cpt[i,1], -10, my_cpt[i,2], 10, col = adjustcolor("black", alpha.f=0.2), border = "transparent", lty = par("lty"), lwd = par("lwd"))

# s(SAT) - It was our second driver, you can change it for whichever driver you found significant
plot(YourData$Year,predictgam$fit[,2],type="l",ylim=c(min(low[,2]), max(high[,2])),xlab="Year",ylab="SAT",col="black",lwd=2,family="Helvetica")
polygon.x <- c(YourData$Year,rev(YourData$Year))
polygon.y <- c(low[,2], rev(high[,2]))
polygon(x=polygon.x, y=polygon.y, col=adjustcolor("black", alpha.f=0.4), border=NA)
abline(h = 0, col="black", lty = "dotted")
for (i in 1:nrow(my_cpt)) rect(my_cpt[i,1], -10, my_cpt[i,2], 10, col = adjustcolor("black", alpha.f=0.2), border = "transparent", lty = par("lty"), lwd = par("lwd"))

##### Synthesis Plot ####
# This code allows you to reproduce Figures 4 and 5 from our paper Bruel et al 2018, Freshwater Biology
# Prepare your data
y = rev(YourData$Resp1)
time = rev(YourData$Year)
Driver1 = rev(YourData[,4])
Driver2 = rev(YourData[,5])

# Prepare the gam
selected.gam <- gam.no.c

env.gam<-as.data.frame(cbind(time,Driver1,Driver2)) # Select only the driver variable
colnames(env.gam) <- c("Year", "TP", "SAT")# Important to put the same name as the "selected.gam"
predictgam<-predict.gam(selected.gam,env.gam,type="terms",se.fit=TRUE)
predictgam.selected.model<-predictgam
low<-predictgam$fit-1.96*predictgam$se.fit
high<-predictgam$fit+1.96*predictgam$se.fit

# Generate early warning signals, Dakos et al. 2008, 2012
EWS_Response <- generic_ews(ts(y,time),winsize=20,detrending="loess",bandwidth=15, interpolate=TRUE,AR_n=FALSE) 
# generating an interpolated time series (loess priviledged based on comments to the Wang et al, 2012 paper)
yint <- loess(y~time)
y4<-predict(yint, time)
res.y<-y-y4
# residuals with standardization (res value multiplied by sqrt of the time period since the more averaging, the least variation you expect)
interval<-diff(time)
res.ystd<-(y[2:length(time)]-y4[2:length(time)])*sqrt(interval)
EWS_Response_res <- generic_ews(ts(res.ystd,time),winsize=20,detrending="loess",bandwidth=15, interpolate=TRUE,AR_n=FALSE) 

description = "PC1" # Add here something to better describe your data
if(is.null(description)) description=Sys.Date()

# Generate the output plot
pdf(paste(getwd(),"/Output/Figures/",ID,"_",description,"_GAM_plot.pdf",sep=""), width=10, height=10.5,family="Helvetica")
color_export = T #Change to T/F if plot color or black/white for MAAT
color_export2 = F #Change to T/F if plot color or black/white for TP
color_export3 = F #Change to T/F if plot color or black for legend axis 4 on panel c
driver_state_phase = "Driver1" #Change to "Driver1"/"Driver2" for state phase plot.
# Driver 1 will be the second column in the env.gam dataframe generated earlier.
# Driver 2 will be the third colum in the env.gam dataframe generated earlier.
#vec.pch <- c(rep(18,7),rep(16,16),rep(17,9),rep(2,42)) # with full pch for recent
vec.pch <- 20
cex_1 = .9 # cex text axis+legend
cex_2 = 1.2 # cex writing inside panel i.e. panel b
cex_3 = 1.2 # cex.axis
cex_4 = 1.2 # cex legends
layout(matrix(c(1,0,2,2,3,0,4,5,6,0,7,7), 3, 4, byrow = TRUE),width=c(1.2,0.05,.95,0.35))
par(mar=c(4.1,4.1,4.1,4.1)) #par(mar=c(5.1,4.1,4.1,2.1))

# Plot your response data with the apparent breaks
plot(time, y, xlab="", ylab="", type='n',xlim=c(min(time), max(time)),axes=F)
for (i in 1:nrow(my_cpt)) rect(my_cpt[i,1], -10, my_cpt[i,2], 10, col = adjustcolor("black", alpha.f=0.2), border = "transparent", lty = par("lty"), lwd = par("lwd"))
lines(time, y, lty=2, col="grey")
par(new=T)
plot(time, y, xlab="", ylab="", pch=vec.pch,axes=F, cex=1.3)
axis(1, at=c(min(time)*.8, max(time)*1.2))
axis(1, cex.axis=cex_3)
mtext("Year",side = 1, line=2.5,cex=cex_1)
axis(2, at=c(min(y)-abs(max(y)), max(y)+abs(max(y))))
axis(2, cex.axis=cex_3)
mtext("Response",side = 2, line=2.5,cex=cex_1)
mtext("a", 3, adj=0, line=2,col="black",font=2, family="Helvetica")

# F density plot
if(is.null(my_cpt)) {
  d11 <- density(y, bw = 0.1)
  plot(d11, xlim=c(min(y),max(y)), main="", ylab="", xlab="", lwd=2, axes=F)
} else {
  # Re-arrange your my_cpt matrix
  my_cpt <- my_cpt[order(my_cpt[,1],decreasing = T),order(my_cpt[1,],decreasing = T)]
  
  # Determine limits plots (this is a badly written code to avoid you to set manually the limits)
  # There is necessary a better way to do it
  myx <- NULL
  myy <- NULL
  d11 <- density(y[time>=max(my_cpt[1,])], bw = 0.1)
  myx <- c(myx, min(d11$x),max(d11$x))
  myy <- c(myx, min(d11$y),max(d11$y))
  for (i in 1:nrow(my_cpt)) {
    if(i<nrow(my_cpt)) {
      d11 <- density(y[time>=max(my_cpt[i+1,])&time<=min(my_cpt[i,])], bw = 0.1)
      myx <- c(myx, min(d11$x),max(d11$x))
      myy <- c(myx, min(d11$y),max(d11$y))
    } else {
      d11 <- density(y[time<=min(my_cpt[1,])], bw = 0.1) 
      myx <- c(myx, min(d11$x),max(d11$x))
      myy <- c(myx, min(d11$y),max(d11$y))
    }
  }
  
  d11 <- density(y[time>=max(my_cpt[1,])], bw = 0.1)
  # You may want to change the y scale for your plot here, in the ylim term.
  plot(d11, xlim=c(min(myy)*.8,max(myy)*1.3), ylim=c(min(myx),max(myx)*1.2), main="", ylab="", xlab="", lwd=nrow(my_cpt)+1, axes=F)
  par(xpd=T)
  text(median(d11$x),max(d11$y),paste("> ", max(my_cpt[1,])), pos = 3,cex=cex_2)
  par(xpd=F)
  
  for (i in 1:nrow(my_cpt)) {
    if(i<nrow(my_cpt)) {
      d12 <- density(y[time>=max(my_cpt[i+1,])&time<=min(my_cpt[i,])], bw = 0.1)
      lines(d12,lwd=nrow(my_cpt)+1-i)
      par(xpd=T)
      text(median(d12$x),max(d12$y),paste(min(my_cpt[i,]),max(my_cpt[i+1,]),sep="-"), pos = 3,cex=cex_2) 
      par(xpd=F)
    } else {
      d12 <- density(y[time<=min(my_cpt[1,])], bw = 0.1) 
      lines(d12,lwd=nrow(my_cpt)+1-i)
      par(xpd=T)
      text(median(d12$x),max(d12$y),paste("< ",min(my_cpt[i,]),sep=""), pos = 3,cex=cex_2)
      par(xpd=F)
    }
  }
  rm(myx,myy)
}
axis(1, at=c(min(y)-abs(max(y)), max(y)+abs(max(y))))
axis(1, cex.axis=cex_3)
mtext("Distribution response scores",side = 1, line=2.5,cex=cex_1)
axis(2, at=c(-5,20), cex.axis=cex_3)
axis(2, cex.axis=cex_3)
mtext("Density",side = 2, line=2.5,cex=cex_1)
mtext("b", 3, adj=0, line=2,col="black",font=2, family="Helvetica")

# Contributions plots
# Plot the background (breaks)
plot(time,predictgam$fit[,1],type="n",xlab="",ylab="",xlim=c(min(time), max(time)), ylim=c(min(low), max(high)),axes=F)
abline(h = 0, col="black", lty = "dotted")
for (i in 1:nrow(my_cpt)) rect(my_cpt[i,1], -10, my_cpt[i,2], 10, col = adjustcolor("black", alpha.f=0.2), border = "transparent", lty = par("lty"), lwd = par("lwd"))

# Second driver (plotted in the background)
par(new=T)
if (color_export) {plot(time,predictgam$fit[,2],type="l",xlim=c(min(time), max(time)),ylim=c(min(low), max(high)),xlab="",ylab="",col="darkorange",lwd=2,family="Helvetica",axes=F)} else {plot(time,predictgam$fit[,2],type="l",xlab="",xlim=c(min(time), max(time)),ylim=c(min(low), max(high)),ylab="",col="black",lwd=1,family="Helvetica",axes=F)}
if (color_export3) {
  axis(4,at=c(min(low)*2,max(high)*2),cex.axis=cex_3, col="darkorange")
  axis(4,cex.axis=cex_3, col="darkorange")
  mtext("s(Driver 2)",side = 4, line=2.5,cex=cex_1, col="darkorange")
} else {
  axis(4,at=c(min(low)*2,max(high)*2),cex.axis=cex_3)
  axis(4,cex.axis=cex_3)
  mtext("s(Driver 2)",side = 4, line=2.5,cex=cex_1)
  legend("topleft",col = c("black", "darkorange"), lwd=2, legend = c("s(First driver)", "s(Second driver)"), cex=cex_4, bty = "n")
}
polygon.x <- c(time,rev(time))
polygon.y <- c(low[,2], rev(high[,2]))
if (color_export) {polygon(x=polygon.x, y=polygon.y, col=adjustcolor("darkorange", alpha.f=0.3), border=NA)} else {lines(time,low[,2],lty=3);lines(YourData$Year,high[,2],lty=3)}
# First driver 
par(new=T)
if (color_export2) {plot(time,predictgam$fit[,1],type="l",xlim=c(min(time), max(time)),ylim=c(min(low), max(high)),xlab="",ylab="",col="forestgreen",lwd=2,family="Helvetica",axes=F)} else {plot(time,predictgam$fit[,1],type="l",xlim=c(min(time), max(time)),ylim=c(min(low), max(high)),xlab="",ylab="",col="black",lwd=2,family="Helvetica",axes=F)}
axis(1, at=c(min(time)*.8, max(time)*1.2))
axis(1, cex.axis=cex_3)
axis(2,at=c(min(low)*2,max(high)*2),cex.axis=cex_3)
axis(2,cex.axis=cex_3)
mtext("Year",side = 1, line=2.5,cex=cex_1)
mtext("s(First driver)",side = 2, line=2.5,cex=cex_1)
polygon.x <- c(time,rev(time))
polygon.y <- c(low[,1], rev(high[,1]))
if (color_export2) {polygon(x=polygon.x, y=polygon.y, col=adjustcolor("forestgreen", alpha.f=0.2), border=NA)
} else {polygon(x=polygon.x, y=polygon.y, col=adjustcolor("black", alpha.f=0.2), border=NA)}

mtext("c", 3, adj=0, line=2,col="black",font=2, family="Helvetica")

# Response vs. forcing
cex_sITP <- abs(predictgam$fit[,1])*2.3
cool = rainbow(10, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(10, start=rgb2hsv(col2rgb('firebrick'))[1], end=rgb2hsv(col2rgb('orange'))[1])
mypalette = c(rev(cool), rev(warm))  #; mypalette <- colorRampPalette(cols)(255)
mypalette <- colorRampPalette(rev(warm))(19)
col_sDriver2 <- round(((predictgam$fit[,2]+4)), digits=0)
col_sDriver2 <- mypalette[col_sDriver2]
par(mar=c(4.1,4.1,4.1,0.7))
if (driver_state_phase == "Driver2") {
  plot(env.gam[,3], YourData$Resp1, pch=vec.pch,col=col_sDriver2,xlab="",ylab="", cex= cex_sITP, axes=F)
  lines(env.gam[,3],YourData$Resp1, lty=2, col=adjustcolor("black", alpha.f=0.5))
  axis(1, at=c(min(env.gam[,3])-min(env.gam[,3]), max(env.gam[,3])+max(env.gam[,3])), cex.axis=cex_3)
  axis(1, cex.axis=cex_3)
  mtext("Driver 2",side = 1, line=2.5,cex=cex_1)
} else {
  plot(env.gam[,2], YourData$Resp1, pch=vec.pch,col=col_sDriver2,xlab="",ylab="", cex= cex_sITP, axes=F)
  lines(env.gam[,2], YourData$Resp1, lty=2, col=adjustcolor("black", alpha.f=0.5))
  axis(1, at=c(min(env.gam[,2])-abs(max(env.gam[,2])), max(env.gam[,2])+max(env.gam[,2])), cex.axis=cex_3)
  axis(1, cex.axis=cex_3)
  mtext("Driver 1",side = 1, line=2.5,cex=cex_1)
}
axis(2, at=c(min(y)-abs(max(y)), max(y)+max(y)), cex.axis=cex_3)
axis(2, cex.axis=cex_3)
mtext("Response",side = 2, line=2.5,cex=cex_1)

mtext("d", 3, adj=0, line=2,col="black",font=2, family="Helvetica")

legend_image <- as.raster(matrix(rev(mypalette), ncol=1))
par(mar=c(2.1,0.7,2.1,0.7))
plot(c(-0.3,2.3),c(0.1,1.1),type = 'n', axes = F,xlab = '', ylab = '')
points(c(0.5,0.5,0.5),c(0.65,0.75,0.85), pch=16,cex=c(0.8,1.5,2.2))
text(x=1.5, y = c(0.65,0.75,0.85), labels = c("Low","Moderate","High"),cex=cex_4)
text(0.6,0.9, labels=c("s(Driver 1)"), pos = 3, cex=cex_4+0.1, font=2)
rasterImage(legend_image, 0, 0.2, 1,0.4)
text(x=1.5, y = c(0.2,0.4), labels = c("Low","High"),cex=cex_4)
text(0.6,0.48, labels=c("s(Driver 2)"), cex=cex_4+0.1, font=2)
par(mar=c(4.1,4.1,4.1,4.1))

#EWS
#changes in Autocorrelation 
plot(time[(length(time)-length(EWS_Response$acf1)+1):(length(time))],EWS_Response$acf1,col="black",type="l",main="",xlab="",ylab="",xlim=c(min(time), max(time)),lwd=2, axes=F)
for (i in 1:nrow(my_cpt)) rect(my_cpt[i,1], -10, my_cpt[i,2], 10, col = adjustcolor("black", alpha.f=0.2), border = "transparent", lty = par("lty"), lwd = par("lwd"))
axis(1, at=c(min(time)*.8, max(time)*1.2))
axis(1, cex.axis=cex_3)
soustraction = max(abs(max(EWS_Response$acf1)),abs(min(EWS_Response$acf1)))
axis(2, at=c(min(EWS_Response$acf1)-soustraction, max(EWS_Response$acf1)+soustraction))
axis(2, cex.axis=cex_3)
mtext("Year",side = 1, line=2.5,cex=cex_1)
mtext("AR(1) of raw response scores",side = 2, line=2.5,cex=cex_1)
#residuals standardized by resolution changes in autocorrelation 
par(new=T)
plot(time[(length(time)-length(EWS_Response_res$acf1)+1):(length(time))],EWS_Response_res$acf1,col="black",type="l",main="",xlab="",ylab="",xlim=c(min(time), max(time)),lwd=1, axes=F)
soustraction = max(abs(max(EWS_Response_res$acf1)),abs(min(EWS_Response_res$acf1)))
axis(4, at=c(min(EWS_Response_res$acf1)-soustraction, max(EWS_Response_res$acf1)+soustraction))
axis(4, cex.axis=cex_3)
mtext("AR(1) of raw residuals response scores",side = 4, line=2.5,cex=cex_1)
legend("topleft",legend=c("Raw", "Residuals"),lwd=c(2, 1), col=c("black", "black"), cex=cex_4, bty = "n")
mtext("e", 3, adj=0, line=2,col="black",font=2, family="Helvetica")

#changes in SD 
plot(time[(length(time)-length(EWS_Response$sd)+1):(length(time))],EWS_Response$sd,col="black",type="l",main="",xlab="",ylab="",xlim=c(min(time), max(time)),lwd=2, axes=F)
for (i in 1:nrow(my_cpt)) rect(my_cpt[i,1], -10, my_cpt[i,2], 10, col = adjustcolor("black", alpha.f=0.2), border = "transparent", lty = par("lty"), lwd = par("lwd"))
axis(1, at=c(min(time)*.8, max(time)*1.2))
axis(1, cex.axis=cex_3)
soustraction = max(abs(max(EWS_Response$sd)),abs(min(EWS_Response$sd)))
axis(2, at=c(min(EWS_Response$sd)-soustraction, max(EWS_Response$sd)+soustraction))
axis(2, cex.axis=cex_3)
mtext("Year",side = 1, line=2.5,cex=cex_1)
mtext("SD of raw response scores",side = 2, line=2.5,cex=cex_1)
#residuals standardized by resolution changes in SD 
par(new=T)
plot(time[(length(time)-length(EWS_Response_res$sd)+1):(length(time))],EWS_Response_res$sd,col="black",type="l",main="",xlab="",ylab="", xlim=c(min(time), max(time)),lwd=1, axes=F)
soustraction = max(abs(max(EWS_Response_res$sd)),abs(min(EWS_Response_res$sd)))
axis(4, at=c(min(EWS_Response_res$sd)-soustraction, max(EWS_Response_res$sd)+soustraction))
axis(side = 4, cex.axis=cex_3)
mtext(side = 4, 'SD of residuals response scores', line=2.5,cex=cex_1)
legend("topleft",legend=c("Raw", "Residuals"),lwd=c(2, 1), col=c("black", "black"), cex=cex_4, bty="n")
mtext("f", 3, adj=0, line=2,col="black",font=2, family="Helvetica")

dev.off()

#### End of the code ####
