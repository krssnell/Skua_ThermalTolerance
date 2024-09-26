library(MASS)  # used to fit negative binomial distribution
library(pscl)  # used to fit zero-inflated and hurdle models
library(lmtest) # used to perform likelihood ratio tests
library(multcomp) # used to perform multiple comparisons
library(RColorBrewer) # for plots

##############
# USE DATA df from 0 DataWrangling.R
load("df.RData")
OFFdf <- df[,c(6:8,12)] #Temp, wind, sun exposure, time off the nest per obs block


###########################################
##  MODEL COMPARISON AND GOODNESS-OF-FIT ##
###########################################

# FIT COUNT MODELS WITH DIFFERENT DISTRIBUTIONAL ASSUMPTIONS
# Below a Poisson, negative binomial, zero-inflated Poisson, zero-inflated negative binomial,
# Poisson logit hurdle and negative binomial logit hurdle model are fitted, respectively.
# By using the . for the predictor function, all available predictors 
# (TEMP, WIND, SEXP) are included as main effects.

model1p<- glm(OFF ~ .,data= OFFdf,family=poisson)
model1nb<- glm.nb(OFF  ~ . , data=OFFdf)
model1zip<-zeroinfl(OFF~.|.,data = OFFdf, link = "logit", dist = "poisson")
model1zinb<-zeroinfl(OFF~.|.,data = OFFdf, link = "logit", dist = "negbin")
model1plh<-hurdle(OFF~.|.,data=OFFdf,zero.dist="binomial",link="logit",dist="poisson")
model1nblh<-hurdle(OFF~.|.,data=OFFdf,zero.dist="binomial",link="logit",dist="negbin", family="binomial")
summary(model1nblh)
# CALCULATE FOR EACH MODEL THE PREDICTED FREQUENCIES FOR EACH OUTCOME
estobsk<-function(k){
cbind(sum(dpois(k,fitted(model1p))),sum(dnbinom(k,mu=fitted(model1nb),size=model1nb$theta)),sum(predict(model1zip,type="prob")[,1+k]),sum(predict(model1zinb,type="prob")[,1+k]),
sum(predict(model1plh,type="prob")[,1+k]),sum(predict(model1nblh,type="prob")[,1+k]))}

table(OFFdf$OFF)
# predicted frequencies are calculated for outcomes from 0 to 15,

estobs<-sapply(0:15,estobsk) # max time off nest was 9 mins 

# CREATE A HISTOGRAM CONTRASTING OBSERVED WITH PREDICTED FREQUENCIES
# This is done for each of the fitted models above,
# The breaks in the hist-function were explicitly specified as typically this function bins zero and one together.
# Alternatively one can plot the frequencies using the following command: plot(table(OFFdf$OFF))

hist(OFFdf$OFF,xlab="Number mins off the nest observed",ylab="Absolute Frequency",breaks=seq(-0.5,15.5,by=1),main="",cex=1.2)
points(0:15-0.1,estobs[1,],col=1,lwd=2,pch=2,lty=2)
lines(0:15-0.1,estobs[1,],col=1,lwd=2,pch=2,lty=2)
points(0:15-0.05,estobs[2,],col=1,lwd=2,pch=3,lty=3)
lines(0:15-0.05,estobs[2,],col=1,lwd=2,pch=3,lty=3)
points(0:15+0.05,estobs[3,],col=1,lwd=2,pch=1,lty=2)
lines(0:15+0.05,estobs[3,],col=1,lwd=2,pch=1,lty=2)
points(0:15+0.1,estobs[4,],col=1,lwd=2,pch=4,lty=3)
lines(0:15+0.1,estobs[4,],col=1,lwd=2,pch=4,lty=3)
points(0:15+0.05,estobs[5,],col=3,lwd=2,pch=1,lty=2)
lines(0:15+0.05,estobs[5,],col=3,lwd=2,pch=1,lty=2)
points(0:15+0.1,estobs[6,],col=3,lwd=2,pch=4,lty=3)
lines(0:15+0.1,estobs[6,],col=3,lwd=2,pch=4,lty=3)
# mtext(">10",side=1,at=11)
# legend(4,150,c("Poisson","Negative Binomial","Zero-Inflated Poisson","Zero-Inflated Negative Binomial"),lty=c(2,3,2,3),pch=c(2,3,1,4),lwd=2,bty="n",cex=1.2)
legend(4,400,c("Poisson","Negative Binomial","Zero-Inflated Poisson","Zero-Inflated Negative Binomial", "Poisson logit hurdle", "negative binomial logit hurdle model"),lty=c(2,3,2,3,2,3),pch=c(2,3,1,4,1,4), col=c(1,1,1,1,3,3),lwd=2,bty="n",cex=1.2)



# PERFORM MODEL COMPARISONS USING LRT, VUONG'S TEST AND AIC 

lrtest(model1p,model1nb)
vuong(model1nb,model1zinb)
t(AIC(model1p,model1nb,model1zip,model1zinb,model1plh,model1nblh))

###############################################
##  MODELLING AND INTERPRETING MAIN EFFECTS  ## 
###############################################

## FIT MAIN EFFECT MODELS FOR TEMP  ##
# As from the previous step the zero-inflated negative binomial yielded the best fit,
# all further models use these underlying models.
# Note that the notation y~x1|z1 implies the count component model y~x1 conditional on (|) the hurdle
# or zero-inflated model y~z1 for the zero-component, in our example we always use x1=z1

# NULL MODEL (WITH NO PREDICTORS)
model0zinb<-zeroinfl(OFF~1|1,data = OFFdf, link = "logit", dist = "negbin")
model0nblh<-hurdle(OFF~1|1,data = OFFdf,zero.dist="binomial",link="logit",dist="negbin")

# MODEL WITH MAIN EFFECT OF TEMP 
zinb_TEMP<-zeroinfl(OFF~TEMP|TEMP,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_TEMP)
nblh_TEMP<-hurdle(OFF~TEMP|TEMP,data = OFFdf,zero.dist="binomial",link="logit",dist="negbin")
summary(nblh_TEMP)

# MODEL WITH MAIN EFFECT OF WIND
# similar approach can be used to explore to effect of WIND speed alone
zinb_WIND<-zeroinfl(OFF~WIND|WIND,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_WIND)
nblh_WIND<-hurdle(OFF~WIND|WIND,data = OFFdf,zero.dist="binomial",link="logit",dist="negbin")
summary(nblh_WIND)

# MODEL WITH MAIN EFFECT OF SUN EXPOSURE
# similar approach can be used to explore to effect of WIND speed alone
zinb_SEXP<-zeroinfl(OFF~SEXP|SEXP,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_SEXP)
nblh_SEXP<-hurdle(OFF~SEXP|SEXP,data = OFFdf,zero.dist="binomial",link="logit",dist="negbin")
summary(nblh_SEXP)

############################
## MULTIPLE TESTING ISSUE ##
############################

# Does TEMP have a significant impact on the time OFF the nest?

# Use LRT to jointly test that both the impact on the zero-component and count-component are zero
lrtest(zinb_TEMP,model0zinb)
lrtest(nblh_TEMP,model0nblh)


# Similarly, one can assess the effect of WIND on the time OFF the nest

lrtest(zinb_WIND,model0zinb)
lrtest(nblh_WIND,model0nblh)

# Similarly, one can assess the effect of sun exposure on the PANT-outcome
lrtest(zinb_SEXP,model0zinb)
lrtest(nblh_SEXP,model0nblh)

##############################
##############################
tiff(file = "SM Fig4a.tiff", units="cm",width=19,height=6,res=500)
par(mfrow = c(1, 3), xpd=NA,cex.axis=0.8, cex.lab=1, mgp=c(1.1,.25,0), tcl=-0.2, mar=c(2,2,2,1), ps=10)


plot(OFFdf$TEMP,OFFdf$OFF, xlab="Air temperature (\u00b0C)", ylab="Nest absence", pch=20, cex=0.6, col=alpha("black", 0.2))
mtext(text="A", side = 3, adj = -0.05, line = 0, cex=0.8)

plot(OFFdf$WIND,OFFdf$OFF, xlab=expression(Wind~speed~(ms^-1)), ylab="", pch=20, cex=0.6, col=alpha("black", 0.2))
mtext(text="B", side = 3, adj = -0.05, line = 0, cex=0.8)


plot(jitter(OFFdf$SEXP),jitter(OFFdf$OFF), xlab="Sun exposure index", ylab="", pch=20, cex=0.6, col=alpha("black", 0.2))
mtext(text="C", side = 3, adj = -0.05, line = 0, cex=0.8)



dev.off()

##############################
##############################

################################################################
## PRESENTING AND INTERPRETING INTERACTIONS IN MIXTURE MODELS ##
################################################################

## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN TEMP AND WIND
zinb_TEMPWIND2<-zeroinfl(OFF~TEMP*WIND|TEMP*WIND,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_TEMPWIND2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_TEMPWIND<-zeroinfl(OFF~TEMP+WIND|TEMP+WIND,data = OFFdf, link = "logit", dist = "negbin")
lrtest(zinb_TEMPWIND2,zinb_TEMPWIND)



## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN TEMP AND SUN EXPOSURE
zinb_TEMPSEXP2<-zeroinfl(OFF~TEMP*SEXP|TEMP*SEXP,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_TEMPSEXP2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_TEMPSEXP<-zeroinfl(OFF~TEMP+SEXP|TEMP+SEXP,data = OFFdf, link = "logit", dist = "negbin")
lrtest(zinb_TEMPSEXP2,zinb_TEMPSEXP)



## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN WIND AND SUN EXPOSURE
zinb_WINDSEXP2<-zeroinfl(OFF~WIND*SEXP|WIND*SEXP,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_WINDSEXP2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_WINDSEXP<-zeroinfl(OFF~WIND+SEXP|WIND+SEXP,data = OFFdf, link = "logit", dist = "negbin")
lrtest(zinb_WINDSEXP2,zinb_WINDSEXP)

#Full model includnig TEMP, WIND, and SEXP
zinb_FULL <- zeroinfl(OFF~.|.,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_FULL)
# Backward stepwise variable selection with significance level alpha=0.01.
library("mpath")
fitbe <- be.zeroinfl(zinb_FULL, data=OFFdf, dist="negbin", alpha=0.01, trace=FALSE)
summary(fitbe)

# graphical presentation of effect on both components of the zero-inflated model and the mean 
newdata<-expand.grid(TEMP=seq(from=min(OFFdf$TEMP),to=max(OFFdf$TEMP), length.out=100),WIND=seq(from=min(OFFdf$WIND),to=max(OFFdf$WIND), length.out=100))

yhatmean<-predict(zinb_TEMPWIND2,newdata=newdata,type="response")
yhatcount<-predict(zinb_TEMPWIND2,newdata=newdata,type="count")
yhatzero<-predict(zinb_TEMPWIND2,newdata=newdata,type="zero")

par(mfrow=c(1,3))
plot(newdata$TEMP,yhatcount,xlab="Air Temp",ylab=expression(paste(mu[i])),pch=20,col="lightgrey",cex=1.2)
plot(newdata$TEMP,1-yhatzero,xlab="Air Temp",ylab=expression(paste(1-p[i])),pch=20,col="lightgrey",cex=1.2)
plot(newdata$TEMP,yhatmean,xlab="Air Temp",ylab="Predicted Overall Mean",pch=20,col="lightgrey",cex=1.2)


# graphical presentation of effect on TEMP component and WIND (mean, +1SD, -1SD) of the zero-inflated model and the mean
WINDmean <- mean(OFFdf$WIND)
WINDminSD <- WINDmean - sd(OFFdf$WIND)
WINDplusSD <- WINDmean + sd(OFFdf$WIND)
newdata<-expand.grid(TEMP=seq(from=min(OFFdf$TEMP),to=max(OFFdf$TEMP), length.out=100),WIND=c(WINDminSD,WINDmean,WINDplusSD))

yhatmean<-predict(zinb_TEMPWIND2,newdata=newdata,type="response")
yhatcount<-predict(zinb_TEMPWIND2,newdata=newdata,type="count")
yhatzero<-predict(zinb_TEMPWIND2,newdata=newdata,type="zero")

cols <- brewer.pal(n=3, name="Blues")

par(mfrow=c(1,3))
plot(newdata$TEMP,yhatcount,xlab="Air Temp",ylab=expression(paste(mu[i])),pch="",col="lightgrey",cex=1.2)
lines(newdata$TEMP[c(101:200)],yhatcount[c(101:200)],lty=1, col=cols[2])#meanwind
lines(newdata$TEMP[c(1:100)],yhatcount[c(1:100)],lty=1, col=cols[1])#low wind
lines(newdata$TEMP[c(201:300)],yhatcount[c(201:300)],lty=1, col=cols[3])#high wind
plot(newdata$TEMP,1-yhatzero,xlab="Air Temp",ylab=expression(paste(1-p[i])),pch="",col="lightgrey",cex=1.2)
lines(newdata$TEMP[c(101:200)],1-yhatzero[c(101:200)],lty=1, col=cols[2])#meanwind
lines(newdata$TEMP[c(1:100)],1-yhatzero[c(1:100)],lty=1, col=cols[1])#low wind
lines(newdata$TEMP[c(201:300)],1-yhatzero[c(201:300)],lty=1, col=cols[3])#high wind
plot(newdata$TEMP,yhatmean,xlab="Air Temp",ylab="Predicted Overall Mean",pch="",col="lightgrey",cex=1.2)
lines(newdata$TEMP[c(101:200)],yhatmean[c(101:200)],lty=1, col=cols[2])#meanwind
lines(newdata$TEMP[c(1:100)],yhatmean[c(1:100)],lty=1, col=cols[1])#low wind
lines(newdata$TEMP[c(201:300)],yhatmean[c(201:300)],lty=1, col=cols[3])#high wind
legend("topleft",lty=1,col=c(cols[3], cols[2], cols[1]), legend=c(" +1 SD", " Wind Mean", " -1 SD"),bty="n",cex=1.2)


plot(OFFdf$TEMP,OFFdf$OFF, xlab="Temp", ylab="off")
lines(newdata$TEMP[c(101:200)],yhatmean[c(101:200)],lty=1)#meanwind


#PANTING as a predictor of time off nest 
# USE DATA df from 0 DataWrangling.R
OFFdf <- df[,c(6:8,10,12)] #Temp, wind, sun exposure, panting per obs block


# NULL MODEL (WITH NO PREDICTORS)
model0zinb<-zeroinfl(OFF~1|1,data = OFFdf, link = "logit", dist = "negbin")

# MODEL WITH MAIN EFFECT OF TEMP 
zinb_TEMP<-zeroinfl(OFF~TEMP|TEMP,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_TEMP)

# MODEL WITH MAIN EFFECT OF PANT 
zinb_PANT<-zeroinfl(OFF~PANT|PANT,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_PANT) #near significant
lrtest(zinb_PANT,model0zinb)
# 

################################################################
## PRESENTING AND INTERPRETING INTERACTIONS IN MIXTURE MODELS ##
################################################################

## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN TEMP AND WIND
zinb_TEMPWIND2<-zeroinfl(OFF~TEMP*PANT|TEMP*PANT,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_TEMPWIND2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_TEMPWIND<-zeroinfl(OFF~TEMP+PANT|TEMP+PANT,data = OFFdf, link = "logit", dist = "negbin")
lrtest(zinb_TEMPWIND2,zinb_TEMPWIND)

# n.s.


#######################
# Including TIMEOFDAY, SEASON in the model
OFFdf <- df[,c(6:8,12,14:15)] #Temp, wind, sun exposure, TIMEOFDAY, SEASON and time off the nest per obs block

#Full model includnig TEMP, WIND, SEXP, DLGT, DAYn
zinb_FULL <- zeroinfl(OFF~.|.,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_FULL)
lrtest(zinb_FULL,model0zinb)
# Backward stepwise variable selection with significance level alpha=0.01.
library("mpath")
fitbe <- be.zeroinfl(zinb_FULL, data=OFFdf, dist="negbin", alpha=0.01, trace=FALSE)
summary(fitbe)


# MODEL WITH MAIN EFFECT OF DLGT 
zinb_DLGT<-zeroinfl(OFF~DLGT|DLGT,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_DLGT)

# MODEL WITH INTERCTION OF DLGT AND TEMP 
zinb_FULL<-zeroinfl(OFF~TEMP+WIND+SEXP+DLGT|TEMP+WIND+SEXP+DLGT,data = PANTdf, link = "logit", dist = "negbin")
zinb_FULL<-zeroinfl(OFF~TEMP+WIND+SEXP+(TEMP*DLGT)|TEMP+WIND+SEXP+(TEMP*DLGT),data = OFFdf, link = "logit", dist = "negbin")
zinb_TEMPDLGT<-zeroinfl(OFF~TEMP*DLGT|TEMP*DLGT,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_TEMPDLGT)
lrtest(zinb_TEMPDLGT,model0zinb)

# MODEL WITH MAIN EFFECT OF DAYn
zinb_DAYn<-zeroinfl(OFF~DAYn|DAYn,data = OFFdf, link = "logit", dist = "negbin")
summary(zinb_DAYn)
zinb_FULL<-zeroinfl(OFF~TEMP+WIND+SEXP+(TEMP*DAYn)|TEMP+WIND+SEXP+(TEMP*DAYn),data = OFFdf, link = "logit", dist = "negbin")
