library(MASS)  # used to fit negative binomial distribution
library(pscl)  # used to fit zero-inflated and hurdle models
library(lmtest) # used to perform likelihood ratio tests
library(multcomp) # used to perform multiple comparisons
library(boot) # used for bootstrapping CI of zinb models
library(scales) # sets transparency

set.seed(1408)
##############
# USE DATA df from 0 DataWrangling.R
load("df.RData")
PANTdf <- df[,c(6:8,10)] #Temp, wind, sun exposure, panting per obs block


###########################################
##  MODEL COMPARISON AND GOODNESS-OF-FIT ##
###########################################

# FIT COUNT MODELS WITH DIFFERENT DISTRIBUTIONAL ASSUMPTIONS
# Below a Poisson, negative binomial, zero-inflated Poisson, zero-inflated negative binomial,
# Poisson logit hurdle and negative binomial logit hurdle model are fitted, respectively.
# By using the . for the predictor function, all available predictors 
# (education level and anxious attachment level) are included as main effects.

model1p<- glm(PANT ~ .,data= PANTdf,family=poisson)
model1nb<- glm.nb(PANT  ~ . , data=PANTdf)
model1zip<-zeroinfl(PANT~.|.,data = PANTdf, link = "logit", dist = "poisson")
model1zinb<-zeroinfl(PANT~.|.,data = PANTdf, link = "logit", dist = "negbin")
model1plh<-hurdle(PANT~.|.,data=PANTdf,zero.dist="binomial",link="logit",dist="poisson")
model1nblh<-hurdle(PANT~.|.,data=PANTdf,zero.dist="binomial",link="logit",dist="negbin", family="binomial")
summary(model1nblh)
# CALCULATE FOR EACH MODEL THE PREDICTED FREQUENCIES FOR EACH OUTCOME
estobsk<-function(k){
cbind(sum(dpois(k,fitted(model1p))),sum(dnbinom(k,mu=fitted(model1nb),size=model1nb$theta)),sum(predict(model1zip,type="prob")[,1+k]),sum(predict(model1zinb,type="prob")[,1+k]),
sum(predict(model1plh,type="prob")[,1+k]),sum(predict(model1nblh,type="prob")[,1+k]))}

estobs<-sapply(0:15,estobsk)

# CREATE A HISTOGRAM CONTRASTING OBSERVED WITH PREDICTED FREQUENCIES (FIGURE 1 in TUTORIAL PAPER)
# This is done for each of the fitted models above,
# except for the fitted hurdle models that yield (almost) identical predicted frequencies.
# The breaks in the hist-function were explicitly specified as typically this function bins zero and one together.


hist(PANTdf$PANT,xlab="Number mins panting observed",ylab="Absolute Frequency",breaks=seq(-0.5,15.5,by=1),main="",cex=1.2)
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

#for data with 15 removed 
estobs<-sapply(0:14,estobsk)

# CREATE A HISTOGRAM CONTRASTING OBSERVED WITH PREDICTED FREQUENCIES (FIGURE 1 in TUTORIAL PAPER)
# This is done for each of the fitted models above,
# except for the fitted hurdle models that yield (almost) identical predicted frequencies.
# The breaks in the hist-function were explicitly specified as typically this function bins zero and one together.

hist(PANTdf$PANT,xlab="Number mins panting observed",ylab="Absolute Frequency",breaks=seq(-0.5,15.5,by=1),main="",cex=1.2)
points(0:14-0.1,estobs[1,],col=1,lwd=2,pch=2,lty=2)
lines(0:14-0.1,estobs[1,],col=1,lwd=2,pch=2,lty=2)
points(0:14-0.05,estobs[2,],col=1,lwd=2,pch=3,lty=3)
lines(0:14-0.05,estobs[2,],col=1,lwd=2,pch=3,lty=3)
points(0:14+0.05,estobs[3,],col=1,lwd=2,pch=1,lty=2)
lines(0:14+0.05,estobs[3,],col=1,lwd=2,pch=1,lty=2)
points(0:14+0.1,estobs[4,],col=1,lwd=2,pch=4,lty=3)
lines(0:14+0.1,estobs[4,],col=1,lwd=2,pch=4,lty=3)
points(0:14+0.05,estobs[5,],col=3,lwd=2,pch=1,lty=2)
lines(0:14+0.05,estobs[5,],col=3,lwd=2,pch=1,lty=2)
points(0:14+0.1,estobs[6,],col=3,lwd=2,pch=4,lty=3)
lines(0:14+0.1,estobs[6,],col=3,lwd=2,pch=4,lty=3)
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

## FIT MAIN EFFECT MODELS FOR EDUCATION LEVEL  ##
# As from the previous step the zero-inflated negative binomial yielded the best fit,
# all further models use these underlying models.
# Note that the notation y~x1|z1 implies the count component model y~x1 conditional on (|) the hurdle
# or zero-inflated model y~z1 for the zero-component, in our example we always use x1=z1

# NULL MODEL (WITH NO PREDICTORS)
model0zinb<-zeroinfl(PANT~1|1,data = PANTdf, link = "logit", dist = "negbin")
model0nblh<-hurdle(PANT~1|1,data = PANTdf,zero.dist="binomial",link="logit",dist="negbin")

# MODEL WITH MAIN EFFECT OF TEMP 
zinb_TEMP<-zeroinfl(PANT~TEMP|TEMP,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_TEMP)
nblh_TEMP<-hurdle(PANT~TEMP|TEMP,data = PANTdf,zero.dist="binomial",link="logit",dist="negbin")
summary(nblh_TEMP)

# MODEL WITH MAIN EFFECT OF WIND
# similar approach can be used to explore to effect of WIND speed alone
zinb_WIND<-zeroinfl(PANT~WIND|WIND,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_WIND)
nblh_WIND<-hurdle(PANT~WIND|WIND,data = PANTdf,zero.dist="binomial",link="logit",dist="negbin")
summary(nblh_WIND)

# MODEL WITH MAIN EFFECT OF SUN EXPOSURE
# similar approach can be used to explore to effect of WIND speed alone
zinb_SEXP<-zeroinfl(PANT~SEXP|SEXP,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_SEXP)
nblh_SEXP<-hurdle(PANT~SEXP|SEXP,data = PANTdf,zero.dist="binomial",link="logit",dist="negbin")
summary(nblh_WIND)

############################
## MULTIPLE TESTING ISSUE ##
############################

# Use LRT to jointly test that both the impact on the zero-component and count-component are zero
lrtest(zinb_TEMP,model0zinb)
lrtest(nblh_TEMP,model0nblh)

# Use the multiplicity procedure by HOTHORN ET AL.
# first define the contrast matrix for the hypotheses of interest
# this needs to be done manually here
contrastmat<-matrix(0,nrow=2,ncol=length(coef(zinb_TEMP)))
colnames(contrastmat)<-names(coef(zinb_TEMP))
rownames(contrastmat)<-names(coef(zinb_TEMP))[c(2,4)]
contrastmat[1,2]<-contrastmat[2,4]<-1
mc_zinb_TEMP<-glht(zinb_TEMP,linfct=contrastmat)
mc_nblh_ed<-glht(nblh_ed,linfct=contrastmat)
summary(mc_zinb_TEMP)


# Similarly, one can assess the effect of wind on the PANT-outcome

lrtest(zinb_WIND,model0zinb)

# Similarly, one can assess the effect of sun exposure on the PANT-outcome
lrtest(zinb_SEXP,model0zinb)

contrastmat<-matrix(0,nrow=2,ncol=length(coef(zinb_anx)))
colnames(contrastmat)<-names(coef(zinb_anx))
rownames(contrastmat)<-names(coef(zinb_anx))[c(2,4)]
contrastmat[1,2]<-contrastmat[2,4]<-1
mc_zinb_anx<-glht(zinb_anx,linfct=contrastmat)
mc_nblh_anx<-glht(nblh_anx,linfct=contrastmat)
summary(mc_zinb_anx)
summary(mc_nblh_anx)

################################################################
## PRESENTING AND INTERPRETING INTERACTIONS IN MIXTURE MODELS ##
################################################################

## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN TEMP AND WIND
zinb_TEMPWIND2<-zeroinfl(PANT~TEMP*WIND|TEMP*WIND,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_TEMPWIND2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_TEMPWIND<-zeroinfl(PANT~TEMP+WIND|TEMP+WIND,data = PANTdf, link = "logit", dist = "negbin")
lrtest(zinb_TEMPWIND2,zinb_TEMPWIND)



## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN TEMP AND SUN EXPOSURE
zinb_TEMPSEXP2<-zeroinfl(PANT~TEMP*SEXP|TEMP*SEXP,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_TEMPSEXP2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_TEMPSEXP<-zeroinfl(PANT~TEMP+SEXP|TEMP+SEXP,data = PANTdf, link = "logit", dist = "negbin")
lrtest(zinb_TEMPSEXP2,zinb_TEMPSEXP)



## ZERO-INFLATED MODEL WITH INTERACTION BETWEEN WIND AND SUN EXPOSURE
zinb_WINDSEXP2<-zeroinfl(PANT~WIND*SEXP|WIND*SEXP,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_WINDSEXP2)

# perform LRT against model with only main effects to assess overall significance of the interaction terms
zinb_WINDSEXP<-zeroinfl(PANT~WIND+SEXP|WIND+SEXP,data = PANTdf, link = "logit", dist = "negbin")
lrtest(zinb_WINDSEXP2,zinb_WINDSEXP)


#Full model includnig TEMP, WIND, and SEXP
zinb_FULL <- zeroinfl(PANT~.|.,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_FULL)
lrtest(zinb_FULL,model0zinb)
# Backward stepwise variable selection with significance level alpha=0.01.
library("mpath")
fitbe <- be.zeroinfl(zinb_FULL, data=PANTdf, dist="negbin", alpha=0.01, trace=FALSE)
summary(fitbe)


# graphical presentation of effect TEMP on the zero-inflated model and the mean (FIGURE 3 IN TUTORIAL PAPER)
newdataTEMP <- seq(from=min(PANTdf$TEMP),to=max(PANTdf$TEMP), length.out=100)

yhatmean<-predict(zinb_TEMP,newdata=data.frame(TEMP = newdataTEMP),type="response")
yhatcount<-predict(zinb_TEMP,newdata=data.frame(TEMP = newdataTEMP),type="count")
yhatzero<-predict(zinb_TEMP,newdata=data.frame(TEMP = newdataTEMP),type="zero")

yhatTEMP <- cbind.data.frame(newdataTEMP,yhatmean,yhatcount,yhatzero)
#CI need to be calculated using bootsrapping https://stackoverflow.com/questions/58002334/how-to-add-confidence-intervals-for-a-zero-inflated-negative-binomial-function-i
stat <- function(df, inds) {
  model <- formula(PANT~TEMP|TEMP)
  predict(
    zeroinfl(model, dist = "negbin", link = "logit", data = PANTdf[inds, ]),
    newdata = PANTdf)
}

res <- boot(PANTdf, stat, R = 100)
all.equal(res$t0, predict(zinb_TEMP))#verify predicted mean are the same
yhatCI <- setNames(as.data.frame(t(sapply(1:nrow(PANTdf), function(row)
  boot.ci(res, conf = 0.95, type = "basic", index = row)$basic[, 4:5]))),
  c("lower", "upper"))
yhatCI <- cbind.data.frame(PANTdf, response = predict(zinb_TEMP), yhatCI)
TEMPyhatCI <- yhatCI[order(yhatCI$TEMP),]


par(mfrow=c(1,3))
plot(newdataTEMP,yhatcount,xlab="Air Temp",ylab=expression(paste(mu[i])),pch="",cex=1.2)
lines(newdataTEMP,yhatcount,lty=1)
plot(newdataTEMP,1-yhatzero,xlab="Air Temp",,ylab=expression(paste(1-p[i])),pch="",cex=1.2)
lines(newdataTEMP,1-yhatzero,lty=1)
plot(newdataTEMP,yhatmean,xlab="Air Temp",ylab="Predicted Overall Mean",pch="",cex=1.2)
lines(newdataTEMP,yhatmean,lty=1)

plot(PANTdf$TEMP,PANTdf$PANT, xlab="Temp", ylab="Pant", pch=20, cex=0.6)
lines(newdataTEMP,yhatmean,lty=1)
lines(yhatCI$TEMP,yhatCI$upper,lty=3, col="darkgrey")
lines(yhatCI$TEMP,yhatCI$lower,lty=3, col="darkgrey")


# graphical presentation of effect WIND on the zero-inflated model and the mean (FIGURE 3 IN TUTORIAL PAPER)
newdataWIND <- seq(from=min(PANTdf$WIND),to=max(PANTdf$WIND), length.out=100)

yhatmean<-predict(zinb_WIND,newdata=data.frame(WIND = newdataWIND),type="response")
yhatcount<-predict(zinb_WIND,newdata=data.frame(WIND = newdataWIND),type="count")
yhatzero<-predict(zinb_WIND,newdata=data.frame(WIND = newdataWIND),type="zero")

yhatWIND <- cbind.data.frame(newdataWIND,yhatmean,yhatcount,yhatzero)
#CI by bootstrapping
stat <- function(df, inds) {
  model <- formula(PANT~WIND|WIND)
  predict(
    zeroinfl(model, dist = "negbin", link = "logit", data = PANTdf[inds, ]),
    newdata = PANTdf)
}
res <- boot(PANTdf, stat, R = 100)
all.equal(res$t0, predict(zinb_WIND))#verify predicted mean are the same
yhatCI <- setNames(as.data.frame(t(sapply(1:nrow(PANTdf), function(row)
  boot.ci(res, conf = 0.95, type = "basic", index = row)$basic[, 4:5]))),
  c("lower", "upper"))
yhatCI <- cbind.data.frame(PANTdf, response = predict(zinb_WIND), yhatCI)
WINDyhatCI <- yhatCI[order(yhatCI$WIND),]



par(mfrow=c(1,3))
plot(newdataWIND,yhatcount,xlab="WIND",ylab=expression(paste(mu[i])),pch="",cex=1.2)
lines(newdataWIND,yhatcount,lty=1)
plot(newdataWIND,1-yhatzero,xlab="WIND",,ylab=expression(paste(1-p[i])),pch="",cex=1.2)
lines(newdataWIND,1-yhatzero,lty=1)
plot(newdataWIND,yhatmean,xlab="WIND",ylab="Predicted Overall Mean",pch="",cex=1.2)
lines(newdataWIND,yhatmean,lty=1)


plot(PANTdf$WIND,PANTdf$PANT, xlab="Wind", ylab="Pant", pch=20, cex=0.6)
lines(newdataWIND,yhatmean,lty=1)
lines(yhatCI$WIND,yhatCI$upper,lty=3, col="darkgrey")
lines(yhatCI$WIND,yhatCI$lower,lty=3, col="darkgrey")

# graphical presentation of effect SEXP on the zero-inflated model and the mean (FIGURE 3 IN TUTORIAL PAPER)
newdataSEXP <- seq(from=min(PANTdf$SEXP),to=max(PANTdf$SEXP), length.out=100)

yhatmean<-predict(zinb_SEXP,newdata=data.frame(SEXP = newdataSEXP),type="response")
yhatcount<-predict(zinb_SEXP,newdata=data.frame(SEXP = newdataSEXP),type="count")
yhatzero<-predict(zinb_SEXP,newdata=data.frame(SEXP = newdataSEXP),type="zero")

yhatSEXP <- cbind.data.frame(newdataSEXP,yhatmean,yhatcount,yhatzero)
#CI by bootstrapping
stat <- function(df, inds) {
  model <- formula(PANT~SEXP|SEXP)
  predict(
    zeroinfl(model, dist = "negbin", link = "logit", data = PANTdf[inds, ]),
    newdata = PANTdf)
}
res <- boot(PANTdf, stat, R = 100)
all.equal(res$t0, predict(zinb_SEXP))#verify predicted mean are the same
yhatCI <- setNames(as.data.frame(t(sapply(1:nrow(PANTdf), function(row)
  boot.ci(res, conf = 0.95, type = "basic", index = row)$basic[, 4:5]))),
  c("lower", "upper"))
yhatCI <- cbind.data.frame(PANTdf, response = predict(zinb_SEXP), yhatCI)
SEXPyhatCI <- yhatCI[order(yhatCI$SEXP),]

par(mfrow=c(1,3))
plot(newdataSEXP,yhatcount,xlab="SEXP",ylab=expression(paste(mu[i])),pch="",cex=1.2)
lines(newdataSEXP,yhatcount,lty=1)
plot(newdataSEXP,1-yhatzero,xlab="SEXP",,ylab=expression(paste(1-p[i])),pch="",cex=1.2)
lines(newdataSEXP,1-yhatzero,lty=1)
plot(newdataSEXP,yhatmean,xlab="SEXP",ylab="Predicted Overall Mean",pch="",cex=1.2)
lines(newdataSEXP,yhatmean,lty=1)


plot(jitter(PANTdf$SEXP),jitter(PANTdf$PANT), xlab="SEXP", ylab="Pant", cex=0.6, pch=20)
lines(newdataSEXP,yhatmean,lty=1)
lines(yhatCI$SEXP,yhatCI$upper,lty=3, col="darkgrey")
lines(yhatCI$SEXP,yhatCI$lower,lty=3, col="darkgrey")


##############################
##############################
tiff(file = "Fig1a.tiff", units="cm",width=19,height=6,res=500)
par(mfrow = c(1, 3), xpd=NA,cex.axis=0.8, cex.lab=1, mgp=c(1.1,.25,0), tcl=-0.2, mar=c(2,2,2,1), ps=10)


plot(TEMPyhatCI$TEMP,TEMPyhatCI$PANT, xlab="Air temperature (\u00b0C)", ylab="Panting", pch=20, cex=0.6, col=alpha("black", 0.2))
lines(TEMPyhatCI$TEMP,TEMPyhatCI$response,lty=1, lwd=1.5)
lines(TEMPyhatCI$TEMP,TEMPyhatCI$upper,lty=3, col="darkgrey", lwd=1.5)
lines(TEMPyhatCI$TEMP,TEMPyhatCI$lower,lty=3, col="darkgrey", lwd=1.5)
mtext(text="A", side = 3, adj = -0.05, line = 0, cex=0.8)

plot(WINDyhatCI$WIND,WINDyhatCI$PANT, xlab=expression(Wind~speed~(ms^-1)), ylab="", pch=20, cex=0.6, col=alpha("black", 0.2))
lines(WINDyhatCI$WIND,WINDyhatCI$response,lty=1, lwd=1.5)
lines(WINDyhatCI$WIND,WINDyhatCI$upper,lty=3, col="darkgrey", lwd=1.5)
lines(WINDyhatCI$WIND,WINDyhatCI$lower,lty=3, col="darkgrey", lwd=1.5)
mtext(text="B", side = 3, adj = -0.05, line = 0, cex=0.8)

plot(jitter(SEXPyhatCI$SEXP),jitter(SEXPyhatCI$PANT), xlab="Sun exposure index", ylab="", pch=20, cex=0.6, col=alpha("black", 0.2))
lines(SEXPyhatCI$SEXP,SEXPyhatCI$response,lty=1, lwd=1.5)
lines(SEXPyhatCI$SEXP,SEXPyhatCI$upper,lty=3, col="darkgrey", lwd=1.5)
lines(SEXPyhatCI$SEXP,SEXPyhatCI$lower,lty=3, col="darkgrey", lwd=1.5)
mtext(text="C", side = 3, adj = -0.05, line = 0, cex=0.8)

dev.off()

##############################
##############################

tiff(file = "SM Fig2.tiff", units="cm",width=19,height=22,res=500)
par(mfrow = c(3, 3), xpd=NA,cex.axis=0.8, cex.lab=1, mgp=c(1.1,.25,0), tcl=-0.2, mar=c(2,2,2,1), ps=10)

plot(yhatTEMP$newdataTEMP,yhatTEMP$yhatcount,xlab="Air temperature (\u00b0C)",ylab=expression(paste(mu[i])),pch="",cex=1.2)
lines(yhatTEMP$newdataTEMP,yhatTEMP$yhatcount,lty=1)
mtext(text="Ai", side = 3, adj = -0.05, line = 0.2, cex=0.8)
plot(yhatTEMP$newdataTEMP,1-yhatTEMP$yhatzero,xlab="Air temperature (\u00b0C)",,ylab=expression(paste(1-p[i])),pch="",cex=1.2)
lines(yhatTEMP$newdataTEMP,1-yhatTEMP$yhatzero,lty=1)
mtext(text="Aii", side = 3, adj = -0.05, line = 0.2, cex=0.8)
plot(yhatTEMP$newdataTEMP,yhatTEMP$yhatmean,xlab="Air temperature (\u00b0C)",ylab="Panting (Predicted Overall Mean)",pch="",cex=1.2)
lines(yhatTEMP$newdataTEMP,yhatTEMP$yhatmean,lty=1)
mtext(text="Aiii", side = 3, adj = -0.05, line = 0.2, cex=0.8)

plot(yhatWIND$newdataWIND,yhatWIND$yhatcount,xlab=expression(Wind~speed~(ms^-1)),ylab=expression(paste(mu[i])),pch="",cex=1.2)
lines(yhatWIND$newdataWIND,yhatWIND$yhatcount,lty=1)
mtext(text="Bi", side = 3, adj = -0.05, line = 0.2, cex=0.8)
plot(yhatWIND$newdataWIND,1-yhatWIND$yhatzero,xlab=expression(Wind~speed~(ms^-1)),ylab=expression(paste(1-p[i])),pch="",cex=1.2)
lines(yhatWIND$newdataWIND,1-yhatWIND$yhatzero,lty=1)
mtext(text="Bii", side = 3, adj = -0.05, line = 0.2, cex=0.8)
plot(yhatWIND$newdataWIND,yhatWIND$yhatmean,xlab=expression(Wind~speed~(ms^-1)),ylab="Panting (Predicted Overall Mean)",pch="",cex=1.2)
lines(yhatWIND$newdataWIND,yhatWIND$yhatmean,lty=1)
mtext(text="Biii", side = 3, adj = -0.05, line = 0.2, cex=0.8)

plot(newdataSEXP,yhatcount,xlab="Sun exposure index",ylab=expression(paste(mu[i])),pch="",cex=1.2)
lines(newdataSEXP,yhatcount,lty=1)
mtext(text="Ci", side = 3, adj = -0.05, line = 0.2, cex=0.8)
plot(newdataSEXP,1-yhatzero,xlab="Sun exposure index",,ylab=expression(paste(1-p[i])),pch="",cex=1.2)
lines(newdataSEXP,1-yhatzero,lty=1)
mtext(text="Cii", side = 3, adj = -0.05, line = 0.2, cex=0.8)
plot(newdataSEXP,yhatmean,xlab="Sun exposure index",ylab="Panting (Predicted Overall Mean)",pch="",cex=1.2)
lines(newdataSEXP,yhatmean,lty=1)
mtext(text="Ciii", side = 3, adj = -0.05, line = 0.2, cex=0.8)

dev.off()


##############################
##############################

# graphical presentation of effect WIND on the negative binomaial logit hurdle model model and the mean (FIGURE 3 IN TUTORIAL PAPER)
newdataWIND <- seq(from=min(couple$WIND),to=max(couple$WIND), length.out=100)

yhatmean<-predict(nblh_WIND,newdata=data.frame(WIND = newdataWIND),type="response")
yhatcount<-predict(nblh_WIND,newdata=data.frame(WIND = newdataWIND),type="count")
yhatzero<-predict(nblh_WIND,newdata=data.frame(WIND = newdataWIND),type="zero")

par(mfrow=c(1,3))
plot(newdataWIND,yhatcount,xlab="WIND",ylab=expression(paste(mu[i],"*")),pch="",cex=1.2)
lines(newdataWIND,yhatcount,lty=1)
plot(newdataWIND,1-yhatzero,xlab="WIND",,ylab=expression(paste(1-p[i],"*")),pch="",cex=1.2)
lines(newdataWIND,1-yhatzero,lty=1)
plot(newdataWIND,yhatmean,xlab="WIND",ylab="Predicted Overall Mean",pch="",cex=1.2)
lines(newdataWIND,yhatmean,lty=1)


plot(couple$WIND,couple$PANT, xlab="Wind", ylab="Pant")
lines(newdataWIND,yhatmean,lty=1)



#######################
# Including TIMEOFDAY, SEASON in the model
PANTdf <- df[,c(6:8,10,14:15)] #Temp, wind, sun exposure, TIMEOFDAY, SEASON and panting per obs block
PANTdf <- df[,c(6:8,10,14)] #Temp, wind, sun exposure, TIMEOFDAY, and panting per obs block

#Full model includnig TEMP, WIND, SEXP, DLGT, DAYn
zinb_FULL <- zeroinfl(PANT~.|.,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_FULL)
lrtest(zinb_FULL,model0zinb)
# Backward stepwise variable selection with significance level alpha=0.01.
library("mpath")
fitbe <- be.zeroinfl(zinb_FULL, data=PANTdf, dist="negbin", alpha=0.01, trace=FALSE)
summary(fitbe)
lrtest(fitbe,model0zinb)

# MODEL WITH MAIN EFFECT OF DLGT 
zinb_DLGT<-zeroinfl(PANT~DLGT|DLGT,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_DLGT)

# MODEL WITH INTERCTION OF DLGT AND TEMP 
zinb_FULL<-zeroinfl(PANT~TEMP+WIND+SEXP+DLGT|TEMP+WIND+SEXP+DLGT,data = PANTdf, link = "logit", dist = "negbin")
zinb_FULL<-zeroinfl(PANT~TEMP+WIND+SEXP+(TEMP*DLGT)|TEMP+WIND+SEXP+(TEMP*DLGT),data = PANTdf, link = "logit", dist = "negbin")
zinb_TEMPDLGT<-zeroinfl(PANT~(TEMP*DLGT)|(TEMP*DLGT),data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_TEMPDLGT)
lrtest(zinb_TEMPDLGT,model0zinb)

# MODEL WITH MAIN EFFECT OF DAYn
zinb_DAYn<-zeroinfl(PANT~DAYn|DAYn,data = PANTdf, link = "logit", dist = "negbin")
summary(zinb_DAYn)




#######################
## pseudo-R2 measures
# McFadden’s pseudo-R squared can be also used with (2.) binomial response logistic regression model BUT not with random effects as far as I can find. 

pR2(zinb_FULL)

pR2(model0zinb)# NULL MODEL (WITH NO PREDICTORS)
pR2(zinb_TEMP)# MODEL WITH MAIN EFFECT OF TEMP 
pR2(zinb_WIND)# MODEL WITH MAIN EFFECT OF WIND
pR2(zinb_SEXP)# MODEL WITH MAIN EFFECT OF SUN EXPOSURE

library(performance)
r2_zeroinflated(zinb_TEMP)
r2_zeroinflated(zinb_WIND)
r2_zeroinflated(zinb_SEXP)
