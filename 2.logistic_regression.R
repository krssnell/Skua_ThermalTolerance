# Analysis from Kasper. 
# 2. Logistc regressionmodels with binary response - panting or not panting
# 2.1 logistic regression
# 2.2 Logistic regression with Random factor *** used in manuscript ***
# 2.3 Bionomial regression of proportional data with random factor
library(lme4)



# USE DATA df from 0 DataWrangling.R
load("df.RData")
df15 <- df[,c(6:8,10)] #Temp, wind, sun exposure, panting per obs block


df15$PANT[df15$PANT!=0] <- 1 #panting now binary 
table(df15$PANT)

# 2.1
fit.full <- glm(PANT ~ WIND + TEMP + SEXP,
                data=df15, family=binomial())
summary(fit.full)

fit.WIND <- glm(PANT ~ WIND,
                data=df15, family=binomial())
summary(fit.WIND)

fit.TEMP <- glm(PANT ~ TEMP,
                data=df15, family=binomial())
summary(fit.TEMP)

fit.SEXP <- glm(PANT ~ SEXP,
                data=df15, family=binomial())
summary(fit.SEXP)

fit.TEMPWIND <- glm(PANT ~ WIND + TEMP,
                data=df15, family=binomial())
summary(fit.TEMPWIND)

fit.TEMPSEXP <- glm(PANT ~  TEMP + SEXP,
                data=df15, family=binomial())
summary(fit.TEMPSEXP)

fit.WINDSEXP <- glm(PANT ~ WIND +  SEXP,
                data=df15, family=binomial())
summary(fit.WINDSEXP)

t(AIC(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP ))

#2.2 #logistic regression with random effect 
# USE DATA df from 0 DataWrangling.R
df15 <- df[,c(4,6:8,10)] #Terr, Temp, wind, sun exposure, panting per obs block
df15$PANT[df15$PANT!=0] <- 1 #panting now binary 
table(df15$PANT)


library(lme4)
fit.full <- glmer(PANT ~ WIND + TEMP + SEXP + (1|TERRID),
                data=df15, family=binomial(link = "logit"))
summary(fit.full)

fit.null <- glmer(PANT ~ 1 + (1|TERRID),   #intercept only model
                  data=df15, family=binomial())
summary(fit.null)

fit.WIND <- glmer(PANT ~ WIND + (1|TERRID),
                data=df15, family=binomial())
summary(fit.WIND)

fit.TEMP <- glmer(PANT ~ TEMP+ (1|TERRID),
                data=df15, family=binomial())
summary(fit.TEMP)

fit.SEXP <- glmer(PANT ~ SEXP+ (1|TERRID),
                data=df15, family=binomial())
summary(fit.SEXP)

fit.TEMPWIND <- glmer(PANT ~ WIND + TEMP+ (1|TERRID),
                    data=df15, family=binomial())
summary(fit.TEMPWIND)

fit.TEMPSEXP <- glmer(PANT ~  TEMP + SEXP+ (1|TERRID),
                    data=df15, family=binomial())
summary(fit.TEMPSEXP)

fit.WINDSEXP <- glmer(PANT ~ WIND +  SEXP+ (1|TERRID),
                    data=df15, family=binomial())
summary(fit.WINDSEXP)

t(AIC(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP, fit.null ))

library(AICcmodavg)
models <- list(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP, fit.null )
aictab(cand.set = models)

#model predictions: 
newdataTEMP <- seq(from=min(df15$TEMP),to=max(df15$TEMP), length.out=100)
pp <- lapply(newdataTEMP, function(j) {
  df15$TEMP <- j
  predict(fit.full, newdata = df15, type = "response")
})
sapply(pp[c(1, 20, 40, 60, 80, 100)], mean)
plotdat <- t(sapply(pp, function(x) {  # get the means with lower and upper quartiles
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
}))
plotdatTEMP <- as.data.frame(cbind(plotdat, newdataTEMP)); colnames(plotdatTEMP) <- c("PredictedProbability", "lower", "upper", "TEMP") # add in predictor values and convert to data frame, and rename

newdataWIND <- seq(from=min(df15$WIND),to=max(df15$WIND), length.out=100)
pp <- lapply(newdataWIND, function(j) {
  df15$WIND <- j
  predict(fit.full, newdata = df15, type = "response")
})
sapply(pp[c(1, 20, 40, 60, 80, 100)], mean)
plotdat <- t(sapply(pp, function(x) {  # get the means with lower and upper quartiles
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
}))
plotdatWIND <- as.data.frame(cbind(plotdat, newdataWIND)); colnames(plotdatWIND) <- c("PredictedProbability", "lower", "upper", "WIND") # add in predictor values and convert to data frame, and rename




newdataSEXP <- seq(from=min(df15$SEXP),to=max(df15$SEXP), length.out=100)
pp <- lapply(newdataSEXP, function(j) {
  df15$SEXP <- j
  predict(fit.full, newdata = df15, type = "response")
})
sapply(pp[c(1, 20, 40, 60, 80, 100)], mean)
plotdat <- t(sapply(pp, function(x) {  # get the means with lower and upper quartiles
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
}))
plotdatSEXP <- as.data.frame(cbind(plotdat, newdataSEXP)); colnames(plotdatSEXP) <- c("PredictedProbability", "lower", "upper", "SEXP") # add in predictor values and convert to data frame, and rename

tiff(file = "SM Fig3a.tiff", units="cm",width=19,height=6,res=500)
par(mfrow = c(1, 3), xpd=NA,cex.axis=0.8, cex.lab=1, mgp=c(1.1,.25,0), tcl=-0.2, mar=c(2,2,2,1), ps=10)

plot(df15$TEMP,df15$PANT, xlab="Air temperature (\u00b0C)", ylab="Panting", pch=20, cex=0.6, col=alpha("black", 0.2))
lines(plotdatTEMP$TEMP,plotdatTEMP$PredictedProbability,lty=1, lwd=1.5)
lines(plotdatTEMP$TEMP,plotdatTEMP$upper,lty=3, col="darkgrey", lwd=1.5)
lines(plotdatTEMP$TEMP,plotdatTEMP$lower,lty=3, col="darkgrey", lwd=1.5)
mtext(text="A", side = 3, adj = -0.05, line = 0, cex=0.8)

plot(df15$WIND,df15$PANT, xlab=expression(Wind~speed~(ms^-1)), ylab="", pch=20, cex=0.6, col=alpha("black", 0.2))
lines(plotdatWIND$WIND,plotdatWIND$PredictedProbability,lty=1, lwd=1.5)
lines(plotdatWIND$WIND,plotdatWIND$upper,lty=3, col="darkgrey", lwd=1.5)
lines(plotdatWIND$WIND,plotdatWIND$lower,lty=3, col="darkgrey", lwd=1.5)
mtext(text="B", side = 3, adj = -0.05, line = 0, cex=0.8)

plot(jitter(df15$SEXP),df15$PANT, xlab="Sun exposure index", ylab="", pch=20, cex=0.6, col=alpha("black", 0.2))
lines(plotdatSEXP$SEXP,plotdatSEXP$PredictedProbability,lty=1, lwd=1.5)
lines(plotdatSEXP$SEXP,plotdatSEXP$upper,lty=3, col="darkgrey", lwd=1.5)
lines(plotdatSEXP$SEXP,plotdatSEXP$lower,lty=3, col="darkgrey", lwd=1.5)
mtext(text="C", side = 3, adj = -0.05, line = 0, cex=0.8)

dev.off()



# 2.3 # Bionomial regression of proportional data with random factor
df15 <- df[,c(4,6:8,10)] #Terr, Temp, wind, sun exposure, panting per obs block
df15$PANT <- df15$PANT / 15 #panting now scaled 0 to 1 as proportion of time panting 
table(df15$PANT)



fit.full <- glmer(PANT ~ WIND + TEMP + SEXP + (1|TERRID),
                  data=df15, family=binomial(link = "logit"))
summary(fit.full)

fit.WIND <- glmer(PANT ~ WIND + (1|TERRID),
                  data=df15, family=binomial())
summary(fit.WIND)

fit.TEMP <- glmer(PANT ~ TEMP+ (1|TERRID),
                  data=df15, family=binomial())
summary(fit.TEMP)

fit.SEXP <- glmer(PANT ~ SEXP+ (1|TERRID),
                  data=df15, family=binomial())
summary(fit.SEXP)

fit.TEMPWIND <- glmer(PANT ~ WIND + TEMP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.TEMPWIND)

fit.TEMPSEXP <- glmer(PANT ~  TEMP + SEXP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.TEMPSEXP)

fit.WINDSEXP <- glmer(PANT ~ WIND +  SEXP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.WINDSEXP)

t(AIC(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP ))


#plot glm - see https://www.statology.org/plot-logistic-regression-in-r/ for plotting in base r
library(ggplot2)
theme_set(theme_bw())
theme_update(text = element_text(size=12),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
ggplot(df15, aes(x=TEMP, y=PANT)) + geom_point() + 
  geom_smooth(method="glm", method.args= list(family ="binomial"))

ggplot(df15, aes(x=WIND, y=PANT)) + geom_point() + 
  geom_smooth(method="glm", method.args= list(family ="binomial"))


ggplot(df15, aes(x=SEXP, y=PANT)) + geom_point() + 
  geom_smooth(method="glm", method.args= list(family ="binomial"))




# USE DATA df from 0 DataWrangling.R
df15 <- df[,c(4,6:8,12)] #Temp, wind, sun exposure, panting per obs block
str(df15)



##nest off - treat also as binary

df15$OFF[df15$OFF!=0] <- 1 #incubation now binary 
table(df15$OFF)


fit.full <- glm(OFF ~ WIND + TEMP + SEXP,
                data=df15, family=binomial())
summary(fit.full)

fit.TEMPWIND <- glm(OFF ~ WIND + TEMP,
                data=df15, family=binomial())
summary(fit.TEMPWIND)

fit.WIND <- glm(OFF ~ WIND,
                data=df15, family=binomial())
summary(fit.WIND)

fit.TEMP <- glm(OFF ~ TEMP,
                data=df15, family=binomial())
summary(fit.TEMP)

fit.SEXP <- glm(OFF ~ SEXP,
                data=df15, family=binomial())
summary(fit.SEXP)

fit.TEMPSEXP <- glm(OFF ~ SEXP + TEMP,
                    data=df15, family=binomial())
summary(fit.TEMPSEXP)

fit.WINDSEXP <- glm(OFF ~ WIND + SEXP,
                    data=df15, family=binomial())
summary(fit.TEMPWIND)

t(AIC(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP ))




library(lme4)
fit.full <- glmer(OFF ~ WIND + TEMP + SEXP + (1|TERRID),
                  data=df15, family=binomial(link = "logit"))
summary(fit.full)

fit.WIND <- glmer(OFF ~ WIND + (1|TERRID),
                  data=df15, family=binomial())
summary(fit.WIND)

fit.TEMP <- glmer(OFF ~ TEMP+ (1|TERRID),
                  data=df15, family=binomial())
summary(fit.TEMP)

fit.SEXP <- glmer(OFF ~ SEXP+ (1|TERRID),
                  data=df15, family=binomial())
summary(fit.SEXP)

fit.TEMPWIND <- glmer(OFF ~ WIND + TEMP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.TEMPWIND)

fit.TEMPSEXP <- glmer(OFF ~  TEMP + SEXP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.TEMPSEXP)

fit.WINDSEXP <- glmer(OFF ~ WIND +  SEXP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.WINDSEXP)

t(AIC(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP ))


#2.3 off nest as proportion 
df15 <- df[,c(4,6:8,12)] #Temp, wind, sun exposure, panting per obs block
str(df15)

##nest off - treat also as binary

df15$OFF <- df15$OFF/15 #incubation now proportional 
table(df15$OFF)

fit.full <- glmer(OFF ~ WIND + TEMP + SEXP + (1|TERRID),
                  data=df15, family=binomial(link = "logit"))
summary(fit.full)

fit.WIND <- glmer(OFF ~ WIND + (1|TERRID),
                  data=df15, family=binomial())
summary(fit.WIND)

fit.TEMP <- glmer(OFF ~ TEMP+ (1|TERRID),
                  data=df15, family=binomial())
summary(fit.TEMP)

fit.SEXP <- glmer(OFF ~ SEXP+ (1|TERRID),
                  data=df15, family=binomial())
summary(fit.SEXP)

fit.TEMPWIND <- glmer(OFF ~ WIND + TEMP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.TEMPWIND)

fit.TEMPSEXP <- glmer(OFF ~  TEMP + SEXP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.TEMPSEXP)

fit.WINDSEXP <- glmer(OFF ~ WIND +  SEXP+ (1|TERRID),
                      data=df15, family=binomial())
summary(fit.WINDSEXP)

t(AIC(fit.full,fit.TEMP, fit.WIND,fit.SEXP,fit.TEMPWIND, fit.TEMPSEXP, fit.WINDSEXP ))


#plot glm - see https://www.statology.org/plot-logistic-regression-in-r/ for plotting in base r
library(ggplot2)
theme_set(theme_bw())
theme_update(text = element_text(size=12),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
ggplot(df15, aes(x=TEMP, y=OFF)) + geom_point() + 
  geom_smooth(method="glm", method.args= list(family ="binomial"))

ggplot(df15, aes(x=WIND, y=OFF)) + geom_point() + 
  geom_smooth(method="glm", method.args= list(family ="binomial"))



######################
###########
#Panting as a predictor of time off the nest

# USE DATA df from 0 DataWrangling.R
df15 <- df[,c(6:8,10,12)] #Temp, wind, sun exposure, panting per obs block

##nest off - treat also as binary

df15$OFF[df15$OFF!=0] <- 1 #panting now binary 
table(df15$OFF)

fit.full <- glm(OFF ~ . ,
                data=df15, family=binomial())
summary(fit.full)
AIC(fit.full)

# 2.3 terr ID as random and GLMM using proportional data
# USE DATA df from 0 DataWrangling.R
df15 <- df[,c(4,6:8,10,12)] #Temp, wind, sun exposure, panting per obs block

df15$OFF <- df15$OFF/15 #incubation now proportional 
table(df15$OFF)

fit.full <- glmer(OFF ~ TEMP + WIND + SEXP + PANT + (1|TERRID) ,
                data=df15, family=binomial())
summary(fit.full)
AIC(fit.full)

# 2.3 terr ID as random and GLMM using proportional data
# USE DATA df from 0 DataWrangling.R
df15 <- df[,c(4,6:8,10,12)] #Temp, wind, sun exposure, panting per obs block

df15$PANT <- df15$PANT/15 #incubation now proportional 
table(df15$PANT)

fit.full <- glmer(PANT ~ TEMP + WIND + SEXP + OFF + (1|TERRID) ,
                  data=df15, family=binomial())
summary(fit.full)
AIC(fit.full)
