## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo = FALSE, warning=FALSE---------------------------------------------
library(cdfquantreg)
library(MASS)

data(cdfqrExampleData)

## ---- warning=FALSE-----------------------------------------------------------
# Ambulance-arrival vs walk-in split shows 683 walk-ins and 211 ambulance-arrivals
table(yoon$Ambulance)
#
# Length of stay (in hours) is skewed, but the log of it isn't:
hist(yoon$LOSh, breaks = 50, main = "", xlab = "LOS", col = "gray")
hist(yoon$LOSh, breaks = 50, main = "", xlab = "ln(LOS)", col = "gray")
hist(log(yoon$LOSh), breaks = 50, main = "", xlab = "ln(LOS)", col = "gray")
loglosh <- log(yoon$LOSh)
# 

## ---- warning=FALSE-----------------------------------------------------------
# Examine alternative distribution models' goodness-of-fit:
fit1 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "arcsinh", sd = "arcsinh", inner = FALSE, version = "V", data = yoon); logLik(fit1)
fit2 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "arcsinh", sd = "cauchy", inner = FALSE, version = "V", data = yoon); logLik(fit2)
fit3 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "cauchit", sd = "arcsinh", inner = FALSE, version = "V", data = yoon); logLik(fit3)
fit4 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "cauchit", sd = "arcsinh", inner = FALSE, version = "W", data = yoon); logLik(fit4)
fit5 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "cauchit", sd = "cauchy", inner = FALSE, version = "V", data = yoon); logLik(fit5)
fit6 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "cauchit", sd = "cauchy", inner = FALSE, version = "W", data = yoon); logLik(fit6)
fit7 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "t2", sd = "t2", inner = FALSE, version = "V", data = yoon); logLik(fit7)
fit8 <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "t2", sd = "t2", inner = FALSE, version = "W", data = yoon); logLik(fit8)
# 
# The Cauchit-ArcSinh outer-W model fits best.
summary(fit4)
#
# Does a model without the effects in the dispersion submodel fit as well?
# The answer is no:
fit4b <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|1, fd = "cauchit", sd = "arcsinh", inner = FALSE, version = "W", data = yoon); anova(fit4b, fit4)
#
# How about a 3-parameter model?
fit4c <- cdfquantregFT(pregptriage ~ Ambulance + loglosh|Ambulance + loglosh, fd = "cauchit", sd = "arcsinh", mu.fo ~ Ambulance + loglosh, inner = FALSE, version = "W", data = yoon)
c(AIC(fit4),AIC(fit4c))
# The AIC values are nearly identical.  Moroever, there are no significant $\mu$ effects:
summary(fit4c)
#
# Finally, what about a 2-parameter model with moderator effects?
# No significant improvement in model fit:
fit4d <- cdfquantregFT(pregptriage ~ Ambulance*loglosh|Ambulance*loglosh, fd = "cauchit", sd = "arcsinh", inner = FALSE, version = "W", data = yoon); anova(fit4, fit4d)
# 

## ---- warning=FALSE-----------------------------------------------------------
# Parameter-estimate correlations
cov2cor(vcov(fit4))
#
# Fit the marginal distributions to the walk-ins:
uniqdat <- sort(unique(yoon$pregptriage[yoon$Ambulance==0]))
densfit <- c(rep(0,length(uniqdat)))
for (i in 1:length(uniqdat)) {
  densfit[i] <- pdfft(uniqdat[i], sigma = exp(coef(fit4)[4] + mean(loglosh)*coef(fit4)[6]), theta = coef(fit4)[1] + mean(loglosh)*coef(fit4)[3], fd = "cauchit", sd = "arcsinh", mu = NULL, inner = FALSE, version = "W")
}
truehist(yoon$pregptriage[yoon$Ambulance==0], nbins = 100, main = "Walk-ins", xlab = "proportion", ylab = "pdf")
lines(uniqdat,densfit,lwd = 2)
# Fit the marginal distributions to the ambulance-arrivals:
uniqdat <- sort(unique(yoon$pregptriage[yoon$Ambulance==1]))
densfit <- c(rep(0,length(uniqdat)))
for (i in 1:length(uniqdat)) {
  densfit[i] <- pdfft(uniqdat[i], sigma = exp(coef(fit4)[4] + coef(fit4)[5] + mean(loglosh)*coef(fit4)[6]), theta = coef(fit4)[1] + coef(fit4)[2] + mean(loglosh)*coef(fit4)[3], fd = "cauchit", sd = "arcsinh", mu = NULL, inner = FALSE, version = "W")
}
truehist(yoon$pregptriage[yoon$Ambulance==1], nbins = 100, main = "Ambulance-arrivals", xlab = "proportion", ylab = "pdf")
lines(uniqdat,densfit,lwd = 2)
#
# The fitted cdf should be monotonically related to the dependent variable.
# How well does our model do?
fittedcdf <- c(rep(0,length(yoon$pregptriage)))
for (i in 1: length(yoon$pregptriage)) {
fittedcdf[i] <- cdfft(yoon$pregptriage[i], sigma = exp(coef(fit4)[4] + yoon$Ambulance[i]*coef(fit4)[5] + loglosh[i]*coef(fit4)[6]), theta = coef(fit4)[1] + yoon$Ambulance[i]*coef(fit4)[2] + loglosh[i]*coef(fit4)[3], fd = "cauchit", sd = "arcsinh", mu = NULL, inner = FALSE, version = "W")
   }
cor(fittedcdf, yoon$pregptriage, method = "spearman")
#

