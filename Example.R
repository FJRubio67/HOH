## ----message=FALSE------------------------------------------------------------------------------------------------------------------
#######################################################################################
# Data preparation
#######################################################################################

rm(list=ls())

# Required packages
library(deSolve)
library(survival)
library(ggplot2)
#library(devtools)
#install_github("FJRubio67/HazReg")
library(HazReg)
library(spBayes)
library(Rtwalk)
library(knitr)
library(bshazard)


## -----------------------------------------------------------------------------------------------------------------------------------
source("routinesHO.R")


## -----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Data preparation
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


head(rotterdam)

dim(rotterdam)

# Kaplan-Meier estimator for the survival times
km <- survfit(Surv(rotterdam$dtime/365.24, rotterdam$death) ~ 1)

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")



# New data frame: logical status, time in years, survival times sorted
df <- data.frame(time = rotterdam$dtime, status = rotterdam$death)
df$status <- as.logical(rotterdam$death)
df$time <- df$time/365.24

df <- df[order(df$time),]

# Required quantities
status <- as.logical(df$status)
t_obs <- df$time[status]
survtimes <- df$time


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Maximum Likelihood Analysis
#==================================================================================================


#--------------------------------------------------------------------------------------------------
# Fitting a Weibull distribution
#--------------------------------------------------------------------------------------------------

# Initial value
initW <- c(0,0)

# Optimisation step
OPTW <- GHMLE(initW, survtimes, status, hstr = "baseline", dist = "Weibull", method = "nlminb", maxit = 10000)

MLEW <- exp(OPTW$OPT$par)

# Fitted Weibull hazard
fithw <- Vectorize(function(t) hweibull( t, exp(OPTW$OPT$par[1]), exp(OPTW$OPT$par[2]) ) )

# Fitted Weibull cumulative hazard
fitchw <- Vectorize(function(t) chweibull( t, exp(OPTW$OPT$par[1]), exp(OPTW$OPT$par[2]) ) )

# Fitted Weibull survival
fitsw <- Vectorize(function(t) exp(-chweibull( t, exp(OPTW$OPT$par[1]), exp(OPTW$OPT$par[2]) ) ))

# AIC
AICW <- 2*OPTW$OPT$objective + 2*length(OPTW$OPT$par)

# BIC
BICW <- 2*OPTW$OPT$objective + length(OPTW$OPT$par)*log(length(survtimes))


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Maximum Likelihood Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Fitting a PGW distribution
#--------------------------------------------------------------------------------------------------

# Initial value
initPGW <- c(0,0,0)

# Optimisation step
OPTPGW <- GHMLE(initPGW, survtimes, status, hstr = "baseline", dist = "PGW", 
                method = "nlminb", maxit = 10000)

MLEPGW <- exp(OPTPGW$OPT$par)

# Fitted Weibull hazard
fithpgw <- Vectorize(function(t) hpgw( t, MLEPGW[1], MLEPGW[2], MLEPGW[3] ) )

# Fitted Weibull cumulative hazard
fitchpgw <- Vectorize(function(t) chpgw( t, MLEPGW[1], MLEPGW[2], MLEPGW[3] ) )

# Fitted Weibull survival
fitspgw <- Vectorize(function(t) exp(-chpgw( t, MLEPGW[1], MLEPGW[2], MLEPGW[3] ) ))

# AIC
AICPGW <- 2*OPTPGW$OPT$objective + 2*length(OPTPGW$OPT$par)

# BIC
BICPGW <- 2*OPTPGW$OPT$objective + length(OPTPGW$OPT$par)*log(length(survtimes))



## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# MLE Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Harmonic Oscillator model for the hazard function: Solver solution
#--------------------------------------------------------------------------------------------------

# Initial conditions

# Survival at time 0
S0 <- 1
# Survival at 1 month
St <- 0.999
# Survival at 2 months
Stt <- 0.998
# dt = 1 month
dt <- 1/12

# Approximation of the initial conditions
Sp <- -(St-S0)/dt
Spp <- (Stt - 2*St + S0)/(dt^2)

h00 <- Sp/St
r00 <- h00^2 - Spp/St

c(h00,r00)

# Initial point
initHO <- c(0.1,0,0)

# Optimisation step
OPTHO <- optim(initHO, log_likHO, control = list(maxit = 1000))
OPTHO

MLEHO <- exp(OPTHO$par)

# AIC
AICHO <- 2*OPTHO$value + 2*length(OPTHO$par)

# BIC
BICHO <- 2*OPTHO$value + length(OPTHO$par)*log(length(survtimes))


## -----------------------------------------------------------------------------------------------------------------------------------
# AIC comparison
AICs <- c(AICW, AICPGW, AICHO)

AICs

# Best model (Harmonic oscillator)
which.min(AICs)

# BIC comparison
BICs <- c(BICW, BICPGW, BICHO)

BICs

# Best model (Harmonic oscillator)
which.min(BICs)

# Fitted hazard, cumulative hazard, and survival
fithHO <- Vectorize(function(t) hHO(t,exp(OPTHO$par[1]), exp(OPTHO$par[2]), exp(OPTHO$par[3]), h00, r00))
fitchHO <- Vectorize(function(t) chHO(t,exp(OPTHO$par[1]), exp(OPTHO$par[2]), exp(OPTHO$par[3]), h00, r00))
fitsHO <- Vectorize(function(t) exp(-chHO(t,exp(OPTHO$par[1]), exp(OPTHO$par[2]), exp(OPTHO$par[3]), h00, r00)))

# Comparison: hazard functions
curve(fithHO, 0, max(survtimes), xlim= c(0,max(survtimes)), ylim = c(0,0.125), lwd = 2,
     xlab = "Time", ylab = "Hazard", main = "Harmonic Oscillator", cex.axis = 1.5, cex.lab = 1.5)
curve(fithw, 0, max(survtimes), lwd= 2, lty = 2, col = "gray", add = TRUE)
curve(fithpgw, 0, max(survtimes), lwd= 2, lty = 3, col = "gray", add = TRUE)
legend("topleft", legend = c("Harmonic Oscillator","Weibull", "PGW"), lty = c(1,2,3), 
       lwd = c(2,2,2), col = c("black","gray","gray"))

# Comparison: cumulative hazard functions
curve(fitchHO, 0, max(survtimes), xlim = c(0,max(survtimes)), ylim = c(0,2.25),  lwd = 2,
     xlab = "Time", ylab = "Cumulative Hazard", main = "Harmonic Oscillator", cex.axis = 1.5, cex.lab = 1.5)
curve(fitchw, 0, max(survtimes), lwd= 2, lty = 2, col = "gray", add = TRUE)
curve(fitchpgw, 0, max(survtimes), lwd= 2, lty = 3, col = "gray", add = TRUE)
legend("topleft", legend = c("Harmonic Oscillator","Weibull","PGW"), lty = c(1,2,3), 
       lwd = c(2,2,2), col = c("black","gray","gray"))

# Comparison: survival functions
curve(fitsHO, 0, max(survtimes), xlim = c(0,max(survtimes)), ylim = c(0,1), lwd = 2,
     xlab = "Time", ylab = "Survival", main = "Harmonic Oscillator", cex.axis = 1.5, cex.lab = 1.5)
curve(fitsw, 0, max(survtimes), lwd= 2, lty = 2, col = "gray", add = TRUE)
curve(fitspgw, 0, max(survtimes), lwd= 2, lty = 3, col = "gray", add = TRUE)
points(km$time, km$surv, type = "l", col = "green", lwd = 2, lty = 1)
legend("topright", legend = c("Harmonic Oscillator","Weibull","PGW","KM"), lty = c(1,2,3,1), 
       lwd = c(2,2,2,2), col = c("black","gray","gray","green"))



## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(1,2))
p_sigmaW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_sigmaW,0,15, n = 1000, xlab = ~sigma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_nuW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_nuW,0,15, n = 1000, xlab = ~nu, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)


par(mfrow = c(1,1))

#--------------------------------------------------------------------------------------------------
# Weibull model for the hazard function
#--------------------------------------------------------------------------------------------------

# Support
SupportW <- function(x) {   TRUE }

# Random initial points
X0W <- function(x) { OPTW$OPT$par + runif(2,-0.01,0.01) }

# twalk for analytic solution
set.seed(1234)
infoW <- Runtwalk( dim=2,  Tr=110000,  Obj=log_postW, Supp=SupportW, 
                   x0=X0W(), xp0=X0W(), PlotLogPost = FALSE) 


# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries
summW <- apply(exp(infoW$output[ind,]),2,summary)
colnames(summW) <- c("sigma","nu")
kable(summW, digits = 3)

# KDEs
sigmapW <- exp(infoW$output[,1][ind])
nupW <- exp(infoW$output[,2][ind])

plot(density(sigmapW), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_sigmaW,13,17, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(nupW), main = "", xlab = expression(nu), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_nuW,1.1,1.5, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(1,3))
p_sigmaPGW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_sigmaPGW,0,15, n = 1000, xlab = ~sigma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_nuPGW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_nuPGW,0,15, n = 1000, xlab = ~nu, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_gammaPGW <- Vectorize(function(t) dgamma(t, shape = 2, scale = 2))
curve(p_gammaPGW,0,15, n = 1000, xlab = ~gamma, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

par(mfrow = c(1,1))

#--------------------------------------------------------------------------------------------------
# Weibull model for the hazard function
#--------------------------------------------------------------------------------------------------

# Support
SupportPGW <- function(x) {   TRUE }

# Random initial points
X0PGW <- function(x) { OPTPGW$OPT$par + runif(3,-0.01,0.01) }


# twalk for analytic solution
set.seed(1234)
infoPGW <- Runtwalk( dim=3,  Tr=110000,  Obj=log_postPGW, Supp=SupportPGW, 
                     x0=X0PGW(), xp0=X0PGW(), PlotLogPost = FALSE) 


# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries
summPGW <- apply(exp(infoPGW$output[ind,]),2,summary)
colnames(summPGW) <- c("sigma","nu","gamma")
kable(summPGW, digits = 3)

# KDEs
sigmapPGW <- exp(infoPGW$output[,1][ind])
nupPGW <- exp(infoPGW$output[,2][ind])
gammapPGW <- exp(infoPGW$output[,3][ind])

plot(density(sigmapPGW), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_sigmaPGW,2,5, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(nupPGW), main = "", xlab = expression(nu), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_nuPGW,1.5,2.7, n = 1000, xlab = ~kappa, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(gammapPGW), main = "", xlab = expression(gamma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_gammaPGW,0,15, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Priors
#--------------------------------------------------------------------------------------------------

par(mfrow = c(2,2))
p_eta <- Vectorize(function(t) dgamma(t, shape = 4, scale = 0.25))
curve(p_eta,0,2, n = 1000, xlab = ~eta, ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_w0 <- Vectorize(function(t) dgamma(t, shape = 10, scale = 0.1))
curve(p_w0,0,10, n = 1000, xlab = expression(w_0), ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)

p_hb <- Vectorize(function(t) dgamma(t, shape = 0.1, scale = 1))
curve(p_hb,0,0.1, n = 1000, xlab = expression(h_b), ylab = "Prior Density",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2)


par(mfrow = c(1,1))


#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE model for the hazard function: Solver solution
#--------------------------------------------------------------------------------------------------


n.batch <- 1100
batch.length <- 100

lp <- function(par) -log_postHO(par)

inits <- OPTHO$par


## ----include=FALSE------------------------------------------------------------------------------------------------------------------
set.seed(1234)
infoHO<- adaptMetropGibbs(ltd=lp, starting=inits, accept.rate=0.44, batch=n.batch, 
                           batch.length=batch.length, report=100, verbose=FALSE)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------
## set.seed(1234)
## infoHO <- adaptMetropGibbs(ltd=lp, starting=inits, accept.rate=0.44, batch=n.batch,
##                            batch.length=batch.length, report=100, verbose=FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------
chainHO <- infoHO$p.theta.samples[,1:3]

# Burning and thinning the chain
burn <- 1e4
thin <- 100
NS <- n.batch*batch.length
ind <- seq(burn,NS,thin)



# Summaries
summHR <- apply(exp(chainHO[ind,1:3]),2,summary)
colnames(summHR) <- c("eta","w_0","h_b")
kable(summHR, digits = 3)

# KDEs
etap <- exp(chainHO[,1][ind])
w0p <- exp(chainHO[,2][ind])
hbp <- exp(chainHO[,3][ind])


plot(density(etap), main = "", xlab = expression(eta), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_eta,0,1, n = 1000, xlab = ~lambda, ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)


plot(density(w0p), main = "", xlab = expression(w_0), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_w0,0.5,2.5, n = 1000, xlab = expression(w_0), ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)

plot(density(hbp), main = "", xlab = expression(h_b), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
curve(p_hb,0,0.1, n = 1000, xlab = expression(h_b), ylab = "Prior Density", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE)



## -----------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Predictive Hazard functions
#---------------------------------------------------------------------------------------

# Predictive Weibull hazard
predhW <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)){
    num[i] <- exp(-chweibull( t, sigmapW[i], nupW[i]))*hweibull( t, sigmapW[i], nupW[i])
    den[i] <- exp(-chweibull( t, sigmapW[i], nupW[i]))
  }
  return(mean(num)/mean(den))
})

# Predictive PGW
predhPGW <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)) num[i] <- exp(-chpgw( t, sigmapPGW[i], nupPGW[i], gammapPGW[i]))*
      hpgw( t, sigmapPGW[i], nupPGW[i], gammapPGW[i])
  for(i in 1:length(ind)) den[i] <- exp(-chpgw( t, sigmapPGW[i], nupPGW[i], gammapPGW[i]))
  return(mean(num)/mean(den))
})



# Creating the credible envelopes
tvec <- seq(0,20,by = 0.01)
ntvec <- length(tvec)

# Weibull
hCIW <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIW[j,k ] <- hweibull( tvec[k],sigmapW[j], nupW[j]) 
  }
} 

hW <-  predhW(tvec)


hCIWL <- apply(hCIW, 2, ql)
hCIWU <- apply(hCIW, 2, qu)

# PGW
hCIPGW <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIPGW[j,k ] <- hpgw( tvec[k],sigmapPGW[j], nupPGW[j], gammapPGW[j]) 
  }
} 

hPGWs <-  predhPGW(tvec)

hCIPGWL <- apply(hCIPGW, 2, ql)
hCIPGWU <- apply(hCIPGW, 2, qu)

# Harmonic Oscillator

Ap <- vector()
phip <- vector()

for(i in 1:length(ind)) {
  newpar <- rep2ap(etap[i], w0p[i], hbp[i], h00, r00)
  
  Ap[i] <- as.numeric(newpar$A)
  phip[i] <- as.numeric(newpar$phi)
}


hist(Ap)
hist(phip)


# Predictive HO hazard
predhHO <- Vectorize(function(t){
  num <- den <- temp <- vector()
  for(i in 1:length(ind)){
    chi <- chHO( t, etap[i], w0p[i], hbp[i], h00, r00 )
    num[i] <- exp(-chi)*hHO( t, etap[i], w0p[i], hbp[i], h00, r00 )
    den[i] <- exp(-chi)
  }
  return(mean(num)/mean(den))
})


# Creating the credible envelopes
tvec <- seq(0,20,by = 0.01)
ntvec <- length(tvec)

hHOs <-  predhHO(tvec)

hCIHO<- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    hCIHO[j,k ] <- hHO( tvec[k],etap[j], w0p[j], hbp[j], h00, r00) 
  }
} 

hCIHOL <- apply(hCIHO, 2, ql)
hCIHOU <- apply(hCIHO, 2, qu)

# Splines model
fit <- bshazard(Surv(df$time, df$status) ~ 1, data = df, nbin = 100, degree = 3, verbose = FALSE)



# Plots
curve(predhHO, 0, 20, n = 1000, xlab = "Time", ylab = "Predictive Hazard", 
      cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.08))


plot(tvec,  hHOs, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.1))
points(tvec,  hCIHOL, col = "gray", type = "l")
points(tvec,  hCIHOU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCIHOL[order(tvec)], rev(hCIHOU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  hHOs, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1)


## -----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------- 
# Predictive Survival functions 
#--------------------------------------------------------------------------------------- 

# Predictive Weibull survival 
predsW <- Vectorize(function(t){ 
  temp <- vector() 
  for(i in 1:length(ind)) temp[i] <- exp(-chweibull(t,sigmapW[i], nupW[i]) )  
  return(mean(temp)) 
}) 

# Predictive PGW survival 
predsPGW <- Vectorize(function(t){ 
  temp <- vector() 
  for(i in 1:length(ind)) temp[i] <- exp(-chpgw( t,sigmapPGW[i], nupPGW[i], gammapPGW[i]) )  
  return(mean(temp)) 
}) 






## ----message=FALSE------------------------------------------------------------------------------------------------------------------
# fitting the estimator
fit <- bshazard(Surv(df$time, df$status) ~ 1, data = df, nbin = 100, degree = 3, verbose = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------
# Hazard 


plot(tvec,  hHOs, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.1))
points(tvec,  hCIHOL, col = "gray", type = "l")
points(tvec,  hCIHOU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(hCIHOL[order(tvec)], rev(hCIHOU[order(tvec)])),
        col = "gray", border = NA)
points(fit$time, fit$hazard, type='l', lwd = 2, lty = 2, ylim = c(0,0.25), col = "darkgray")
lines(fit$time, fit$lower.ci, lty = 2, lwd = 1, col = "gray")
lines(fit$time, fit$upper.ci, lty = 2, lwd = 1, col = "gray")
points(tvec,  hHOs, type = "l", xlab = "Time", ylab = "Predictive Hazard", 
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1)
legend("bottomright", legend = c("HO", "BS"), lty = c(1,2),  
       lwd = c(2,2), col = c("black","darkgray"))

