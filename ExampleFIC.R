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
library(Rtwalk)
library(knitr)
library(matrixStats)


## ----include=FALSE----------------------------------------------------------------------------------------

#source("/Users/FJRubio/Dropbox/ODESurv/HO/codes/routines/routinesHO.R")
#source("C:/Users/Javier/Dropbox/ODESurv/HO/codes/routines/routinesHO.R")
#source("/Users/javierrubio/Dropbox/ODESurv/HO/codes/routines/routinesHO.R")
source("routinesHO.R")

## ----eval=FALSE-------------------------------------------------------------------------------------------
## source("routines.R")



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





## ---------------------------------------------------------------------------------------------------------
#==================================================================================================
# MLE Analysis
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# Harmonic Oscillator ODE model for the hazard function: Analytic solution
#--------------------------------------------------------------------------------------------------

# Initial conditions

St <- 0.999
Stt <- 0.998
S0 <- 1
dt <- 1/12

Sp <- -(St-S0)/dt

h0s <- Sp/St

h0s

Spp <- (Stt - 2*St + S0)/(dt^2)
h0p <- h0s^2 - Spp/St

h0p

h00 <- h0s
r00 <- h0p


# Initial point
initHO <- c(0.1,0,0)

OPT <- optim(initHO, log_likHOFICS, control = list(maxit = 1000))
OPT

temph <- Vectorize(function(t) hHO(t,exp(OPT$par[1]), exp(OPT$par[2]), exp(OPT$par[3]), h00, r00))
curve(temph,0,20, lwd = 2)

## ---------------------------------------------------------------------------------------------------------
#==================================================================================================
# Bayesian Analysis
#==================================================================================================


#--------------------------------------------------------------------------------------------------
# Harmonic Oscillator ODE model for the hazard function: Analytic solution
#--------------------------------------------------------------------------------------------------


# Random initial points
X0HO <- function(x) { OPT$par + runif(3,-0.01,0.01) }

NMC = 110000

# twalk for analytic solution
set.seed(1234)
infoHO <- Runtwalk( dim=3,  Tr=NMC,  Obj=log_postHOFICS, Supp=SupportHOFICS, 
                    x0=X0HO(), xp0=X0HO(),PlotLogPost = TRUE) 

ind = ind=seq(5000,NMC,100) 


postHO <- cbind(exp(infoHO$output[ind,1:3]))

etap <- as.vector(postHO[,1])
w0p <- as.vector(postHO[,2])
hbp <- as.vector(postHO[,3])

Ap <- vector()
phip <- vector()

for(i in 1:length(ind)) {
  newpar <- rep2ap(etap[i], w0p[i], hbp[i], h00, r00)
  
  Ap[i] <- as.numeric(newpar$A)
  phip[i] <- as.numeric(newpar$phi)
}


hist(etap)
hist(w0p)
hist(hbp)
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
fit <- bshazard(Surv(df$time, df$status) ~ 1, data = df, nbin = 100, degree = 3, verbose = FALSE)

curve(predhHO, 0, 20, n = 1000, xlab = "Time", ylab = "Predictive Hazard", 
      cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1, ylim = c(0,0.08))


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



# Predictive HO survival 
predSHO <- Vectorize(function(t){ 
  temp <- vector() 
  for(i in 1:length(ind)) temp[i] <- exp(-chshoode( t,taup[i], w2p[i], hbp[i], Ap[i], phip[i]) )  
  return(mean(temp)) 
}) 



SCIHO <- matrix(0, ncol = ntvec, nrow = length(ind)) 

for(j in 1:length(ind)){ 
  for(k in 1:ntvec){ 
    SCIHO[j,k ] <- exp(-chshoode( tvec[k],taup[j], w2p[j], hbp[j], Ap[j], phip[j])) 
  } 
}


SpredHO <-  predSHO(tvec) 

SCIHOL <- apply(SCIHO, 2, ql) 
SCIHOU <- apply(SCIHO, 2, qu) 


plot(tvec,  SpredHO, type = "l", ylim = c(0,1), xlab = "Time", ylab = "Predictive Survival",  
       cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 1) 
points(tvec,  SCIHOL, col = "gray", type = "l") 
points(tvec,  SCIHOU, col = "gray", type = "l") 
polygon(c(tvec, rev(tvec)), c(SCIHOL[order(tvec)], rev(SCIHOU[order(tvec)])), 
        col = "gray", border = NA) 
points(tvec,  SpredHO,type = "l", col = "black", lwd = 2, lty =1) 
points(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 2, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
legend("topright", legend = c("HO", "KM"), lty = c(1,2),  
       lwd = c(2,2), col = c("black","black")) 
