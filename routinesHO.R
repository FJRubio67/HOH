# Quantile functions
ql <- function(p){
  quantile(p, 0.025)
}
qu <- function(p){
  quantile(p, 0.975)
}

######################################################################################
######################################################################################
######################################################################################
# Weibull
######################################################################################
######################################################################################
######################################################################################


#--------------------------------------------------------------------------------------------------
# Weibull -log-posterior function
#--------------------------------------------------------------------------------------------------

log_postW <- function(par) {
  sigma <- exp(par[1])
  nu <- exp(par[2])
  
  
  # Terms in the log log likelihood function
  ll_haz <- sum(hweibull(t_obs, sigma, nu, log = TRUE))
  
  ll_chaz <- sum(chweibull(survtimes, sigma, nu))
  
  log_lik <- -ll_haz + ll_chaz
  
  # Log prior
  
  log_prior <- -dgamma(sigma,
                       shape = 2,
                       scale = 2,
                       log = TRUE) -
    dgamma(nu,
           shape = 2,
           scale = 2,
           log = TRUE)
  
  # Log Jacobian
  
  log_jacobian <- -par[1] - par[2]
  
  # log posterior
  
  log_post0 <- log_lik + log_prior + log_jacobian
  
  return(as.numeric(log_post0))
}

######################################################################################
######################################################################################
######################################################################################
# Power Generalised Weibull
######################################################################################
######################################################################################
######################################################################################


#--------------------------------------------------------------------------------------------------
# Power Generalised Weibull -log-posterior function
#--------------------------------------------------------------------------------------------------

log_postPGW <- function(par) {
  sigma <- exp(par[1])
  nu <- exp(par[2])
  gamma <- exp(par[3])
  
  # Terms in the log log likelihood function
  ll_haz <- sum(hpgw(t_obs, sigma, nu, gamma, log = TRUE))
  
  ll_chaz <- sum(chpgw(survtimes, sigma, nu, gamma))
  
  log_lik <- -ll_haz + ll_chaz
  
  # Log prior
  
  log_prior <- -dgamma(sigma,
                       shape = 2,
                       scale = 2,
                       log = TRUE) -
    dgamma(nu,
           shape = 2,
           scale = 2,
           log = TRUE) -
    dgamma(gamma,
           shape = 2,
           scale = 2,
           log = TRUE)
  
  # Log Jacobian
  
  log_jacobian <- -par[1] - par[2] - par[3]
  
  # log posterior
  
  log_post0 <- log_lik + log_prior + log_jacobian
  
  return(as.numeric(log_post0))
}


######################################################################################
######################################################################################
######################################################################################
# Harmonic Oscillator
######################################################################################
######################################################################################
######################################################################################


#-----------------------------------------------------------------------------
# Reparameterisations
#-----------------------------------------------------------------------------


# reparameterisation to A and phi
rep2ap <- function(eta, w0, hb, h0, r0) {
  test <- r0 + w0 * eta * (h0 - hb)
  testr <- h0 - hb
  w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
  
  if (testr == 0) {
    phi <- 0
    A <- r0 / w1
  }
  
  if (test == 0) {
    phi <- sign(h0 - hb) * pi / 2
    A <- h0 - hb
  }
  
  if (test != 0 & testr != 0) {
    phi <- atan((w1 * (h0 - hb)) / (r0 + w0 * eta * (h0 - hb)))
    A <- (h0 - hb) / sin(phi)
  }
  return(list(A = A, phi = phi))
}


#-----------------------------------------------------------------------------
# Harmonic Oscillator ODE Hazard Function: Analytic solution
#-----------------------------------------------------------------------------
# t: time (positive)
# eta: shape parameter (positive)
# w0 : shape parameter (positive)
# hb: shift parameter (positive)
# h0 : hazard initial condition (positive)
# r0 : hazard derivative initial condition (real)

hHO <- function(t, eta, w0, hb, h0, r0) {
  # Notation
  w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
  
  # Under-damped
  if (eta < 1) {
    rep = rep2ap(eta, w0, hb, h0, r0)
    A <- rep$A
    phi <- rep$phi
    haz <-  hb + A * exp(-w0 * eta * t) * sin(w1 * t + phi)
    }
  
  # Over-damped
  if (eta > 1) {
    mu = w1 / (w0 * eta)
    a = (h0 - hb) / mu + r0 / w1
    
haz <- hb + exp(-w0 * eta * t) * (0.5 * (h0 - hb + a) * exp(w1 * t) + 
                                 0.5 * (h0 - hb - a) * exp(-w1 * t)) 
  }
  
  # Critically damped
  if (eta == 1) {
    haz <- hb + (h0 - hb + t * (r0 + w0 * (h0 - hb))) * exp(-w0 * tv)
  }
  

    return(haz)
  
}


#-----------------------------------------------------------------------------
# Harmonic Oscillator ODE Cumulative Hazard Function: Analytic solution
#-----------------------------------------------------------------------------
# t: time (positive)
# eta: shape parameter (positive)
# w0 : shape parameter (positive)
# hb: shift parameter (positive)
# h0 : hazard initial condition (positive)
# r0 : hazard derivative initial condition (real)

chHO <- function(t, eta, w0, hb, h0, r0) {
  # Notation
  w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
  mu = w1 / (w0 * eta)
  
  # Under-damped
  if (eta < 1) {
    rep = rep2ap(eta, w0, hb, h0, r0)
    A <- rep$A
    phi <- rep$phi
    chaz <- hb * t +
      (A / (w0 * eta)) * (sin(phi) + mu * cos(phi) -
                            exp(-w0 * eta * t) * (sin(w1 * t + phi) + mu * cos(w1 * t +
                                                                                 phi))) / (mu ^ 2 + 1)
  }
  
  # Over-damped
  if (eta > 1) {
    a = (h0 - hb) / mu + r0 / w1
    
    chaz <- hb * t + (0.5 * (h0 - hb + a) * (exp((w1 - w0 * eta) * t) -
                                               1) / (w1 - w0 * eta) -
                        0.5 * (h0 - hb - a) * (exp((-w1 - w0 * eta) * t) -
                                                 1) / (w1 + w0 * eta))
  }
  
  # Critically-damped
  if (eta == 1) {
    chaz <- hb * t + (h0 - hb) * (1 - exp(-w0 * t)) / w0 +
      (r0 + w0 * (h0 - hb)) * (exp(w0 * t) - w0 * t - 1) * exp(-w0 *
                                                                 t) / (w0 ^ 2)
  }
  return(chaz)
}

#-----------------------------------------------------------------------------
# Shifted Harmonic Oscillator ODE Probability Density Function: Analytic solution
#-----------------------------------------------------------------------------
# t: time (positive)
# tau: shape parameter
# w2 : shape parameter
# hb: shift parameter
# A: Amplitude
# phi: Phase
dHO  <- function(t, eta, w0, hb, h0, r0) {
  den <-  hHO(t, eta, w0, hb, h0, r0) *
    exp(-chHO(t, eta, w0, hb, h0, r0))
  return(den)
}

#-----------------------------------------------------------------------------
# Shifted Harmonic Oscillator ODE random number generation: Analytic solution
#-----------------------------------------------------------------------------
# n: number of simulations
# tau: shape parameter
# w2 : shape parameter
# hb: shift parameter
# A: Amplitude
# phi: Phase
# interval: interval to find the roots

rHO <- function(n, interval, eta, w0, hb, h0, r0) {
  u <- runif(n)
  out <- vector()
  
  for (i in 1:n) {
    tempf <- Vectorize(function(t)
      chHO(t, eta, w0, hb, h0, r0) + log(u[i]))
    out[i] <- uniroot(tempf, interval = interval)$root
  }
  
  return(as.vector(out))
}




######################################################################################################
# Function to plot the hazard function
######################################################################################################
# par: parameters
# xl: lower limit for plot - x axis
# xu: upper limit for plot - y axis
# yl: lower limit for plot - y axis
# yu: upper limit for plot - y axis
# param: parametrization. "original" = c(tau, w2, hb, h0, r0), repar = c(tau, w2, hb, A, phi)
# n: number of points

plotSHOhaz <- function(par,
                       xl,
                       xu,
                       yl,
                       yu,
                       add.plot = FALSE,
                       col = "black",
                       lty = 1,
                       lwd = 2,
                       param = "original",
                       n = 100) {
  if (param == "original") {
    haz <- Vectorize(function(t)
      hHO(t, par[1], par[2], par[3], par[4], par[5]))
  }
  if (param == "repar") {
    reparam <- rep2ap(par[1], par[2], par[3], par[4], par[5])
    A <- reparam$A
    phi <- reparam$phi
    haz <- Vectorize(function(t)
      hHO(t, par[1], par[2], par[3], A, phi))
  }
  curve(
    haz,
    xl,
    xu,
    ylim = c(yl, yu),
    col = col,
    cex.axis = 1.25,
    cex.lab = 1.25,
    ylab = "h(t)",
    xlab = "t",
    lty = lty,
    lwd = lwd,
    add = add.plot,
    n = n
  )
}




######################################################################################
######################################################################################
######################################################################################
# Shifted Harmonic Oscillator
######################################################################################
######################################################################################
######################################################################################


# Support function on the log scale
# SupportHO <- function(par) {
#   # Parameters
#   eta <- exp(par[1])
#   w0 <- exp(par[2])
#   hb <- exp(par[3])
#   
#   h0 <- exp(par[4])
#   r0 <- par[5]
#   
#   # Additional required parameters
#   w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
#   mu <- w1 / (w0 * eta)
#   phi <- rep2ap(eta, w0, hb, h0, r0)$phi
#   
#   # HO Hazard function with parameters given by par
#   temph <- Vectorize(function(t)
#     hHO(t, eta, w0, hb, h0, r0))
#   
#   if (eta < 1) {
#     ts <- (atan(mu) - phi + seq(0, 1, by = 1) * pi) / w1
#     hs <- temph(ts)
#     out <- ifelse(any(hs <= 0), FALSE, TRUE)
#   }
#   
#   if (eta > 1) {
#     if (r0 >= 0)
#       out <- TRUE
#     if (r0 < 0) {
#       test <- (h0 - hb - a) / (h0 - hb + a)
#       if (test <= 0)
#         out = FALSE
#       if (test > 0) {
#         t_hat <- 0.5 * w0 * eta * log(test) / mu
#         
#         hs <- temph(c(0.5 * t_hat, t_hat))
#         
#         out <- ifelse(any(hs <= 0), FALSE, TRUE)
#       }
#     }
#     
#   }
#   if (eta == 1){
#     out <- FALSE
#   }
#   
#   return(out)
#   
# }


# Support function on the log scale
# Fixed initial conditions
SupportHOFICS <- function(par) {
  # Parameters
  eta <- exp(par[1])
  w0 <- exp(par[2])
  hb <- exp(par[3])
  
  h0 <- h00
  r0 <- r00
  
  # Additional required parameters
  w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
  mu <- w1 / (w0 * eta)
  phi <- rep2ap(eta, w0, hb, h0, r0)$phi
  
  # HO Hazard function with parameters given by par
  temph <- Vectorize(function(t)
    hHO(t, eta, w0, hb, h0, r0))
  
  if (eta < 1) {
    ts <- (atan(mu) - phi + seq(0, 1, by = 1) * pi) / w1
    hs <- temph(ts)
    out <- ifelse(any(hs <= 0), FALSE, TRUE)
  }
  
  if (eta > 1) {
    mu <- w1/(w0*eta)
    a <- (h0-hb)/mu + r0/w1
    test <- ((h0 - hb - a) * (w1+w0*eta)) / ((h0 - hb + a)*(w1-w0*eta))
    
    if (r0 >= 0){
      if (test <= 0){
        out = FALSE
      }
      if (test > 0) {
        t_hat <- 0.5 * log(test) / w1
        hs <- temph(seq(1,10)*t_hat/10)
        out <- ifelse(any(hs <= 0), FALSE, TRUE)
      }
      
    }
    if (r0 < 0) {
      if (test <= 0){
        out = FALSE
      }
      
      if (test > 0) {
        t_hat <- 0.5 * log(test) / w1
        hs <- temph(seq(1,10)*t_hat/10)
        out <- ifelse(any(hs <= 0), FALSE, TRUE)
      }
    }
  }
  
  
  if (eta == 1){
    out <- FALSE
  }
  
  return(out)
  
}
# 
# #--------------------------------------------------------------------------------------------------
# # Shifted Harmonic Oscillator ODE -log-likelihood function: Analytic solution
# #--------------------------------------------------------------------------------------------------
# 
# log_likHO <- function(par) {
#   tau <- exp(par[1])
#   w2 <- exp(par[2])
#   hb <- exp(par[3])
#   
#   h0 <- exp(par[4])
#   r0 <- par[5]
#   
#   newpar <- reo2ap(tau, w2, hb, h0, r0)
#   
#   A <- as.numeric(newpar$A)
#   phi <- as.numeric(newpar$phi)
#   
#   cond <- SupportHO(par)
#   
#   if (!cond)
#     log_lik <- Inf
#   
#   if (cond) {
#     # Terms in the log log likelihood function
#     ll_haz <- sum(log(hHO(t_obs, tau, w2, hb, A, phi)))
#     
#     ll_chaz <- sum(chHO(survtimes, tau, w2, hb, A, phi))
#     
#     log_lik <- -ll_haz + ll_chaz
#   }
#   
#   return(log_lik)
# }
# 
# 
# 
# #--------------------------------------------------------------------------------------------------
# # Shifted Harmonic Oscillator ODE -log-posterior function: Analytic solution
# #--------------------------------------------------------------------------------------------------
# 
# log_postHO <- function(par) {
#   tau <- exp(par[1])
#   w2 <- exp(par[2])
#   hb <- exp(par[3])
#   
#   h0 <- exp(par[4])
#   r0 <- par[5]
#   
#   newpar <- rep2ap(tau, w2, hb, h0, r0)
#   
#   A <- as.numeric(newpar$A)
#   phi <- as.numeric(newpar$phi)
#   
#   # Terms in the log log likelihood function
#   ll_haz <- sum(log(hHO(t_obs, tau, w2, hb, A, phi)))
#   
#   ll_chaz <- sum(chHO(survtimes, tau, w2, hb, A, phi))
#   
#   log_lik <- -ll_haz + ll_chaz
#   
#   # Log prior
#   
#   
#   log_prior <- -dgamma(tau,
#                        shape = 2,
#                        scale = 2,
#                        log = TRUE) -
#     dgamma(w2,
#            shape = 2,
#            scale = 2,
#            log = TRUE) -
#     dgamma(hb,
#            shape = 2,
#            scale = 2,
#            log = TRUE) -
#     dgamma(h0,
#            shape = 2,
#            scale = 2,
#            log = TRUE) -
#     dnorm(r0,
#           mean = 0,
#           sd = 100,
#           log = TRUE)
#   
#   # Log-Jacobian
#   
#   log_jacobian <- -log(tau) - log(w2) - log(hb) - log(h0)
#   
#   # log posterior
#   
#   log_post0 <- log_lik + log_prior + log_jacobian
#   
#   return(as.numeric(log_post0))
# }



#--------------------------------------------------------------------------------------------------
# Shifted Harmonic Oscillator ODE -log-likelihood function: Analytic solution
# Fixed initial conditions
#--------------------------------------------------------------------------------------------------

log_likHOFICS <- function(par) {
  eta <- exp(par[1])
  w0 <- exp(par[2])
  hb <- exp(par[3])
  
  
  cond <- SupportHOFICS(par)
  
  if (!cond)
    log_lik <- Inf
  
  if (cond) {
    # Terms in the log log likelihood function
    ll_haz <- sum(log(hHO(t_obs, eta, w0, hb, h00, r00)))
    
    ll_chaz <- sum(chHO(survtimes, eta, w0, hb, h00, r00))
    
    log_lik <- -ll_haz + ll_chaz
  }
  
  return(log_lik)
}



#--------------------------------------------------------------------------------------------------
# Shifted Harmonic Oscillator ODE -log-posterior function: Analytic solution
# Fixed initial conditions
#--------------------------------------------------------------------------------------------------

log_postHOFICS <- function(par) {
  eta <- exp(par[1])
  w0 <- exp(par[2])
  hb <- exp(par[3])
  
  
  cond <- SupportHOFICS(par)
  
  if (!cond)
    log_lik <- Inf
  
  if (cond) {
    # Terms in the log log likelihood function
    ll_haz <- sum(log(hHO(t_obs, eta, w0, hb, h00, r00)))
    
    ll_chaz <- sum(chHO(survtimes, eta, w0, hb, h00, r00))
    
    log_lik <- -ll_haz + ll_chaz
  }
  
  # Log prior
  
  
  log_prior <- -dgamma(eta,
                       shape = 0.5,
                       scale = 1,
                       log = TRUE) -
    dgamma(w0,
           shape = 2,
           scale = 1,
           log = TRUE) -
    dgamma(hb,
           shape = 0.1,
           scale = 0.1,
           log = TRUE)
  
  # Log-Jacobian
  
  log_jacobian <- -log(eta) - log(w0) - log(hb)
  
  # log posterior
  
 log_post0 <- log_lik + log_prior + log_jacobian
   
  return(as.numeric(log_post0))
}
