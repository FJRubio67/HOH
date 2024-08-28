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
                       shape = 0.001,
                       scale = 1/0.001,
                       log = TRUE) -
    dgamma(nu,
           shape = 0.001,
           scale = 1/0.001,
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
                       shape = 0.001,
                       scale = 1/0.001,
                       log = TRUE) -
    dgamma(nu,
           shape = 0.001,
           scale = 1/0.001,
           log = TRUE) -
    dgamma(gamma,
           shape = 0.001,
           scale = 1/0.001,
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

######################################################################################
######################################################################################
######################################################################################
# Harmonic Oscillator ODE
######################################################################################
######################################################################################
######################################################################################


#-----------------------------------------------------------------------------
# Hazard function
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
    haz <- hb + (h0 - hb + t * (r0 + w0 * (h0 - hb))) * exp(-w0 * t)
  }
  
    return(haz)
}


#-----------------------------------------------------------------------------
# Cumulative hazard function
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
# Probability Density Function
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
# Random number generation
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

#-----------------------------------------------------------------------------
# Critical value, test for critical value, and hazard at critical value
#-----------------------------------------------------------------------------

crit <- function(eta,w0,hb,h0,r0){
  w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
  mu <- w1 / (w0 * eta)
  a <- (h0-hb)/mu + r0/w1
 test <- sign((h0 - hb - a)  / ((h0 - hb + a)*(w1-w0*eta)))
 crit <- 0.5*log( ((h0 - hb - a)*(w1+w0*eta))  / ((h0 - hb + a)*(w1-w0*eta)) )/w1
 hcrit <- hHO(crit,eta,w0,hb,h0,r0)
 out = list(test = test, crit = crit, hcrit = hcrit)
 return(out)
}

#-----------------------------------------------------------------------------
# Support function on the log scale
# Fixed initial conditions
#-----------------------------------------------------------------------------
# par = (log(eta), log(w0), log(hb))

SupportHO <- function(par) {
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
    a <- (h0-hb)/mu + r0/w1
    test <- (h0 - hb - a)  / ((h0 - hb + a)*(w1-w0*eta))
    
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


#--------------------------------------------------------------------------------------------------
# Harmonic Oscillator ODE -log-likelihood function
# Fixed initial conditions
#--------------------------------------------------------------------------------------------------

log_likHO <- function(par) {
  eta <- exp(par[1])
  w0 <- exp(par[2])
  hb <- exp(par[3])
  
  
  cond <- SupportHO(par)
  
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
# Harmonic Oscillator ODE -log-posterior function: Analytic solution
# Fixed initial conditions
#--------------------------------------------------------------------------------------------------

log_postHO <- function(par) {
  eta <- exp(par[1])
  w0 <- exp(par[2])
  hb <- exp(par[3])
  
  
  cond <- SupportHO(par)
  
  if (!cond)
    log_lik <- Inf
  
  if (cond) {
    # Terms in the log log likelihood function
    ll_haz <- sum(log(hHO(t_obs, eta, w0, hb, h00, r00)))
    
    ll_chaz <- sum(chHO(survtimes, eta, w0, hb, h00, r00))
    
    log_lik <- -ll_haz + ll_chaz
  }
  

  log_prior <- -dgamma(eta,
                       shape = 0.001,
                       scale = 1/0.001,
                       log = TRUE) -
    dgamma(w0,
           shape = 0.001,
           scale = 1/0.001,
           log = TRUE) -
    dgamma(hb,
           shape = 0.001,
           scale = 1/0.001,
           log = TRUE)
  
  # Log-Jacobian
  
  log_jacobian <- -log(eta) - log(w0) - log(hb)
  
  # log posterior
  
 log_post0 <- log_lik + log_prior + log_jacobian
   
  return(as.numeric(log_post0))
}


