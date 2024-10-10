source("routinesHO.R")


###############################################################################
# Parameter values
###############################################################################

eta <- 0.5
w0 <- 1.2
hb <- 0.5
h0 <- 0.1
r0 <- 0.5
w1 <- w0*sqrt(abs(eta^2-1))
A <- rep2ap(eta, w0, hb, h0, r0)$A
phi <- rep2ap(eta, w0, hb, h0, r0)$phi
mu <- w1/(w0*eta)


###############################################################################
# Under-damped case (eta < 1)
###############################################################################


K1 <- (abs(A)/(w0*eta))*(1+mu)/(1+mu^2)
K2 <- w0*eta
K = exp( -(A/(w0*eta))*( sin(phi) + mu*cos(phi) )/(mu^2+1) )

bound <- Vectorize( function(t) K*exp(-hb*t + K1*exp(-K2*t)) )

survf <- Vectorize(function(t) exp(-chHO(t, eta, w0, hb, h0, r0))   )

curve(bound,0,5, lwd = 2, n = 1000)
curve(survf, 0, 5, lwd = 2, lty = 2, col = "red", add = TRUE, n = 1000)


curve(bound,20,25, lwd = 2)
curve(survf, 20, 25, lwd = 1, lty = 2, col = "red", add = TRUE)

curve(bound,20,21, lwd = 2)
curve(survf, 20, 21, lwd = 1, lty = 2, col = "red", add = TRUE)




###############################################################################
# Parameter values
###############################################################################

eta <- 1.5
w0 <- 1.2
hb <- 0.5
h0 <- 0.1
r0 <- 0.5
w1 <- w0*sqrt(abs(eta^2-1))
A <- rep2ap(eta, w0, hb, h0, r0)$A
phi <- rep2ap(eta, w0, hb, h0, r0)$phi
mu <- w1/(w0*eta)
a <- (h0-hb)/mu + r0/w1

###############################################################################
# Over-damped case (eta > 1)
###############################################################################


K1 <- max( c((hb-h0-a)/(w1-w0*eta), (h0-hb-a)/(w1+w0*eta)  ))
K2 <- w0*eta - w1
K = exp( -0.5*(hb-h0-a)/(w1-w0*eta) - 0.5*(h0-hb-a)/(w1+w0*eta)  )

bound <- Vectorize( function(t) K*exp(-hb*t + K1*exp(-K2*t)) )

survf <- Vectorize(function(t) exp(-chHO(t, eta, w0, hb, h0, r0))   )

curve(bound,0,5, lwd = 2, n = 1000, ylim=c(0,2))
curve(survf, 0, 5, lwd = 2, lty = 2, col = "red", add = TRUE, n = 1000)


curve(bound,20,25, lwd = 2)
curve(survf, 20, 25, lwd = 1, lty = 2, col = "red", add = TRUE)

curve(bound,20,21, lwd = 2)
curve(survf, 20, 21, lwd = 1, lty = 2, col = "red", add = TRUE)


###############################################################################
# Parameter values
###############################################################################

eta <- 1
w0 <- 1.2
hb <- 0.5
h0 <- 0.1
r0 <- 0.5
w1 <- w0*sqrt(abs(eta^2-1))
A <- rep2ap(eta, w0, hb, h0, r0)$A
phi <- rep2ap(eta, w0, hb, h0, r0)$phi
mu <- w1/(w0*eta)
a <- (h0-hb)/mu + r0/w1


###############################################################################
# Critically-damped case (eta = 1)
###############################################################################


K1 <- max( c((h0-hb)/w0, (r0 + w0*(h0-hb))/(w0^2)  ))
K2 <- w0
K3 <- (r0 + w0*(h0-hb))/w0
K = exp( -(h0-hb)/w0 -(r0 + w0*(h0-hb))/(w0^2) )

bound <- Vectorize( function(t) K*exp(-hb*t + K1*exp(-K2*t) + K3*t*exp(-K2*t)  ) )

survf <- Vectorize(function(t) exp(-chHO(t, eta, w0, hb, h0, r0))   )

curve(bound,0,5, lwd = 2, n = 1000, ylim=c(0,2))
curve(survf, 0, 5, lwd = 2, lty = 2, col = "red", add = TRUE, n = 1000)


curve(bound,20,25, lwd = 2)
curve(survf, 20, 25, lwd = 1, lty = 2, col = "red", add = TRUE)

curve(bound,20,21, lwd = 2)
curve(survf, 20, 21, lwd = 1, lty = 2, col = "red", add = TRUE)

