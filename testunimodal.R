eta = 2.4
w0 = 1.4
hb = 2.5
h0 = 2.8
r0 = -1.8
w1 <- w0 * sqrt(abs(eta ^ 2 - 1))
mu <- w1/(w0*eta)
a <- (h0-hb)/mu + r0/w1


((h0 - hb - a) * (w1+w0*eta)) / ((h0 - hb + a)*(w1-w0*eta))
