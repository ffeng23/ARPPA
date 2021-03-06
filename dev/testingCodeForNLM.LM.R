## function to simulate data 
f <- function(TT, tau, N0, a, f0) {
    expr <- expression(N0*exp(-TT/tau)*(1 + a*cos(f0*TT)))
    eval(expr)
}

## helper function for an analytical gradient 
j <- function(TT, tau, N0, a, f0) {
    expr <- expression(N0*exp(-TT/tau)*(1 + a*cos(f0*TT)))
    c(eval(D(expr, "tau")), eval(D(expr, "N0" )),
      eval(D(expr, "a"  )), eval(D(expr, "f0" )))
}

## values over which to simulate data 
TT <- seq(0, 8, length=501)

## parameter values underlying simulated data  
p <- c(tau = 2.2, N0 = 1000, a = 0.25, f0 = 8)

## get data 
Ndet <- do.call("f", c(list(TT = TT), as.list(p)))
## with noise
N <- Ndet +  rnorm(length(Ndet), mean=Ndet, sd=.01*max(Ndet))

## plot the data to fit
par(mfrow=c(2,1), mar = c(3,5,2,1))  
plot(TT, N, bg = "black", cex = 0.5, main="data")

## define a residual function 
fcn     <- function(p, TT, N, fcall, jcall)
    (N - do.call("fcall", c(list(TT = TT), as.list(p))))

## define analytical expression for the gradient 
fcn.jac <- function(p, TT, N, fcall, jcall) 
    -do.call("jcall", c(list(TT = TT), as.list(p)))

## starting values 
guess <- c(tau = 2.2, N0 = 1500, a = 0.25, f0 = 10)

## to use an analytical expression for the gradient found in fcn.jac
## uncomment jac = fcn.jac
out <- nls.lm(par = guess, fn = fcn, jac = fcn.jac,
              fcall = f, jcall = j,
              TT = TT, N = N, control = nls.lm.control(nprint=1))

## get the fitted values 
N1 <- do.call("f", c(list(TT = TT), out$par))   

## add a blue line representing the fitting values to the plot of data 
lines(TT, N1, col="blue", lwd=2)

## add a plot of the log residual sum of squares as it is made to
## decrease each iteration; note that the RSS at the starting parameter
## values is also stored
plot(1:(out$niter+1), log(out$rsstrace), type="b",
main="log residual sum of squares vs. iteration number",
xlab="iteration", ylab="log residual sum of squares", pch=21,bg=2) 

## get information regarding standard errors
summary(out) 
