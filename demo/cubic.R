## Numerical illustration: simulation with a cubic structural function

library(npiv)
library(MASS)

## Simulate the data
n <- 10000
cov.ux <- 0.5
var.u <- 0.1
mu <- c(1,1,0)
Sigma <- matrix(c(1.0,0.85,cov.ux,
                  0.85,1.0,0.0,
                  cov.ux,0.0,1.0),
                3,3,
                byrow=TRUE)
foo <- mvrnorm(n = n,
               mu,
               Sigma)
X <- 2*pnorm(foo[,1],mean=mu[1],sd=sqrt(Sigma[1,1])) -1
W <- 2*pnorm(foo[,2],mean=mu[2],sd=sqrt(Sigma[2,2])) -1
U <- foo[,3]
## Cubic structural function
h0 <- X^3
Y <- h0 + sqrt(var.u)*U
## Create evaluation data for plotting
X.eval <- seq(min(X),max(X),length=100)
h0 <- X.eval^3
d0 <- 3*X.eval^2

## Call the npiv() function with specific arguments
model <- npiv(Y, X, W, X.eval=X.eval)

## Plot the estimated function and uniform confidence bands
plot(model, showdata = TRUE)
## Add true function
lines(X.eval, h0)

## Plot the estimated derivative and uniform confidence bands
plot(model, type = "deriv")
## Add true function
lines(X.eval, d0)
