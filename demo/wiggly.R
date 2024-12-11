## Numerical illustration: simulation with a wiggly regression function

library(npiv)
library(MASS)

## Simulate the data
n <- 10000
X <- runif(n)
U <- rnorm(n)
## Wiggly regression function
h0 <- sin(15*pi*X)*cos(X)
Y <- h0 + U
## Create evaluation data for plotting
X.eval <- seq(min(X),max(X),length=100)
h0 <- sin(15*pi*X.eval)*cos(X.eval)
d0 <- 15*pi*cos(X.eval)*cos(15*pi*X.eval) - sin(X.eval)*sin(15*pi*X.eval)

## Call the npiv() function with specific arguments
model <- npiv(Y, X, X, X.eval=X.eval)

## Plot the estimated function and uniform confidence bands
plot(model, showdata = TRUE)
## Add true function
lines(X.eval, h0)

## Plot the estimated derivative and uniform confidence bands
plot(model, type = "deriv")
## Add true function
lines(X.eval, d0)
