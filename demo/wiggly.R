library(npiv)

#?npiv

## Simulate data using the mvrnorm() function in the MASS package

library(MASS)

## Simulate data using the mvrnorm() function in the MASS package

n <- 10000

X <- runif(n)
U <- rnorm(n)

## h0 is the instrumental DGP function

h0 <- sin(15*pi*X)*cos(X)
Y <- h0 + U
X.eval <- seq(min(X),max(X),length=100)
h0 <- sin(15*pi*X.eval)*cos(X.eval)

## Create evaluation data and instrumental regression function for the
## endogenous predictor (for plotting with lines as this is sorted)



## Call the npiv() function with specific arguments

model <- npiv(Y,
              X,
              X,
              X.eval=X.eval)

## Create a plot of the instrumental regression function and its
## asymptotic standard error bounds

ylim <- c(min(Y,model$fitted-1.96*model$fitted.asy.se,model$fitted+1.96*model$fitted.asy.se),
          max(Y,model$fitted-1.96*model$fitted.asy.se,model$fitted+1.96*model$fitted.asy.se))

plot(X,Y,cex=0.25,
     col="lightgrey",
     ylim=ylim,
     sub=paste("n = ",format(n,format="d", big.mark=','),sep=""),
     xlab="X",
     ylab="Y")

lines(X.eval,h0,lty=1,col=1,lwd=1)
lines(X.eval,model$fitted,lty=2,col=2,lwd=2)

legend("topleft",c("NP-DGP","NPREG"),col=1:2,lty=1:2,lwd=c(1,2),bty="n")
