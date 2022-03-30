library(npiv)

#?npiv

## Simulate data using the mvrnorm() function in the MASS package

library(MASS)

## Simulate data using the mvrnorm() function in the MASS package

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

## h0 is the instrumental DGP function

h0 <- sin(pi*X)
Y <- h0 + sqrt(var.u)*U
X.eval <- seq(min(X),max(X),length=100)
h0 <- sin(pi*X.eval)

## Create evaluation data and instrumental regression function for the
## endogenous predictor (for plotting with lines as this is sorted)



## Call the npiv() function with specific arguments

model <- npiv(Y,
              X,
              W,
              X.eval=X.eval)

## Create a plot of the instrumental regression function and its
## asymptotic standard error bounds

ylim <- c(min(Y,model$h-1.96*model$h.asy.se,model$h+1.96*model$h.asy.se),
          max(Y,model$h-1.96*model$h.asy.se,model$h+1.96*model$h.asy.se))

plot(X,Y,cex=0.25,
     col="lightgrey",
     ylim=ylim,
     sub=paste("n = ",format(n,format="d", big.mark=','),sep=""),
     xlab="X",
     ylab="Y")

lines(X.eval,h0,lty=1,col=1,lwd=1)
lines(X.eval,model$h,lty=2,col=2,lwd=2)

legend("topleft",c("IV-DGP","NPIV"),col=1:2,lty=1:2,lwd=c(1,2),bty="n")
