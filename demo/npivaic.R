## Model selection via Hurvich, Siminoff and Tsai's AIC.c
## criterion. Note we set the degree and number of knots/segments to
## be equal for both the Psi.x and B.w bases. You could, naturally,
## explore more possibilities, perhaps exploiting multiple cores to do
## this simultaneously.

library(MASS)

n <- 2500

cov.uy2 <- 0.5
var.u <- 0.1
mu <- c(1,1,0)

Sigma <- matrix(c(1.0,0.85,cov.uy2,
                  0.85,1.0,0.0,
                  cov.uy2,0.0,1.0),
                3,3,
                byrow=TRUE)

foo <- mvrnorm(n = n,
               mu,
               Sigma)

X <- 2*pnorm(foo[,1],mean=mu[1],sd=sqrt(Sigma[1,1])) -1
W <- 2*pnorm(foo[,2],mean=mu[2],sd=sqrt(Sigma[2,2])) -1
U <- foo[,3]

## h0 is the instrumental DGP function - try changing from **1
## (linear), **2 (quadratic), etc.

h0 <- X**4
Y <- h0 + sqrt(var.u)*U

AIC <- numeric()
model <- list()

for(S in 1:5) {
    model[[S]] <- npiv(Y,
                       X,
                       W,
                       K.w.degree=S,   
                       K.w.segments=S,
                       J.x.degree=S,
                       J.x.segments=S)
    AIC[S] <- model[[S]]$AIC.c
}

## Select the model that minimizes AIC (model index will also be the
## value of S used above)

model <- model[[which.min(AIC)]]

## Create a plot of the instrumental regression function, data, and DGP

plot(X,Y,cex=0.25,
     col="lightgrey",
     sub=paste("n = ",format(n,format="d", big.mark=','),sep=""),
     xlab="X",
     ylab="Y")

lines(X[order(X)],h0[order(X)],lty=1,col=1,lwd=1)
lines(X[order(X)],model$h[order(X)],lty=2,col=2,lwd=2)

legend("topleft",c("DGP",paste("NPIV (K.w.degree = ",model$K.w.degree,
                               ", W.knots = ",model$K.w.segments+1,
                               ", J.x.degree = ", model$J.x.degree,
                               ", X.knots = ",model$J.x.segments+1,")",sep="")),
       lty=1:2,
       col=1:2,
       lwd=c(1,2),
       bty="n",
       cex=0.75)
