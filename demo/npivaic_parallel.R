## Model selection via Hurvich, Siminoff and Tsai's AIC.c
## criterion. Note we set the degree and number of knots/segments to
## be equal for both the Psi.x and B.w bases. You could, naturally,
## explore more possibilities, perhaps exploiting multiple cores to do
## this simultaneously.

library(npiv)
library(doParallel)
library(MASS)

## Setting parallel.cores to NULL and parallel to TRUE will use all
## physical (not logical) cores. On multi-core systems you can
## experiment by setting parallel to TRUE and change parallel.cores
## from 1, 2, to the maximum number of cores. However, you can exhaust
## memory fairly quickly as each instance will demand its own memory
## resources.

parallel <- TRUE # FALSE
parallel.cores <- NULL

n <- 10000

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
## (linear), **2 (quadratic),cos(2*pi*X), sin(2*pi*X) etc.

h0 <- sin(pi*X) # X**4
Y <- h0 + sqrt(var.u)*U

AIC <- numeric()

## Note K >= J is necessary (tested for, i.e. K.w.degree+K.w.segments >=
## J.x.degree+J.x.segments). You could loop over either spline degree or number
## of knots (#knots = #segments+1) for either the W basis or X basis or both,
## _providing_ you ensure that J >= K is met (otherwise the function will test
## for this condition and halt)

t0 <- Sys.time()

S.max <- 10
cl<-makeCluster(if(is.null(parallel.cores)){detectCores(logical=FALSE)}else{parallel.cores})
registerDoParallel(cl)
output <- foreach(S=1:S.max,.verbose=FALSE) %dopar% {
    model <- npiv::npivaic(Y,
                     X,
                     W,
                     K.w.degree=3,   
                     K.w.segments=S,
                     J.x.degree=3,
                     J.x.segments=S)
    list(model)
}
stopCluster(cl)
## Elapsed time
Sys.time()-t0

## Select the model that minimizes AIC (model index will also be the
## value of S used above)

for(S in 1:S.max) AIC[S] <- output[[S]][[1]]$AIC.c
model <- output[[which.min(AIC)]][[1]]

## Create a plot of the instrumental regression function, data, and DGP

plot(X,Y,cex=0.25,
     col="lightgrey",
     sub=paste("n = ",format(n,format="d", big.mark=','),
               " (knots searched over: 1,2,...,",S.max,")",sep=""),
     main="Data-Driven Knot Selection (AIC, Hurvich et al (1998))",
     xlab="X",
     ylab="Y")

lines(X[order(X)],h0[order(X)],lty=1,col=1,lwd=1)
lines(X[order(X)],model$h[order(X)],lty=2,col=2,lwd=2)

legend("topleft",c("DGP (IV Function h0)",paste("NPIV (K.w.degree = ",model$K.w.degree,
                               ", W.knots = ",model$K.w.segments+1,
                               ", J.x.degree = ", model$J.x.degree,
                               ", X.knots = ",model$J.x.segments+1,")",sep="")),
       lty=1:2,
       col=1:2,
       lwd=c(1,2),
       bty="n",
       cex=0.75)
