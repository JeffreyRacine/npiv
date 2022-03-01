\name{npiv}
\alias{npiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Instrumental Variable Estimation
}
\description{
...
}
\usage{
npiv(Y,
     X,
     W,
     X.eval = NULL,
     deriv.index = 1,
     deriv.order = 1,
     W.degree = 3,
     W.segments = 1,
     X.degree = 3,
     X.segments = 1,
     knots = c("uniform", "quantiles"),
     basis = c("tensor", "additive", "glp"),
     check.is.fullrank = FALSE,
     chol.pivot = FALSE,
     lambda = sqrt(.Machine$double.eps))
}

\arguments{
  \item{Y}{
Dependent variable vector 
}
  \item{X}{
Endogenous predictor matrix
}
  \item{W}{
Instrument matrix
}
  \item{X.eval}{
Evaluation points for endogenous predictor matrix
}
  \item{deriv.index}{
Column of endogenous predictor matrix for which to compute derivative
}
  \item{deriv.order}{
Order of derivative to be computed
}
  \item{W.degree}{
B-spline degree for instrument matrix
}
  \item{W.segments}{
B-spline number of segments for instrument matrix
}
  \item{X.degree}{
B-spline degree for endogenous predictor matrix
}
  \item{X.segments}{
B-spline number of segments for endogemous predictor matrix
}
  \item{knots}{
Knot type
}
  \item{basis}{
Basis type
}
  \item{check.is.fullrank}{
Check for full rank input matrices
}
  \item{chol.pivot}{
Argument fed to chol()
}
  \item{lambda}{
Argument fed to chol()
}
}
\details{
This function...
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Simulate data using the mvrnorm() function in the MASS package

library(MASS)

set.seed(42)
n <- 1000

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

## h.0 is the instrumental DGP function

h.0 <- sin(pi*X)
Y <- h.0 + sqrt(var.u)*U

## Create evaluation data and instrumental regression function for the
## endogenous predictor (for plotting with lines as this is sorted)

X.eval <- X[order(X)]
h.0 <- sin(pi*X.eval)

## Call the npiv() function with specific arguments

model <- npiv(Y,
              X,
              W,
              X.eval=X.eval,
              deriv.index=1,
              deriv.order=1,
              W.degree=3,
              W.segments=1,
              X.degree=3,
              X.segments=1)

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

lines(X.eval,h.0,lty=1,col=1,lwd=1)
lines(X.eval,model$h,lty=2,col=2,lwd=2)

lines(X.eval,model$h-1.96*model$h.asy.se,col=4,lty=4,lwd=1)
lines(X.eval,model$h+1.96*model$h.asy.se,col=4,lty=4,lwd=1)

legend("topleft",c("DGP",paste("NPIV (W.degree = ",model$W.degree,
                               ", W.knots = ",model$W.segments+1,
                               ", X.degree = ", model$X.degree,
                               ", X.knots = ",model$X.segments+1,")",sep="")),
               lty=1:2,
               col=1:2,
               lwd=c(1,2),
               bty="n",
               cex=0.75)

## Create a plot of the instrumental regression function's first
## partial derivative and its asymptotic standard error bounds

ylim <- c(min(c(model$h.deriv-1.96*model$h.deriv.asy.se,
                model$h.deriv+1.96*model$h.deriv.asy.se)),
          max(c(model$h.deriv-1.96*model$h.deriv.asy.se,
                model$h.deriv+1.96*model$h.deriv.asy.se)))

plot(X.eval,model$h.deriv,type="l",
     ylim=ylim,
     sub=paste("W.degree = ",model$W.degree,
               ", W.knots = ",model$W.segments+1,
               ", X.degree = ", model$X.degree,
               ", X.knots = ",model$X.segments+1, sep=""),
     xlab="X",
     ylab=paste("Order",model$deriv.order,"Derivative"),
     lty=2,
     lwd=2,
     col=2)

lines(X.eval,model$h.deriv-1.96*model$h.deriv.asy.se,col=4,lty=4)
lines(X.eval,model$h.deriv+1.96*model$h.deriv.asy.se,col=4,lty=4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
