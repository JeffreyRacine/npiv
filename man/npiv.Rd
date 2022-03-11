\name{npiv}
\alias{npiv}
\title{
Nonparametric Instrumental Variable Estimation
}
\description{

This function implements a computationally simple, data-driven procedures for estimation and inference on a structural function \code{h0} and its derivatives in nonparametric models using instrumental variables. The first procedure is a bootstrap-based, data-driven choice of sieve dimension for sieve nonparametric instrumental variables (NPIV) estimators. When implemented with this data-driven choice, sieve NPIV estimators of \code{h0} and its derivatives are adaptive: they converge at the best possible (i.e., minimax) sup-norm rate, without having to know the smoothness of \code{h0}, degree of endogeneity of the regressors, or instrument strength. The second procedure is a data-driven approach for constructing honest and adaptive uniform confidence bands (UCBs) for \code{h0} and its derivatives. The data-driven UCBs guarantee coverage for \code{h0} and its derivatives uniformly over a generic class of data-generating processes (honesty) and contract at, or within a logarithmic factor of, the minimax sup-norm rate (adaptivity). As such, these data-driven UCBs deliver asymptotic efficiency gains relative to UCBs constructed via the usual approach of undersmoothing. In addition, both procedures apply to nonparametric regression as a special case. 

}
\usage{
npiv(Y,
     X,
     W,
     X.eval = NULL,
     deriv.index = 1,
     deriv.order = 1,
     K.w.degree = 3,
     K.w.segments = 1,
     J.x.degree = 3,
     J.x.segments = 1,
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
  \item{K.w.degree}{
B-spline degree for instrument matrix
}
  \item{K.w.segments}{
B-spline number of segments for instrument matrix
}
  \item{J.x.degree}{
B-spline degree for endogenous predictor matrix
}
  \item{J.x.segments}{
B-spline number of segments for endogenous predictor matrix
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

Note - the asymptotic covariance matrix is computed for all \eqn{n} sample realizations  by default, unless an evaluation matrix for \code{X} is provided. This matrix will be \eqn{n\times n} which can quickly exhaust memory in a system. To avoid this simply provide an evaluation matrix for X (\code{X.eval}) that contains, say, a grid of \eqn{100} points (one typically only wants the asymptotic variance for plotting confidence intervals). The covariance matrix will now be of dimension \eqn{100\times 100} thereby allowing one to run the function with a large number of sample realizations.

}
\value{
This function returns a list with the following components:

\item{h}{h returns the estimated instrumental regression function}
\item{h.deriv}{h.deriv returns the estimated derivative of the instrumental regression function }
\item{h.asy.se}{asy.se returns the value of }
\item{h.deriv.asy.se}{asy.se.deriv returns the value of }
\item{deriv.index}{deriv.index returns the index provided and used for computing the derivative}
\item{deriv.order}{deriv.order returns the order of the estimated derivative }
\item{K.w.degree}{K.w.degree returns the value of }
\item{K.w.segments}{K.w.segments returns the value of }
\item{J.x.degree}{J.x.degree returns the value of }
\item{J.x.segments}{J.x.segments returns the value of }
\item{beta}{beta returns the value of }
\item{B.w}{B.w returns the value of }
\item{Psi.x}{Psi.x returns the value of }
\item{Psi.x.deriv}{Psi.x.deriv returns the value of }
\item{residuals}{residuals returns the residuals for the sample data}
\item{AIC.c}{Hurvich, Siminoff and Tsai's (1998) corrected AIC criterion}
}
\references{

Chen, X. and T. Christensen and S. Kankanala (2021), \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV}, https://arxiv.org/pdf/2107.11869.pdf

}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>
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

## h0 is the instrumental DGP function

h0 <- sin(pi*X)
Y <- h0 + sqrt(var.u)*U

## Create evaluation data and instrumental regression function for the
## endogenous predictor (for plotting with lines as this is sorted)

X.eval <- seq(min(X),max(X),length=100)
h0 <- sin(pi*X.eval)

## Call the npiv() function with specific arguments

model <- npiv(Y,
              X,
              W,
              X.eval=X.eval,
              deriv.index=1,
              deriv.order=1,
              K.w.degree=3,
              K.w.segments=1,
              J.x.degree=3,
              J.x.segments=1)

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

lines(X.eval,model$h-1.96*model$h.asy.se,col=4,lty=4,lwd=1)
lines(X.eval,model$h+1.96*model$h.asy.se,col=4,lty=4,lwd=1)

legend("topleft",c("DGP",paste("NPIV (K.w.degree = ",model$K.w.degree,
                               ", W.knots = ",model$K.w.segments+1,
                               ", J.x.degree = ", model$J.x.degree,
                               ", X.knots = ",model$J.x.segments+1,")",sep="")),
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
     sub=paste("K.w.degree = ",model$K.w.degree,
               ", W.knots = ",model$K.w.segments+1,
               ", J.x.degree = ", model$J.x.degree,
               ", X.knots = ",model$J.x.segments+1, sep=""),
     xlab="X",
     ylab=paste("Order",model$deriv.order,"Derivative"),
     lty=2,
     lwd=2,
     col=2)

lines(X.eval,model$h.deriv-1.96*model$h.deriv.asy.se,col=4,lty=4)
lines(X.eval,model$h.deriv+1.96*model$h.deriv.asy.se,col=4,lty=4)
abline(h=0,col="grey")

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
