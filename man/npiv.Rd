\name{npiv}
\alias{npiv}
\alias{npiv.default}
\alias{npiv.formula}
\title{
Nonparametric Instrumental Variable Estimation
}
\description{

This function implements computationally simple, data-driven procedures for estimation and inference on a structural function \eqn{h_0}{h0} and its derivatives in nonparametric models using instrumental variables. The first procedure is a bootstrap-based, data-driven choice of sieve dimension for sieve nonparametric instrumental variables (NPIV) estimators (\code{\link{npiv_choose_J}}). When implemented with this data-driven choice, sieve NPIV estimators of \eqn{h_0}{h0} and its derivatives are adaptive: they converge at the best possible (i.e., minimax) sup-norm rate, without having to know the smoothness of \eqn{h_0}{h0}, degree of endogeneity of the regressors, or instrument strength. The second procedure is a data-driven approach for constructing honest and adaptive uniform confidence bands (UCBs) for \eqn{h_0}{h0} and its derivatives. The data-driven UCBs guarantee coverage for \eqn{h_0}{h0} and its derivatives uniformly over a generic class of data-generating processes (honesty) and contract at, or within a logarithmic factor of, the minimax sup-norm rate (adaptivity). As such, these data-driven UCBs deliver asymptotic efficiency gains relative to UCBs constructed via the usual approach of undersmoothing. In addition, both procedures apply to nonparametric regression as a special case. 

}
\usage{
npiv(\dots)

\method{npiv}{formula}(formula,
     data,
     subset,
     na.action,
     call,
     \dots)

\method{npiv}{default}(Y,
     X,
     W,
     X.eval=NULL,
     alpha=0.05,
     basis=c("tensor","additive","glp"),
     boot.num=99,
     check.is.fullrank=FALSE,
     deriv.index=1,
     deriv.order=1,
     eval.num=100,
     J.x.degree=3,
     J.x.segments=NULL,
     K.w.degree=4,
     K.w.segments=NULL,
     K.w.smooth=1,
     knots=c("uniform","quantiles"),
     progress=TRUE,
     random.seed=42,
     ucb.h=TRUE,
     ucb.deriv=TRUE,
     W.max=NULL,
     W.min=NULL,
     X.min=NULL,
     X.max=NULL,
     \dots)

}

\arguments{
\item{formula}{ 
   a symbolic description of the model to be fit 
}
\item{data}{ 
   an optional data frame containing the variables in the
  model 
}
\item{subset}{
  an optional vector specifying a subset of observations to be used in the fitting process (see additional details about how this argument interacts with data-dependent bases in the \sQuote{Details} section of the \code{\link{model.frame}} documentation)
}
\item{na.action}{
  a function which indicates what should happen when the data contain NAs. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh} default is \code{\link{na.omit}}. Another possible value is \code{NULL}, no action. Value \code{\link{na.exclude}} can be useful.
}
  \item{call}{
    the original function call (this is passed internally by
    \code{\link{npiv}} - it is not recommended that the user set this)
  }

  \item{Y}{
dependent variable vector \code{Y}
}
  \item{X}{
endogenous predictor matrix \code{X}
}
  \item{W}{
instrument matrix \code{W}
}
  \item{X.eval}{
evaluation points for endogenous predictor matrix \code{X}
}
\item{alpha}{
value lying in \code{(0,0.5)} specifying the nominal size of the uniform confidence bounds (defaults to \code{0.05} which produces 95\% coverage bounds )
}
  \item{basis}{
basis type, i.e,, a character string (default \code{basis="tensor"}) indicating whether the additive or tensor product B-spline basis matrix for a multivariate polynomial spline or generalized B-spline polynomial basis should be used
}
  \item{boot.num}{
number of bootstrap replications (for complexity routines where X.eval necessary but not provided)
}
  \item{check.is.fullrank}{
check for full rank input matrices for \code{X} and \code{W}
}
  \item{deriv.index}{
column of endogenous predictor matrix \code{X} for which to compute derivative (scalar)
}
  \item{deriv.order}{
order of derivative to be computed (scalar)
}
  \item{eval.num}{
number of evaluation observations (for complexity routines where X.eval necessary but not provided)
}
  \item{J.x.degree}{
B-spline degree for endogenous predictor matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{X} (default \code{degree=3}, i.e. cubic spline)
}
  \item{J.x.segments}{
B-spline number of segments for endogenous predictor matrix, i.e., , i.e., integer/vector specifying the number of segments of the B-spline basis for the the matrix \code{X}(i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier curve)
}
  \item{K.w.degree}{
B-spline degree for instrument matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{W} (default \code{degree=3}, i.e. cubic spline)
}
  \item{K.w.segments}{
B-spline number of segments for instrument matrix, i.e., integer/vector specifying the number of segments of the B-spline basis for the the matrix \code{W}(i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier curve)
}
  \item{K.w.smooth}{
B-spline smoothness parameter (for complexity routines where X.eval necessary but not provided)
}
  \item{knots}{
knots type, i.e., a character string (default \code{knots="quantiles"}) specifying where knots are to be placed. \sQuote{quantiles} specifies knots placed at equally spaced quantiles (equal number of observations lie in each segment) and \sQuote{uniform} specifies knots placed at equally spaced intervals
}
  \item{progress}{
  logical indicating whether to display progress bar or not (defaults to \code{TRUE})
}
\item{random.seed}{
 random seed (defaults to 42) used for replication of bootstrap resamples in applied settings
}
\item{ucb.h}{
logical indicating whether to compute bootstrap uniform confidence bounds for the instrumental function \code{h} or not (defaults to \code{TRUE})
}
\item{ucb.deriv}{
logical indicating whether to compute bootstrap uniform confidence bounds for the specified derivative of the instrumental function \code{h} or not (defaults to \code{TRUE})
}
  \item{W.min}{
  the lower bound on which to construct the spline for predictors in \code{W} - defaults to \code{min(W)}
}
   \item{W.max}{
   the upper bound on which to construct the spline for predictors in \code{W} - defaults to \code{max(W)}
}
  \item{X.min}{
  the lower bound on which to construct the spline for predictors in \code{X} - defaults to \code{min(X)}
}
   \item{X.max}{
   the upper bound on which to construct the spline for predictors in \code{X} - defaults to \code{max(X)}
}
 \item{\dots}{
   optional arguments
}
}
\details{

This function is written to address nonparametric instrumental variable models of the form
 \eqn{Y=h_0(X)+U,\quad E[U|W]=0\quad{(\rm almost\, surely)}}{Y=h0(X)+U, E[U|W]=0 (almost surely)}
where estimators are known to be particularly sensitive to tuning parameters. Data-driven procedures for selecting these tuning parameters can render much needed transparency to nonparametric estimation and inference procedures. The key tuning parameter is the sieve dimension \code{J}, which is tuned via modifying the number and placement of interior knots in the B-spline basis (defaults to a cubic B-spline formulation with the number of knots/segments selected via a sup-norm data-driven procedure).

Typical usages are (see above and below for a list of options and the example at the bottom of this document)
  \preformatted{
    foo <- npiv(y~x|w)
    foo <- npiv(y~x1+x2|w1+w2)
    foo <- npiv(Y=y,X=x,W=w)
  }

This function can be used in two ways:

  i) manually where one specifies the number of segments (i.e., number of knots-1) by specifying particular values for \code{K.w.segments} and \code{J.x.segments} (not recommended for general use)

  ii) using the data-driven procedure for selecting model complexity (which will be invoked if either or both \code{K.w.segments} and \code{J.x.segments} are unspecified or set explicitly to NULL (default))

In general, \code{X} and \code{W} are not equal, otherwise if \code{X}=\code{W} this is not instrumental regression, rather is is conditional mean regression.

}
\value{

 \code{npiv} returns a \code{npiv} object. The generic function \code{\link{fitted}} extracts the estimated values for the sample (or evaluation data, if provided), while the generic function \code{\link{residuals}} extracts the sample residuals. The generic function \code{\link{summary}} provides a simple model summary.

The function \code{npiv} returns a list with the following components:

\item{h}{h returns the estimated instrumental regression function for the sample data (or evaluation data, if provided) }
\item{residuals}{residuals returns the residuals for the sample data}
\item{deriv}{deriv returns the estimated derivative of the instrumental regression function for the sample data (or evaluation data, if provided) }
\item{asy.se}{asy.se returns the value of asymptotic standard errors for the sample data (or evaluation data, if provided) }
\item{deriv.asy.se}{asy.se.deriv returns the derivative asymptotic standard errors for the sample data (or evaluation data, if provided)  }
\item{deriv.index}{deriv.index returns the index provided and used for computing the derivative}
\item{deriv.order}{deriv.order returns the order of the estimated derivative }
\item{K.w.degree}{K.w.degree returns the value used }
\item{K.w.segments}{K.w.segments returns the value used }
\item{J.x.degree}{J.x.degree returns the value used }
\item{J.x.segments}{J.x.segments returns the value used }
\item{beta}{beta returns the value of the spline coefficients}
}
\references{

Chen, X. and T. Christensen and S. Kankanala (2021), \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV}, \url{https://arxiv.org/pdf/2107.11869.pdf}

}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M. Christensen <timothy.christensen@nyu.edu>
}
\note{
TBD
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{npiv_choose_J}}
}
\examples{
## Simulate data using the mvrnorm() function in the MASS package

library(MASS)

set.seed(42)
n <- 1000

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

## Create evaluation data and instrumental regression function for the
## endogenous predictor (for plotting with lines as this is sorted)

X.eval <- seq(min(X),max(X),length=100)
h0 <- sin(pi*X.eval)

## Call the npiv() function using the formula interface with specific arguments

model <- npiv(Y~X|W,
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

lines(X.eval,model$h-1.96*model$asy.se,col=4,lty=4,lwd=1)
lines(X.eval,model$h+1.96*model$asy.se,col=4,lty=4,lwd=1)

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

ylim <- c(min(c(model$deriv-1.96*model$deriv.asy.se,
                model$deriv+1.96*model$deriv.asy.se)),
          max(c(model$deriv-1.96*model$deriv.asy.se,
                model$deriv+1.96*model$deriv.asy.se)))

plot(X.eval,model$deriv,type="l",
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

lines(X.eval,model$deriv-1.96*model$deriv.asy.se,col=4,lty=4)
lines(X.eval,model$deriv+1.96*model$deriv.asy.se,col=4,lty=4)
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
