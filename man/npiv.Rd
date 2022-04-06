\name{npiv}
\alias{npiv}
\alias{npiv.default}
\alias{npiv.formula}
\title{
Nonparametric Instrumental Variable Estimation and Inference
}
\description{
\code{npiv} performs nonparametric two-stage least-squares estimation of, and constructs uniform confidence bands for, a structural function \code{h0} and its derivatives using a B-spline sieve. Sieve dimensions are determined in a data-dependent way if not provided by the user, via the methods described in Chen, Christensen, and Kankanala (2021). This data-driven choice of sieve dimension ensures estimators of \code{h0} and its derivatives converge at the best possible (i.e., minimax) sup-norm rate. The resulting uniform confidence bands for \code{h0} and its derivatives are honest and adaptive; see Chen, Christensen, and Kankanala (2021).

If sieve dimensions are provided by the user, \code{npiv} implements the bootstrap-based procedure of Chen and Christensen (2018) to construct uniform confidence bands for \code{h0} and its derivatives. 

\code{npiv} may also be used for nonparametric estimation and inference on a regression function as a special case.
}
\usage{
npiv(\dots)

\method{npiv}{formula}(formula,
     data,
     newdata,
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
     K.w.smooth=2,
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
   a symbolic description of the model to be fit.
}
\item{data}{ 
   an optional data frame containing the variables in the model. 
}
\item{newdata}{
  an optional data frame in which to look for variables with which to predict (i.e., predictors in \code{X} passed in \code{X.eval} which must contain identically named variables).
}
\item{subset}{
  an optional vector specifying a subset of observations to be used in the fitting process (see additional details about how this argument interacts with data-dependent bases in the \sQuote{Details} section of the \code{\link{model.frame}} documentation).
}
\item{na.action}{
  a function which indicates what should happen when the data contain NAs. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh} default is \code{\link{na.omit}}. Another possible value is \code{NULL}, no action. Value \code{\link{na.exclude}} can be useful.
}
\item{call}{
  the original function call (this is passed internally by \code{\link{npiv}}). It is not recommended that the user set this.
}
\item{Y}{
  dependent variable vector.
}
  \item{X}{
  matrix of endogenous regressors.
}
  \item{W}{
  matrix of instrumental variables. Set \code{W=X} for nonparametric regression.
}
  \item{X.eval}{
  vector of evaluation point(s). Default uses 50 equally spaced points over the support of each \code{X} variable.
}
\item{alpha}{
  nominal size of the uniform confidence bands. Default is \code{0.05} for 95\% uniform confidence bands.
}
  \item{basis}{
  basis type (if \code{X} or \code{W} are multivariate), a character string. Options are:
  
  \code{tensor} tensor product basis. Default option.
  
  \code{additive} additive basis for additively separable models.
  
  \code{glp} generalized B-spline polynomial basis.
}
  \item{boot.num}{
  number of bootstrap replications.
}
  \item{check.is.fullrank}{
  check that \code{X} and \code{W} have full rank. Default is \code{FALSE}.
}
  \item{deriv.index}{
  integer indicating the column of \code{X} for which to compute the derivative.
}
  \item{deriv.order}{
  integer indicating the order of derivative to be computed.
}
  \item{eval.num}{
  number of evaluation points for each \code{X} variable if \code{X.eval} is not provided.
}
  \item{J.x.degree}{
  B-spline degree (integer or vector of integers of length \code{ncol(X)}) for approximating the structural function. Default is \code{degree=3} (cubic B-spline).
}
  \item{J.x.segments}{
  B-spline number of segments (integer or vector of integers of length \code{ncol(X)}) for approximating the structural function. Default is \code{NULL}. If either \code{J.x.segments=NULL} or \code{K.w.segments=NULL}, these are both chosen automatically using \code{\link{npiv_choose_J}}.
}
  \item{K.w.degree}{
  B-spline degree (integer or vector of integers of lenth \code{ncol(W)}) for estimating the nonparametric first-stage. Default is \code{degree=4} (quartic B-spline).
}
  \item{K.w.segments}{
  B-spline number of segments (integer or vector of integers of length \code{ncol(W)}) estimating the nonparametric first stage. Defulat is \code{NULL}. If either \code{J.x.segments=NULL} or \code{K.w.segments=NULL}, these are both chosen automatically using \code{\link{npiv_choose_J}}.
}
  \item{K.w.smooth}{
  non-negative integer. Basis for the nonparametric first-stage uses \eqn{2^{K.w.smooth}} more B-spline segments for each instrument than the basis approximating the structural function. Default is \code{2}. Setting \code{K.w.smooth=0} uses the same number of segments for \code{X} and \code{W}.
}
  \item{knots}{
  knots type, a character string. Options are:

  \code{quantiles} interior knots are placed at equally spaced quantiles (equal number of observations lie in each segment). Default option. 

  \code{uniform} interior knots are placed at equally spoaced intervals over the support of the variable.
}
  \item{progress}{
  whether to display progress bar or not. Default is \code{TRUE}.
}
\item{random.seed}{
  random seed to initialize bootstrap replications for reproducibility. Default is \code{42}.
}
\item{ucb.h}{
  whether to compute a uniform confidence band for the structural function. Default is \code{TRUE}.
}
\item{ucb.deriv}{
 whether to compute a uniform confidence band for the derivative of the structural function. Default is \code{TRUE}.
}
  \item{W.min}{
  lower bound on the support of each \code{W} variable. Default is \code{min(W)}.
}
  \item{W.max}{
  upper bound on the support of each \code{W} variable. Default is \code{max(W)}.
}
  \item{X.min}{
  lower bound on the support of each \code{X} variable. Default is \code{min(X)}.
}
  \item{X.max}{
  upper bound on the support of each \code{X} variable. Default is \code{max(X)}.
}
 \item{\dots}{
   optional arguments
}
}
\details{

\code{npiv} estimates and constructs uniform confidence bands for a nonparametric structural function \eqn{h_0}{h0} and its derivatives in the model 
 \eqn{Y=h_0(X)+U,\quad E[U|W]=0\quad{(\rm almost\, surely).}}{Y=h0(X)+U, E[U|W]=0 (almost surely).}
Estimation is performed using nonparametric two-stage least-squares with a B-spline sieve. The key tuning parameter is the dimension \eqn{J}{J} of the sieve used to approximate \eqn{h_0}{h0}. The dimension is tuned via modifying the number and placement of interior knots in the B-spline basis  (equivalently, the number of segments of the basis). Sieve dimensions can be user-provided or data-determined using the procedure of Chen, Christensen, and Kankanala (2021).

Typical usages mirror \code{ivreg} (see above and below for a list of options and the example at the bottom of this document)
  \preformatted{
    foo <- npiv(y~x|w)
    foo <- npiv(y~x1+x2|w1+w2)
    foo <- npiv(Y=y,X=x,W=w)
  }

\code{npiv} can be used in two ways:

  1. Data-driven sieve dimension is invoked if either \code{K.w.segments} or \code{J.x.segments} are unspecified or \code{NULL} (the default). Sieve dimensions are chosen automatically using \code{npiv_choose_J}. Uniform confidence bands for \eqn{h_0}{h0} and its derivatives are constructed using the data-driven method of Chen, Christensen, and Kankanala (2021).

  2. The user may specify the sieve dimensions of both bases by specifying values for \code{K.w.segments} and \code{J.x.segments}. Uniform confidence bands for \eqn{h_0}{h0} and its derivatives are constructed using the method of Chen and Christensen (2018).

\code{npiv} can also be used for estimation and inference on a nonparametric regression function by setting \code{W=X}.

}
\value{

 \code{npiv} returns a \code{npiv} object. The generic function \code{\link{fitted}} extracts the estimated values for the sample (or evaluation data, if provided), while the generic function \code{\link{residuals}} extracts the sample residuals. The generic function \code{\link{summary}} provides a simple model summary.

The function \code{npiv} returns a list with the following components:

\item{h}{estimated structural function evaluated at the sample data (or evaluation data, if provided).}
\item{residuals}{residuals for the sample data.}
\item{deriv}{estimated derivative of the structural function evaluated at the sample data (or evaluation data, if provided).}
\item{asy.se}{pre-asymptotic standard errors for the estimator of the structural function evaluated at the sample data (or evaluation data, if provided) }
\item{deriv.asy.se}{pre-asymptotic standard errors for the estimator of the derivative of the structural function evaluated at the sample data (or evaluation data, if provided).}
\item{deriv.index}{index for the estimated derivative.}
\item{deriv.order}{order of the estimated derivative.}
\item{K.w.degree}{value of \code{K.w.degree} used.}
\item{K.w.segments}{value of \code{K.w.segments} used (will be data-determined if not provided).}
\item{J.x.degree}{value of \code{J.x.degree} used.}
\item{J.x.segments}{value of \code{J.x.segments} used (will be data-determined if not provided).}
\item{beta}{vector estimated spline coefficients.}
}
\references{

Chen, X. and T. Christensen (2018). \dQuote{Optimal Sup-norm Rates and Uniform Inference on Nonlinear Functionals of Nonparametric IV Regression.} \emph{Quantitative Economics}, \strong{9}(1), 39-85. doi: \url{https://doi.org/10.3982/QE722}

Chen, X., T. Christensen and S. Kankanala (2021). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV.} \url{https://arxiv.org/pdf/2107.11869.pdf}

}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M. Christensen <timothy.christensen@nyu.edu>
}
%%\note{
%%TBD
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{npiv_choose_J}}
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

## Call the npiv() function using the formula interface npiv.formula()
## with specific arguments, provide evaluation data for X via newdata
## argument (equivalent to npiv(Y,X,W,X.eval,...) npiv.default() call)

model <- npiv(Y~X|W,
              newdata=data.frame(X=X.eval),
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
