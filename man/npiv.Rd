\name{npiv}
\alias{npiv}
\title{
Nonparametric Instrumental Variable Estimation
}
\description{

This function implements computationally simple, data-driven procedures for estimation and inference on a structural function \code{h0()} and its derivatives in nonparametric models using instrumental variables. The first procedure is a bootstrap-based, data-driven choice of sieve dimension for sieve nonparametric instrumental variables (NPIV) estimators (\code{\link{npiv_Jhat_max}}, \code{\link{npiv_choose_J}}). When implemented with this data-driven choice, sieve NPIV estimators of \code{h0()} and its derivatives are adaptive: they converge at the best possible (i.e., minimax) sup-norm rate, without having to know the smoothness of \code{h0()}, degree of endogeneity of the regressors, or instrument strength. The second procedure is a data-driven approach for constructing honest and adaptive uniform confidence bands (UCBs) for \code{h0()} and its derivatives. The data-driven UCBs guarantee coverage for \code{h0()} and its derivatives uniformly over a generic class of data-generating processes (honesty) and contract at, or within a logarithmic factor of, the minimax sup-norm rate (adaptivity). As such, these data-driven UCBs deliver asymptotic efficiency gains relative to UCBs constructed via the usual approach of undersmoothing. In addition, both procedures apply to nonparametric regression as a special case. 

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
     K.w.smooth = 1,
     J.x.degree = 3,
     J.x.segments = 1,
     boot.num = 99,
     eval.num = 100,
     knots = c("uniform", "quantiles"),
     basis = c("tensor", "additive", "glp"),
     X.min = NULL,
     X.max = NULL,
     W.min = NULL,
     W.max = NULL,
     check.is.fullrank = FALSE,
     chol.pivot = FALSE,
     lambda = sqrt(.Machine$double.eps),
     progress = TRUE)
}

\arguments{
  \item{Y}{
Dependent variable vector \code{Y}
}
  \item{X}{
Endogenous predictor matrix \code{X}
}
  \item{W}{
Instrument matrix \code{W}
}
  \item{X.eval}{
Evaluation points for endogenous predictor matrix \code{X}
}
  \item{deriv.index}{
Column of endogenous predictor matrix \code{X} for which to compute derivative (scalar)
}
  \item{deriv.order}{
Order of derivative to be computed (scalar)
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
  \item{J.x.degree}{
B-spline degree for endogenous predictor matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{X} (default \code{degree=3}, i.e. cubic spline)
}
  \item{J.x.segments}{
B-spline number of segments for endogenous predictor matrix, i.e., , i.e., integer/vector specifying the number of segments of the B-spline basis for the the matrix \code{X}(i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier curve)
}
  \item{boot.num}{
Number of bootstrap replications (for complexity routines where X.eval necessary but not provided)
}
  \item{eval.num}{
Number of evaluation observations (for complexity routines where X.eval necessary but not provided)
}

  \item{knots}{
Knots type, i.e., a character string (default \code{knots="quantiles"}) specifying where knots are to be placed. \sQuote{quantiles} specifies knots placed at equally spaced quantiles (equal number of observations lie in each segment) and \sQuote{uniform} specifies knots placed at equally spaced intervals
}
  \item{basis}{
Basis type, i.e,, a character string (default \code{basis="tensor"}) indicating whether the additive or tensor product B-spline basis matrix for a multivariate polynomial spline or generalized B-spline polynomial basis should be used
}
  \item{X.min}{
  the lower bound on which to construct the spline for predictors in \code{X} - defaults to \code{min(X)}
}
   \item{X.max}{
   the upper bound on which to construct the spline for predictors in \code{X} - defaults to \code{max(X)}
}
  \item{W.min}{
  the lower bound on which to construct the spline for predictors in \code{W} - defaults to \code{min(W)}
}
   \item{W.max}{
   the upper bound on which to construct the spline for predictors in \code{W} - defaults to \code{max(W)}
}
  \item{check.is.fullrank}{
Check for full rank input matrices for \code{X} and \code{W}
}
  \item{chol.pivot}{
Argument fed to \code{chol()} (see \code{\link{chol}})
}
  \item{lambda}{
Argument fed to \code{chol()} (see \code{\link{chol}})
}
  \item{progress}{
  logical indicating whether to display progress bar or not (defaults to \code{TRUE})
}  
}
\details{
This function can be used in two ways, i) manually where one specifies the number of segments (i.e., number of knots-1) by specifying particular values for \code{K.w.segments} and \code{J.x.segments}, or ii) using the data-driven procedure for selecting model complexity (which will be invoked if either or both \code{K.w.segments} and \code{J.x.segments} are unspecified or set explicitly to NULL (default)).

In general, \code{X} and \code{W} are not equal, otherwise if \code{X}=\code{W} this is not instrumental regression, rather is is conditional mean regression. If the user passes \code{X}=\code{W} then this is recognized as standard nonlinear regression and the B-spline complexity is chosen by Hurvich, Simonoff and Tsai's (1998) \dQuote{corrected} AIC nonlinear model selection procedure.

}
\value{

 \code{npiv} returns a \code{npiv} object. The generic function \code{\link{fitted}} extracts the estimated values for the sample (or evaluation data, if provided), while the generic function \code{\link{residuals}} extracts the sample residuals.

The function \code{npiv} returns a list with the following components:

\item{fitted}{fitted returns the estimated instrumental regression function for the sample data (or evaluation data, if provided) }
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
\item{B.w}{B.w returns the value of the spline bases for \code{W}}
\item{Psi.x}{Psi.x returns the value of the spline bases for \code{X}}
\item{Psi.x.deriv}{Psi.x.deriv returns the value of the spline bases for the derivatives}
\item{AIC.c}{Hurvich, Simonoff and Tsai's (1998) corrected AIC criterion}
}
\references{

Chen, X. and T. Christensen and S. Kankanala (2021), \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV}, https://arxiv.org/pdf/2107.11869.pdf

Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998), \dQuote{Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion,} Journal of the Royal Statistical Society B, 60, 271-293.

}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M. Christensen <timothy.christensen@nyu.edu>
}
\note{
TBD
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{npiv_choose_J}},  \code{\link{npiv_Jhat_max}},  \code{\link{npivJ}}
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

lines(X.eval,model$fitted-1.96*model$asy.se,col=4,lty=4,lwd=1)
lines(X.eval,model$fitted+1.96*model$asy.se,col=4,lty=4,lwd=1)

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
