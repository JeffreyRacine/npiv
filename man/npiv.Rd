\name{npiv}
\alias{npiv}
\alias{npiv.default}
\alias{npiv.formula}
\title{
Nonparametric Instrumental Variable Estimation and Inference
}
\description{
\code{npiv} performs nonparametric a structural function \code{h0} and its derivatives using a B-spline sieve. It also constructs uniform confidence bands for \code{h0} and its derivative. 

Sieve dimensions are determined in a data-dependent way if not provided by the user, via the methods described in Chen, Christensen, and Kankanala (2024). This data-driven choice of sieve dimension ensures estimators of \code{h0} and its derivatives converge at the optimal sup-norm rate. The resulting uniform confidence bands for \code{h0} and its derivatives also converge at the minimax rate up to log factors; see Chen, Christensen, and Kankanala (2024).

If sieve dimensions are provided by the user, \code{npiv} implements the bootstrap-based procedure of Chen and Christensen (2018) to construct uniform confidence bands based on undersmoothing for \code{h0} and its derivatives. 

The methods in \code{npiv} apply to estimation and inference on a nonparametric regression function as a special case.
}
\usage{
npiv(\dots)

\method{npiv}{formula}(formula,
     data=NULL,
     newdata=NULL,
     subset=NULL,
     na.action="na.omit",
     call,
     \dots)

\method{npiv}{default}(Y,
     X,
     W,
     X.eval=NULL,
     X.grid=NULL,
     alpha=0.05,
     basis=c("tensor","additive","glp"),
     boot.num=99,
     check.is.fullrank=FALSE,
     deriv.index=1,
     deriv.order=1,
     grid.num=50,
     J.x.degree=3,
     J.x.segments=NULL,
     K.w.degree=4,
     K.w.segments=NULL,
     K.w.smooth=2,
     knots=c("uniform","quantiles"),
     progress=TRUE,
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
  optional matrix of evaluation data for the endogenous regressors.
}
  \item{X.grid}{
  optional vector of grid points for \code{X} when determining model complexity. Default (\code{X.grid=NULL}) uses 50 equally spaced points (can be changed in \code{grid.num}) over the support of each \code{X} variable.
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
  \item{grid.num}{
  number of grid points for each \code{X} variable if \code{X.grid} is not provided.
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

  \code{quantiles} interior knots are placed at equally spaced quantiles (equal number of observations lie in each segment). 

  \code{uniform} interior knots are placed at equally spoaced intervals over the support of the variable. Default option.
}
  \item{progress}{
  whether to display progress bar or not. Default is \code{TRUE}.
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
Estimation is performed using nonparametric two-stage least-squares with a B-spline sieve. The key tuning parameter is the dimension \eqn{J}{J} of the sieve used to approximate \eqn{h_0}{h0}. The dimension is tuned via modifying the number and placement of interior knots in the B-spline basis  (equivalently, the number of segments of the basis). Sieve dimensions can be user-provided or data-determined using the procedure of Chen, Christensen, and Kankanala (2024).

Typical usages mirror \code{ivreg} (see above and below for a list of options and the example at the bottom of this document)
  \preformatted{
    foo <- npiv(y~x|w)
    foo <- npiv(y~x1+x2|w1+w2)
    foo <- npiv(Y=y,X=x,W=w)
  }

\code{npiv} can be used in two ways:

  1. Data-driven sieve dimension is invoked if either \code{K.w.segments} or \code{J.x.segments} are unspecified or \code{NULL} (the default). Sieve dimensions are chosen automatically using \code{npiv_choose_J}. Uniform confidence bands for \eqn{h_0}{h0} and its derivatives are constructed using the data-driven method of Chen, Christensen, and Kankanala (2024).

  2. The user may specify the sieve dimensions of both bases by specifying values for \code{K.w.segments} and \code{J.x.segments}. Uniform confidence bands for \eqn{h_0}{h0} and its derivatives are constructed using the method of Chen and Christensen (2018).

\code{npiv} can also be used for estimation and inference on a nonparametric regression function by setting \code{W=X}.

}
\value{

 \code{npiv} returns a \code{npiv} object. The generic function \code{\link{fitted}} extracts the estimated values for the sample (or evaluation data, if provided), while the generic function \code{\link{residuals}} extracts the sample residuals. The generic function \code{\link{summary}} provides a simple model summary. The generic function \code{\link{plot}} also plots the estimated function and derivative, together with uniform confidence bands.

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
\item{beta}{vector of estimated spline coefficients.}
}
\references{

Chen, X. and T. Christensen (2018). \dQuote{Optimal Sup-norm Rates and Uniform Inference on Nonlinear Functionals of Nonparametric IV Regression.} \emph{Quantitative Economics}, \strong{9}(1), 39-85. \doi{10.3982/QE722}

Chen, X., T. Christensen and S. Kankanala (2024). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric Structural Functions and Elasticities.} \emph{Review of Economic Studies}, \strong{forthcoming}. \doi{10.1093/restud/rdae025}

}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy Christensen <timothy.christensen@yale.edu>
}

\seealso{
\code{\link{npiv_choose_J}}
}

\examples{
## load data
data("Engel95", package = "npiv")

## sort on logexp (the regressor) for plotting purposes
Engel95 <- Engel95[order(Engel95$logexp),] 
attach(Engel95)

## Estimate the Engel curve for food using logwages as an instrument
fm1 <- npiv(food ~ logexp | logwages)

## Plot the estimated Engel curve and data-driven uniform confidence bands
plot(logexp,food,
     ylab="Food Budget Share",
     xlab="log(Total Household Expenditure)",
     xlim=c(4.75, 6.25),
     ylim=c(0, 0.4),
     main="",
     type="p",
     cex=.5,
     col="lightgrey")
lines(logexp,fm1$h,col="blue",lwd=2,lty=1)
lines(logexp,fm1$h.upper,col="blue",lwd=2,lty=2)
lines(logexp,fm1$h.lower,col="blue",lwd=2,lty=2)

## Estimate the Engel curve using pre-specified sieve dimension 
## (dimension 5 for logexp, dimension 9 for logwages)
fm2 <- npiv(food ~ logexp | logwages,
            J.x.segments = 2,
            K.w.segments = 5)

## Plot uniform confidence bands based on undersmoothing
lines(logexp,fm2$h.upper,col="red",lwd=2,lty=2)
lines(logexp,fm2$h.lower,col="red",lwd=2,lty=2)

## Plot pointwise confidence bands based on pre-asymptotic standard errors
lines(logexp,fm2$h+1.96*fm2$asy.se,col="red",lwd=2,lty=3)
lines(logexp,fm2$h-1.96*fm2$asy.se,col="red",lwd=2,lty=3)

legend("topright",
       legend=c("Data-driven Estimate",
                "Data-driven UCBs",
                "Undersmoothed UCBs",
                "Pointwise CBs"),
       col=c("blue","blue","red","red"),
       lty=c(1,2,2,3),
       lwd=c(2,2,2,2))

## Plot the data-driven estimate of the derivative of the Engel curve
plot(logexp,fm1$deriv,col="blue",lwd=2,lty=1,type="l",
     ylab="Derivative of Food Budget Share",
     xlab="log(Total Household Expenditure)",
     xlim=c(4.75, 6.25),
     ylim=c(-1,1))

## Plot data-driven uniform confidence bands for the derivative
lines(logexp,fm1$h.upper.deriv,col="blue",lwd=2,lty=2)
lines(logexp,fm1$h.lower.deriv,col="blue",lwd=2,lty=2)

## Plot uniform confidence bands based on undersmoothing
lines(logexp,fm2$h.upper.deriv,col="red",lwd=2,lty=2)
lines(logexp,fm2$h.lower.deriv,col="red",lwd=2,lty=2)

## Plot pointwise confidence bands based on pre-asymptotic standard errors
lines(logexp,fm2$deriv+1.96*fm2$deriv.asy.se,col="red",lwd=2,lty=3)
lines(logexp,fm2$deriv-1.96*fm2$deriv.asy.se,col="red",lwd=2,lty=3)

legend("topright",
       legend=c("Data-driven Estimate",
                "Data-driven UCBs",
                "Undersmoothed UCBs",
                "Pointwise CBs"),
       col=c("blue","blue","red","red"),
       lty=c(1,2,2,3),
       lwd=c(2,2,2,2))
}
