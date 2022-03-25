\name{npiv_choose_J}
\alias{npiv_choose_J}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Instrumental Variables Estimation Model Complexity
}
\description{
This function chooses the data-driven index set \eqn{\hat{\mathcal{J}} =\{J\in \mathcal{T} :0.1(\log\hat J_{\max}) \le J \le \hat J_{\max}}{J ={J∈T :0.1(log J.max) ≤J≤Jmax}} with \eqn{\mathcal{T} = \{2k + 3 : k = 0, 1, 2,\dots, \hat J_{\max}\}}{T={2k + 3 : k = 0, 1, 2,...,J.max}} obtained from \code{\link{npiv_Jhat_max}}
}
\usage{
npiv_choose_J(Y, 
              X,
              W,
              X.eval = NULL,
              J.x.degree = 3,
              K.w.degree = 3,
              K.w.smooth = 1,
              knots = c("uniform", "quantiles"),
              basis = c("tensor", "additive", "glp"),
              X.min = NULL,
              X.max = NULL,
              W.min = NULL,
              W.max = NULL,
              eval.num = 50,
              boot.num = 99,
              check.is.fullrank = FALSE,
              chol.pivot = FALSE,
              lambda = sqrt(.Machine$double.eps),
              progress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{J.x.degree}{
B-spline degree for endogenous predictor matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{X} (default \code{degree=3}, i.e. cubic spline)
}
  \item{K.w.degree}{
B-spline degree for instrument matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{W} (default \code{degree=3}, i.e. cubic spline)
}
  \item{K.w.smooth}{
  Non-negative integer, such that \eqn{2^{K.w.smooth}} is how many times larger the bases for the instruments are than the bases for the endogenous variables (i.e., \code{K.w.smooth=0} uses the same bases for \code{X} and \code{W})
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
  \item{eval.num}{
Number of evaluation observations (for complexity routines where X.eval necessary but not provided)
}
  \item{boot.num}{
Number of bootstrap replications (for complexity routines where X.eval necessary but not provided)
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
TBD
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{J.hat.max}{Description TBD}
\item{J.hat.n}{Description TBD}
\item{J.hat}{Description TBD}
\item{J.tilde}{Description TBD}
\item{J.x.seg}{Description TBD}
\item{K.w.seg}{Description TBD}
\item{theta.star}{Description TBD}
}
\references{
Chen, X. and T. Christensen and S. Kankanala (2021), \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV}, https://arxiv.org/pdf/2107.11869.pdf
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M. Christensen <timothy.christensen@nyu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(MASS)

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

out <- npiv_choose_J(Y,
                     X,
                     W)
out                     
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
