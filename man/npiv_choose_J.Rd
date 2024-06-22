\name{npiv_choose_J}
\alias{npiv_choose_J}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data-driven Choice of Sieve Dimension for Nonparametric Instrumental Variables Estimation and Inference
}
\description{
\code{npiv_choose_J} implements the data-driven choice of sieve dimension developed in Chen, Christensen, and Kankanala (2024) for nonparametric instrumental variables estimation using a B-spline sieve. It applies to nonparametric regression as a special case when predictors are exogenous.
}
\usage{
npiv_choose_J(Y, 
              X,
              W,
              X.grid = NULL,
              J.x.degree = 3,
              K.w.degree = 4,
              K.w.smooth = 2,
              knots = c("uniform", "quantiles"),
              basis = c("tensor", "additive", "glp"),
              X.min = NULL,
              X.max = NULL,
              W.min = NULL,
              W.max = NULL,
              grid.num = 50,
              boot.num = 99,
              random.seed = 42,
              check.is.fullrank = FALSE,
              progress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  dependent variable vector.
}
  \item{X}{
  matrix of endogenous regressors.
}
  \item{W}{
  matrix of instrumental variables. Set \code{W=X} for nonparametric regression.
}
  \item{X.grid}{
  vector of grid point(s). Default uses 50 equally spaced points over the support of each \code{X} variable.
}
  \item{J.x.degree}{
  B-spline degree (integer or vector of integers of length \code{ncol(X)}) for approximating the structural function. Default is \code{degree=3} (cubic B-spline).
}
  \item{K.w.degree}{
  B-spline degree (integer or vector of integers of lenth \code{ncol(W)}) for estimating the nonparametric first-stage. Default is \code{degree=4} (quartic B-spline).
}
  \item{K.w.smooth}{
  non-negative integer. Basis for the nonparametric first-stage uses \eqn{2^{K.w.smooth}} more B-spline segments for each instrument than the basis approximating the structural function. Default is \code{2}. Setting \code{K.w.smooth=0} uses the same number of segments for \code{X} and \code{W}.
}
  \item{knots}{
  knots type, a character string. Options are:

  \code{quantiles} interior knots are placed at equally spaced quantiles (equal number of observations lie in each segment). Default option. 

  \code{uniform} interior knots are placed at equally spoaced intervals over the support of the variable. 
}
  \item{basis}{
  basis type (if \code{X} or \code{W} are multivariate), a character string. Options are:
  
  \code{tensor} tensor product basis. Default option.
  
  \code{additive} additive basis for additively separable models.
  
  \code{glp} generalized B-spline polynomial basis.
}
  \item{X.min}{
  lower bound on the support of each \code{X} variable. Default is \code{min(X)}.
}
  \item{X.max}{
  upper bound on the support of each \code{X} variable. Default is \code{max(X)}.
}
  \item{W.min}{
  lower bound on the support of each \code{W} variable. Default is \code{min(W)}.
}
  \item{W.max}{
  upper bound on the support of each \code{W} variable. Default is \code{max(W)}.
}
  \item{grid.num}{
  number of grid points for each \code{X} variable if \code{X.grid} is not provided.
}
  \item{boot.num}{
  number of bootstrap replications.
}
  \item{random.seed}{
  random seed to initialize bootstrap replications for reproducibility. Default is \code{42}.
}
  \item{check.is.fullrank}{
  check that \code{X} and \code{W} have full rank. Default is \code{FALSE}.
}
  \item{progress}{
  whether to display progress bar or not. Default is \code{TRUE}.
}  
}
%%\details{
%%TBD
%%
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{J.hat.max}{largest element of candidate set of sieve dimensions searched over.}
\item{J.hat.n}{second largest element of candidate set of sieve dimensions searched over.}
\item{J.hat}{bootstrap-based Lepski choice of sieve dimension.}
\item{J.tilde}{data-driven choice of sieve dimension using the method of Chen, Christensen, and Kankanala (2024). Minimum of \code{J.hat} and \code{J.hat.n}.}
\item{J.x.seg}{data-driven number of segments for \code{X} using the method of Chen, Christensen, and Kankanala (2024).}
\item{K.w.seg}{data-driven number of segments for \code{W} using the method of Chen, Christensen, and Kankanala (2024).}
\item{theta.star}{Lepski critical value used in determination of \code{J.hat}.}
}
\references{
Chen, X., T. Christensen and S. Kankanala (2024). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric Structural Functions and Elasticities.} \emph{Review of Economic Studies}, \strong{forthcoming}. doi: \url{https://doi.org/10.1093/restud/rdae025}
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M. Christensen <timothy.christensen@nyu.edu>
}
\note{
TBD
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{npiv}}
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
                3,
                3,
                byrow=TRUE)

foo <- mvrnorm(n = n,
               mu,
               Sigma)

X <- 2*pnorm(foo[,1],mean=mu[1],sd=sqrt(Sigma[1,1])) -1
W <- 2*pnorm(foo[,2],mean=mu[2],sd=sqrt(Sigma[2,2])) -1
U <- foo[,3]

## h0 is the true structural function

h0 <- sin(pi*X)
Y <- h0 + sqrt(var.u)*U

npiv_choose_J(Y,X,W)
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
