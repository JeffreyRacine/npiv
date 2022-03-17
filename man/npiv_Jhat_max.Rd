\name{npiv_Jhat_max}
\alias{npiv_Jhat_max}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Instrumental Variables Estimation Maximum Model Complexity
}
\description{
This function computes the maximum size of the feasible set of seive dimension values \eqn{\hat J_{\max}}{J.max}, i.e., 
\eqn{\mathcal{T} = \{2k + 3 : k = 0, 1, 2,\dots, \hat J_{\max}\}}{T={2k + 3 : k = 0, 1, 2,...,J.max}} for B-spline bases.
}
\usage{
npiv_Jhat_max(X,
              W,
              J.x.degree = 3,
              K.w.degree = 3,
              K.w.smooth = 1,
              knots = c("uniform", "quantiles"),
              basis = c("tensor", "additive", "glp"),
              check.is.fullrank = FALSE,
              chol.pivot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Endogenous predictor matrix \code{X}
}
  \item{W}{
Instrument matrix \code{W}
} 
  \item{J.x.degree}{
B-spline degree for endogenous predictor matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{X} (default \code{degree=3}, i.e. cubic spline)
}
  \item{K.w.degree}{
B-spline degree for instrument matrix, i.e., integer/vector specifying the polynomial degree of the B-spline basis for \code{W} (default \code{degree=3}, i.e. cubic spline)
}
  \item{K.w.smooth}{
TBD
}
  \item{knots}{
Knots type, i.e., a character string (default \code{knots="quantiles"}) specifying where knots are to be placed. \sQuote{quantiles} specifies knots placed at equally spaced quantiles (equal number of observations lie in each segment) and \sQuote{uniform} specifies knots placed at equally spaced intervals
}
  \item{basis}{
Basis type, i.e,, a character string (default \code{basis="tensor"}) indicating whether the additive or tensor product B-spline basis matrix for a multivariate polynomial spline or generalized B-spline polynomial basis should be used
}
   \item{check.is.fullrank}{
Check for full rank input matrices for \code{X} and \code{W}
}
  \item{chol.pivot}{
Argument fed to \code{chol()} (see \code{\link{chol}})
}
}
\details{
Compute \eqn{\hat J_{\max}}{J.max} per Chen, Christensen and Kankanala (2021).
}
\value{
%%  ~Describe the value returned
\item{J.x.segments.set}{Description TBD}
\item{K.w.segments.set}{Description TBD}
\item{J.hat.max}{Description TBD}
\item{alpha.hat}{Description TBD}
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

out <- npiv_Jhat_max(X,
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
