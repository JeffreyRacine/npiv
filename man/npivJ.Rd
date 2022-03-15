\name{npivJ}
\alias{npivJ}

\title{
Selection of Complexity Parameter for Nonparametric Instrumental Variables Estimation
}
\description{
This function computes \eqn{\hat J} per Chen, X. and T. Christensen and S. Kankanala (2021).
}
\usage{
npivJ(Y,
      X,
      W,
      X.eval = NULL,
      K.w.degree = 3,
      J.x.degree = 3,
      J.x.segments.set = c(2, 4, 8, 16, 32, 64) - 1,
      knots = c("uniform", "quantiles"),
      basis = c("tensor", "additive", "glp"),
      eval.num = 50,
      boot.num = 99,
      alpha = 0.5,
      check.is.fullrank = FALSE,
      chol.pivot = FALSE,
      lambda = sqrt(.Machine$double.eps))
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
  \item{J.x.segments.set}{
B-spline segments set for endogenous predictor matrix, i.e., , i.e., integer/vector specifying the number of segments of the B-spline basis for the the matrix \code{X}(i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier curve)
}
  \item{K.w.segments.set}{
B-spline segments set for instrument matrix, i.e., integer/vector specifying the number of segments of the B-spline basis for the the matrix \code{W}(i.e. number of knots minus one) (default \code{segments=1}, i.e. Bezier curve)
}
  \item{boot.num}{
Number of bootstrap replications (for complexity routines where X.eval necessary but not provided)
}
  \item{eval.num}{
Number of evaluation observations (for complexity routines where X.eval necessary but not provided)
}

  \item{basis}{
Basis type, i.e,, a character string (default \code{basis="tensor"}) indicating whether the additive or tensor product B-spline basis matrix for a multivariate polynomial spline or generalized B-spline polynomial basis should be used
}
  \item{knots}{
Knots type, i.e., a character string (default \code{knots="quantiles"}) specifying where knots are to be placed. \sQuote{quantiles} specifies knots placed at equally spaced quantiles (equal number of observations lie in each segment) and \sQuote{uniform} specifies knots placed at equally spaced intervals
}
  \item{eval.num}{
Number of evaluation observations (for complexity routines where X.eval necessary but not provided)
}

  \item{boot.num}{
Number of bootstrap replications (for complexity routines where X.eval necessary but not provided)
}
  \item{alpha}{
TBD
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

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

## Call the npivJ() function

out <- npivJ(Y,
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
