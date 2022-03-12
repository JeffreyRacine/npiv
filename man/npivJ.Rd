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
  \item{K.w.degree}{
B-spline degree for instrument matrix
}
  \item{J.x.degree}{
B-spline degree for endogenous predictor matrix
}
  \item{J.x.segments.set}{
TBD  
}
  \item{knots}{
Knot type
}
  \item{basis}{
Basis type
}
  \item{eval.num}{
Number of evaluation points for grid for \code{X}
}
  \item{boot.num}{
Number of bootstrap replications
}
  \item{alpha}{
TBD
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
library(MASS)

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
