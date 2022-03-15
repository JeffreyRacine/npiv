\name{npiv_choose_J}
\alias{npiv_choose_J}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Instrumental Variables Estimation Model Complexity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
              eval.num = 50,
              boot.num = 99,
              check.is.fullrank = FALSE,
              chol.pivot = FALSE,
              lambda = sqrt(.Machine$double.eps))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{W}{
%%     ~~Describe \code{W} here~~
}
  \item{X.eval}{
%%     ~~Describe \code{X.eval} here~~
}
  \item{J.x.degree}{
%%     ~~Describe \code{J.x.degree} here~~
}
  \item{K.w.degree}{
%%     ~~Describe \code{K.w.degree} here~~
}
  \item{K.w.smooth}{
%%     ~~Describe \code{K.w.smooth} here~~
}
  \item{knots}{
%%     ~~Describe \code{knots} here~~
}
  \item{basis}{
%%     ~~Describe \code{basis} here~~
}
  \item{eval.num}{
%%     ~~Describe \code{eval.num} here~~
}
  \item{boot.num}{
%%     ~~Describe \code{boot.num} here~~
}
  \item{check.is.fullrank}{
%%     ~~Describe \code{check.is.fullrank} here~~
}
  \item{chol.pivot}{
%%     ~~Describe \code{chol.pivot} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
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
%% ~put references to the literature/web site here ~
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

out <- npiv_choose_J(Y,
                     X,
                     W,
                     X.eval=X.eval,
                     J.x.degree=J.x.degree,
                     K.w.degree=K.w.degree,
                     K.w.smooth=K.w.smooth,
                     knots=knots,
                     basis=basis,
                     eval.num=eval.num,
                     boot.num=boot.num,
                     check.is.fullrank=check.is.fullrank,
                     chol.pivot=chol.pivot,
                     lambda=lambda)
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
