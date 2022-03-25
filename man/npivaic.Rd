\name{npivaic}
\alias{npivaic}

\title{
Nonparametric Instrumental Variables Estimation Model Complexity via AIC.c 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npivaic(Y,
        X,
        W,
        K.w.degree = 3,
        K.w.segments = 1,
        J.x.degree = 3,
        J.x.segments = 1,
        knots = c("uniform", "quantiles"),
        basis = c("tensor", "additive", "glp"),
        X.min = NULL,
        X.max = NULL,
        W.min = NULL,
        W.max = NULL,
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
  \item{K.w.degree}{
%%     ~~Describe \code{K.w.degree} here~~
}
  \item{K.w.segments}{
%%     ~~Describe \code{K.w.segments} here~~
}
  \item{J.x.degree}{
%%     ~~Describe \code{J.x.degree} here~~
}
  \item{J.x.segments}{
%%     ~~Describe \code{J.x.segments} here~~
}
  \item{knots}{
%%     ~~Describe \code{knots} here~~
}
  \item{basis}{
%%     ~~Describe \code{basis} here~~
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
Jeffrey S. Racine <racinej@mcmaster.ca>}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Model selection via Hurvich, Siminoff and Tsai's AIC.c
## criterion. Note we set the degree and number of knots/segments to
## be equal for both the Psi.x and B.w bases. You could, naturally,
## explore more possibilities, perhaps exploiting multiple cores to do
## this simultaneously.

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

## h0 is the instrumental DGP function - try changing from X**1
## (linear), X**2 (quadratic), cos(2*pi*X), sin(2*pi*X) etc.

h0 <- sin(pi*X) # X**4
Y <- h0 + sqrt(var.u)*U

AIC <- numeric()
model <- list()

## Note K >= J is necessary (tested for, i.e. K.w.degree+K.w.segments
## >= J.x.degree+J.x.segments). You could loop over either spline
## degree or number of knots (#knots = #segments+1) for either the W
## basis or X basis or both, _providing_ you ensure that K >= J is met
## (otherwise the function will test for this condition and halt)

S.max <- 10
for(S in S.max:1) {
    model[[S]] <- npivaic(Y,
                          X,
                          W,
                          K.w.degree=3,   
                          K.w.segments=S,
                          J.x.degree=3,
                          J.x.segments=S)
    AIC[S] <- model[[S]]$AIC.c
}

## In this illustration degree is set to 3 for Psi and B and we
## consider knots J, K with K >= J from 1,2,...,10 (here S.max is only
## used in the plot legend)

## J.max <- K.max <- S.max <- 10
## S <- 1
## for(J in 1:J.max) {
##   for(K in J:K.max) {
##     model[[S]] <- npivaic(Y,
##                           X,
##                           W,
##                           K.w.degree=3,   
##                           K.w.segments=K,
##                           J.x.degree=3,
##                           J.x.segments=J)
##     AIC[S] <- model[[S]]$AIC.c
##     S <- S+1
##   }
## }

## Select the model that minimizes AIC (model index will also be the
## value of S used above)

model <- model[[which.min(AIC)]]

## Create a plot of the instrumental regression function, data, and DGP

plot(X,Y,cex=0.25,
     col="lightgrey",
     sub=paste("n = ",format(n,format="d", big.mark=','),
               " (knots searched over: 1,2,...,",S.max,")",sep=""),
     main="Data-Driven Knot Selection (AIC, Hurvich et al (1998))",
     xlab="X",
     ylab="Y")

lines(X[order(X)],h0[order(X)],lty=1,col=1,lwd=1)
lines(X[order(X)],model$fitted[order(X)],lty=2,col=2,lwd=2)

legend("topleft",c("DGP (IV Function h0)",paste("NPIV (K.w.degree = ",model$K.w.degree,
                               ", W.knots = ",model$K.w.segments+1,
                               ", J.x.degree = ", model$J.x.degree,
                               ", X.knots = ",model$J.x.segments+1,")",sep="")),
       lty=1:2,
       col=1:2,
       lwd=c(1,2),
       bty="n",
       cex=0.75)
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
