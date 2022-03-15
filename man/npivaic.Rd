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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, X, W, K.w.degree = 3, K.w.segments = 1, J.x.degree = 3, 
    J.x.segments = 1, knots = c("uniform", "quantiles"), basis = c("tensor", 
        "additive", "glp"), check.is.fullrank = FALSE, chol.pivot = FALSE, 
    lambda = sqrt(.Machine$double.eps)) 
{
    basis <- match.arg(basis)
    knots <- match.arg(knots)
    if (missing(Y)) 
        stop(" must provide Y")
    if (missing(X)) 
        stop(" must provide X")
    if (missing(W)) 
        stop(" must provide W")
    if (K.w.degree < 0) 
        stop("K.w.degree must be a non-negative integer")
    if (J.x.degree < 0) 
        stop("J.x.degree must be a non-negative integer")
    if (K.w.segments <= 0) 
        stop("K.w.segments must be a positive integer")
    if (J.x.segments <= 0) 
        stop("J.x.segments must be a positive integer")
    if (check.is.fullrank) {
        if (!is.fullrank(Y)) 
            stop("Y is not of full column rank")
        if (!is.fullrank(X)) 
            stop("X is not of full column rank")
        if (!is.fullrank(W)) 
            stop("W is not of full column rank")
    }
    if (K.w.degree + K.w.segments < J.x.degree + J.x.segments) 
        stop("K.w.degree+K.w.segments must be >= J.x.degree+J.x.segments")
    if (K.w.degree == 0) {
        B.w <- matrix(1, NROW(W), 1)
    }
    else {
        B.w <- prod.spline(x = W, K = cbind(rep(K.w.degree, NCOL(W)), 
            rep(K.w.segments, NCOL(W))), knots = knots, basis = basis)
        if (basis != "tensor") 
            B.w <- cbind(1, B.w)
    }
    if (J.x.degree == 0) {
        Psi.x <- matrix(1, NROW(X), 1)
    }
    else {
        Psi.x <- prod.spline(x = X, K = cbind(rep(J.x.degree, 
            NCOL(X)), rep(J.x.segments, NCOL(X))), knots = knots, 
            basis = basis)
        if (basis != "tensor") 
            Psi.x <- cbind(1, Psi.x)
    }
    Psi.xTB.wB.wTB.w.invB.w <- t(Psi.x) \%*\% B.w \%*\% chol2inv(chol(t(B.w) \%*\% 
        B.w, pivot = chol.pivot)) \%*\% t(B.w)
    h <- Psi.x \%*\% (chol2inv(chol(Psi.xTB.wB.wTB.w.invB.w \%*\% 
        Psi.x + diag(lambda, NCOL(Psi.x)), pivot = chol.pivot)) \%*\% 
        Psi.xTB.wB.wTB.w.invB.w \%*\% Y)
    trH <- dim(Psi.x)[2]
    aic.penalty <- (1 + trH/length(Y))/(1 - (trH + 2)/length(Y))
    aic.c <- ifelse(aic.penalty > 0, log(mean((Y - h)^2)) + aic.penalty, 
        .Machine$double.xmax)
    return(list(fitted = h, K.w.degree = K.w.degree, K.w.segments = K.w.segments, 
        J.x.degree = J.x.degree, J.x.segments = J.x.segments, 
        AIC.c = aic.c))
  }
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
