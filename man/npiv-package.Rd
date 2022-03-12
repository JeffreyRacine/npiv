\name{npiv-package}
\alias{npiv-package}
\docType{package}
\title{
\packageTitle{npiv}
}
\description{
This package implements the nonparametric instrumental variable methods described in Chen, X. and T. Christensen and S. Kankanala (2021).
}
\details{
This package provides a function with a simple interface for conducting nonparametric instrumental variable (IV) estimation, \code{npiv(...)} given a dependent variable vector \code{Y}, matrix of endogenous predictors \code{X}, and matrix of instruments \code{W}.

The procedures contained in this package constitute a computationally simple, data-driven procedures for estimation and inference on a structural function \code{h0} and its derivatives in nonparametric models using instrumental variables. The first procedure is a bootstrap-based, data-driven choice of sieve dimension for sieve nonparametric instrumental variables (NPIV) estimators. When implemented with this data-driven choice, sieve NPIV estimators of \code{h0} and its derivatives are adaptive: they converge at the best possible (i.e., minimax) sup-norm rate, without having to know the smoothness of \code{h0}, degree of endogeneity of the regressors, or instrument strength. The second procedure is a data-driven approach for constructing honest and adaptive uniform confidence bands (UCBs) for \code{h0} and its derivatives. The data-driven UCBs guarantee coverage for \code{h0} and its derivatives uniformly over a generic class of data-generating processes (honesty) and contract at, or within a logarithmic factor of, the minimax sup-norm rate (adaptivity). As such, these data-driven UCBs deliver asymptotic efficiency gains relative to UCBs constructed via the usual approach of undersmoothing. In addition, both procedures apply to nonparametric regression as a special case. 

}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M Christensen <timothy.christensen@nyu.edu>

Maintainer: Jeffrey S. Racine <racinej@mcmaster.ca>
}
\references{

Chen, X. and T. Christensen and S. Kankanala (2021), \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV}, https://arxiv.org/pdf/2107.11869.pdf

}
\keyword{ package }
