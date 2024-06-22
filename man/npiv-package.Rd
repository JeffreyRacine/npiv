\name{npiv-package}
\alias{npiv-package}
\docType{package}
\title{
\packageTitle{npiv}
}
\description{
This package provides tools for data-driven statistical analysis using the nonparametric instrumental variables methods described in Chen, Christensen and Kankanala (2024) and Chen and Christensen (2018): \code{npiv} performs performs nonparametric two-stage least-squares estimation of, and constructs uniform confidence bands for, a structural function \code{h0} and its derivatives using a B-spline sieve, and \code{npiv_choose_J} performs data-driven choice of sieve dimension.
}
\details{
This package provides a function with a simple interface for conducting nonparametric instrumental variable (IV) estimation and inference, \code{npiv(...)} given a dependent variable vector \code{Y}, matrix of endogenous predictors \code{X}, and matrix of instruments \code{W}. All methods in this package apply to nonparametric regression as a special case when \code{X} is exogenous.

The function \code{npiv} performs nonparametric two-stage least-squares estimation of, and constructs uniform confidence bands for, a structural function \code{h0} and its derivatives. Sieve dimensions are determined in a data-dependent way if not provided by the user, via the function \code{npiv_choose_J} which implements the methods described in Chen, Christensen, and Kankanala (2024). This data-driven choice of sieve dimension ensures estimators of \code{h0} and its derivatives converge at the best possible (i.e., minimax) sup-norm rate. The resulting uniform confidence bands for \code{h0} and its derivatives contract within a logarithmic factor of the minimax rate.

If sieve dimensions are provided by the user, \code{npiv} implements the bootstrap-based procedure of Chen and Christensen (2018) to construct uniform confidence bands for \code{h0} and its derivatives. 


}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy M. Christensen <timothy.christensen@nyu.edu>

Maintainer: Jeffrey S. Racine <racinej@mcmaster.ca>
}
\references{

Chen, X. and T. Christensen (2018). \dQuote{Optimal Sup-norm Rates and Uniform Inference on Nonlinear Functionals of Nonparametric IV Regression.} \emph{Quantitative Economics}, \strong{9}(1), 39-85. doi: \url{https://doi.org/10.3982/QE722}

Chen, X., T. Christensen and S. Kankanala (2024). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric Structural Functions and Elasticities.} \emph{Review of Economic Studies}, \strong{forthcoming}. doi: \url{https://doi.org/10.1093/restud/rdae025}

}
\keyword{ package }
