\name{npiv-package}
\alias{npiv-package}
\docType{package}
\title{
\packageTitle{npiv}
}
\description{
This package implements the nonparametric instrumental variables estimation and inference methods described in Chen, Christensen, and Kankanala (2024) and Chen and Christensen (2018). The function \code{npiv} estimates the nonparametric structural function \code{h0} using B-splines and constructs uniform confidence bands for \code{h0}. The function \code{npiv_choose_J} performs data-driven choice of sieve dimension. All methods in this package apply to estimation and inference for nonparametric regression as a special case.}
\details{
This package provides a function \code{npiv(...)} with a simple interface for performing nonparametric instrumental variable estimation and inference. 

Given a dependent variable vector \code{Y}, matrix of endogenous regressors \code{X}, and matrix of instruments \code{W}, \code{npiv} nonparametrically estimates the structural function \code{h0} and its derivative using B-splines. \code{npiv} can also be used for estimting the conditional mean \code{h0} of \code{Y} given \code{X}, as as well as the derivative of the conditional mean function, by nonparametric regression. 

The function \code{npiv} also constructs uniform confidence bands for \code{h0} and its derivative.

Sieve dimensions are determined in a data-dependent way if not provided by the user via the function \code{npiv_choose_J}, which implements the methods described in Chen, Christensen, and Kankanala (2024). This data-driven choice of sieve dimension ensures estimators of \code{h0} and its derivatives converge at the optimal sup-norm rate. The resulting uniform confidence bands for \code{h0} and its derivative contract within a logarithmic factor of the optimal rate. In this way, \code{npiv} facilitates fully data-driven estimation and uniform inference on \code{h0} and its derivative.

If sieve dimensions are provided by the user, \code{npiv} implements the bootstrap-based procedure of Chen and Christensen (2018) to construct uniform confidence bands for \code{h0} and its derivative. 
}
\author{
Jeffrey S. Racine <racinej@mcmaster.ca>, Timothy Christensen <timothy.christensen@yale.edu>

Maintainer: Timothy Christensen <timothy.christensen@yale.edu>
}
\references{

Chen, X. and T. Christensen (2018). \dQuote{Optimal Sup-norm Rates and Uniform Inference on Nonlinear Functionals of Nonparametric IV Regression.} \emph{Quantitative Economics}, \strong{9}(1), 39-85. doi: \url{https://doi.org/10.3982/QE722}

Chen, X., T. Christensen and S. Kankanala (2024). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric Structural Functions and Elasticities.} \emph{Review of Economic Studies}, \strong{forthcoming}. doi: \url{https://doi.org/10.1093/restud/rdae025}

}
\keyword{ package }
