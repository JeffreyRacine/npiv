\name{Engel95}
\docType{data}
\alias{Engel95}
\title{ 1995 British Family Expenditure Survey  }
\description{
This dataset is based on a sample taken from the British Family Expenditure Survey for 1995. It includes households consisting of married or cohabiting couples with an employed head of household, aged between 25 and 55 years, and with at most two children. There are 1655 household-level observations in total.
}
\usage{data("Engel95")}
\format{ A data frame with 10 columns, and 1655 rows.
  \describe{
 \item{food}{ expenditure share on food, of type \code{numeric}}
 \item{catering}{ expenditure share on catering, of type \code{numeric}}
 \item{alcohol}{ expenditure share on alcohol, of type \code{numeric}}
 \item{fuel}{ expenditure share on fuel, of type \code{numeric}}
 \item{motor}{ expenditure share on motor, of type \code{numeric}}
 \item{fares}{ expenditure share on fares, of type \code{numeric}}
 \item{leisure}{ expenditure share on leisure, of type \code{numeric}}
 \item{logexp}{ logarithm of total expenditure, of type \code{numeric}}
 \item{logwages}{ logarithm of total earnings, of type \code{numeric}}
 \item{nkids}{ `0` indicates no children, `1` indicates 1-2 children, of type \code{numeric}}  
  }
}
\source{ Richard Blundell and Dennis Kristensen }
\references{

Blundell, R., X. Chen and D. Kristensen (2007). \dQuote{Semi-Nonparametric IV Estimation of Shape-Invariant Engel Curves.} \emph{Econometrica}, \strong{75}(6), 1613-1669. doi: \doi{10.1111/j.1468-0262.2007.00808.x}
  
Chen, X. and T. Christensen (2018). \dQuote{Optimal Sup-norm Rates and Uniform Inference on Nonlinear Functionals of Nonparametric IV Regression.} \emph{Quantitative Economics}, \strong{9}(1), 39-85. doi: \doi{10.3982/QE722}
    
Chen, X., T. Christensen and S. Kankanala (2024). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric Structural Functions and Elasticities.} \emph{Review of Economic Studies}, \strong{forthcoming}. doi: \doi{10.1093/restud/rdae025}

}

\keyword{datasets}

\examples{
## Load data
data("Engel95", package = "npiv")

## Sort on logexp (the regressor) for plotting purposes
Engel95 <- Engel95[order(Engel95$logexp),] 
attach(Engel95)
logexp.eval <- seq(4.5,6.5,length=100)

## Estimate the Engel curve for food using logwages as an instrument
food_engel <- npiv(food, logexp, logwages, X.eval = logexp.eval)

## Plot the estimated function and uniform confidence bands
plot(food_engel, showdata = TRUE)
}
