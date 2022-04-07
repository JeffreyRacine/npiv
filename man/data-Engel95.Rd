\name{Engel95}
\docType{data}
\alias{Engel95}
\title{ 1995 British Family Expenditure Survey  }
\description{
British cross-section data consisting of a random sample taken from
the British Family Expenditure Survey for 1995. The households consist
of married couples with an employed head-of-household between the ages
of 25 and 55 years. There are 1655 household-level observations in
total.
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
 \item{nkids}{ number of children, of type \code{numeric}}  
  }
}
\source{ Richard Blundell and Dennis Kristensen }
\references{

Blundell, R., X. Chen and D. Kristensen (2007). \dQuote{Semi-Nonparametric IV Estimation of Shape-Invariant Engel Curves.} \emph{Econometrica}, \strong{75}(6), 1613-1669. doi: \url{https://doi.org/10.1111/j.1468-0262.2007.00808.x}
  
Chen, X. and T. Christensen (2018). \dQuote{Optimal Sup-norm Rates and Uniform Inference on Nonlinear Functionals of Nonparametric IV Regression.} \emph{Quantitative Economics}, \strong{9}(1), 39-85. doi: \url{https://doi.org/10.3982/QE722}
    
Chen, X., T. Christensen and S. Kankanala (2021). \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV.} \url{https://arxiv.org/pdf/2107.11869.pdf}

}

\keyword{datasets}

\examples{
## load data
data("Engel95", package = "npiv")

## sort on logexp (the regressor) for plotting purposes
Engel95 <- Engel95[order(Engel95$logexp),] 
attach(Engel95)

## Estimate the Engel curve for food using logwages as an instrument
fm1 <- npiv(food ~ logexp | logwages)

## Plot the estimated Engel curve and data-driven uniform confidence bands
plot(logexp,food,
     ylab="Food Budget Share",
     xlab="log(Total Household Expenditure)",
     xlim=c(4.75, 6.25),
     ylim=c(0, 0.4),
     main="",
     type="p",
     cex=.5,
     col="lightgrey")
lines(logexp,fm1$h,col="blue",lwd=2,lty=1)
lines(logexp,fm1$h.upper,col="blue",lwd=2,lty=2)
lines(logexp,fm1$h.lower,col="blue",lwd=2,lty=2)

legend("topright",
       legend=c("Estimate",
                "95 per cent UCBs"),
       col=c("blue","blue"),
       lty=c(1,2),
       lwd=c(2,2))
}
