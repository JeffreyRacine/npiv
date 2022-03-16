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

  Blundell, R. and X. Chen and D. Kristensen (2007),
  \dQuote{Semi-Nonparametric IV Estimation of Shape-Invariant Engel
  Curves,} Econometrica, 75, 1613-1669.
    
  Chen, X. and T. Christensen and S. Kankanala (2021), \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric IV}, https://arxiv.org/pdf/2107.11869.pdf

}

\keyword{datasets}

\examples{
data(Engel95)

## Sort on logexp (the endogenous regressor) for plotting purposes

Engel95 <- Engel95[order(Engel95$logexp),] 

attach(Engel95)

model.iv <- npiv(Y=food,
                 X=logexp,
                 W=logwages)

## Compute the non-IV regression (i.e. regress Y on X, simply use
## logexp for W)

model.niv <- npiv(Y=food,
                  X=logexp,
                  W=logexp,
                  K.w.segments=model.iv$K.w.segments,
                  J.x.segments=model.iv$J.x.segments)

## For the plots, restrict focal attention to the bulk of the data
## (i.e. for the plotting area trim out 1/4 of one percent from each
## tail of Y and X)

trim <- 0.0025

plot(logexp,food,
     ylab="Food Budget Share",
     xlab="log(Total Expenditure)",
     xlim=quantile(logexp,c(trim,1-trim)),
     ylim=quantile(food,c(trim,1-trim)),
     main="Nonparametric Instrumental Variables",
     type="p",
     cex=.5,
     col="lightgrey")

lines(logexp,model.iv$fitted,col="blue",lwd=2,lty=2)

lines(logexp,model.niv$fitted,col="red",lwd=2,lty=4)

legend("topright",
       c("Nonparametric IV","Nonparametric Regression"),
       lty=c(2,4),
       col=c("blue","red"),
       lwd=c(2,2),
       bty="n")
}
