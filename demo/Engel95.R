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


## Compare with the Engel curve for leisure
fm2 <- npiv(leisure ~ logexp | logwages)

## Plot the estimated Engel curve and data-driven uniform confidence bands
plot(logexp,leisure,
     ylab="Leisure Budget Share",
     xlab="log(Total Household Expenditure)",
     xlim=c(4.75, 6.25),
     ylim=c(0, 0.4),
     main="",
     type="p",
     cex=.5,
     col="lightgrey")
lines(logexp,fm2$h,col="blue",lwd=2,lty=1)
lines(logexp,fm2$h.upper,col="blue",lwd=2,lty=2)
lines(logexp,fm2$h.lower,col="blue",lwd=2,lty=2)

legend("topright",
       legend=c("Estimate",
                "95 per cent UCBs"),
       col=c("blue","blue"),
       lty=c(1,2),
       lwd=c(2,2))
