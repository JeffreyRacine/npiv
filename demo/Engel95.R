data(Engel95)

## Sort on logexp (the endogenous regressor) for plotting purposes

Engel95 <- Engel95[order(Engel95$logexp),] 

attach(Engel95)

model.iv <- npiv(food,
                 logexp,
                 logwages)

model.niv <- npiv(food,
                  logexp,
                  logexp)

## For the plots, restrict focal attention to the bulk of the data
## (i.e. for the plotting area trim out 1/4 of one percent from each
## tail of Y and X)

trim <- 0.025

plot(logexp,food,
     ylab="Food Budget Share",
     xlab="log(Total Expenditure)",
     xlim=quantile(logexp,c(trim,1-trim)),
     ylim=quantile(food,c(trim,1-trim)),
     main="Nonparametric Instrumental Variables",
     type="p",
     cex=.5,
     col="lightgrey")

lines(logexp,model.iv$h,col="blue",lwd=2,lty=2)

lines(logexp,model.niv$h,col="red",lwd=2,lty=4)

legend("topright",
       c("Nonparametric IV","Nonparametric Regression"),
       lty=c(2,4),
       col=c("blue","red"),
       lwd=c(2,2),
       bty="n")
