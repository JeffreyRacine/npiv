## Empirical illustration: Engel curve estimation

library(npiv)

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

## Repeat for leisure
leisure_engel <- npiv(leisure, logexp, logwages, X.eval = logexp.eval)
plot(leisure_engel, showdata = TRUE)

## Estimate the Engel curve for food by regression using the default (knots uniformly spaced)
food_engel_reg_uniform <- npiv(food, logexp, logexp, X.eval = logexp.eval)
plot(food_engel_reg_uniform, showdata = TRUE)

## Repeat with knots at quantiles
food_engel_reg_quantiles <- npiv(food, logexp, logexp, X.eval = logexp.eval, knots = "quantiles")
plot(food_engel_reg_quantiles, showdata = TRUE)
