
# npiv

This is the R package **npiv** (Nonparametric Instrumental Variables
Estimation and Inference) written by Jeffrey S. Racine
(<racinej@mcmaster.ca>) and co-authored and maintained by Timothy
Christensen (<timothy.christensen@yale.edu>).

## Description

This package implements methods introduced in [Chen, Christensen, and
Kankanala (2024)](https://doi.org/10.1093/restud/rdae025) for estimating
and constructing uniform confidence bands for nonparametric structural
functions using instrumental variables, including data-driven choice of
tuning parameters. It also provides functionality to construct uniform
confidence bands using the method of [Chen and Christensen
(2018)](https://doi.org/10.3982/QE722). All methods in this package
apply to nonparametric regression as a special case.

## Installation

You can install by either downloading the [zip
ball](https://github.com/JeffreyRacine/npiv/zipball/main) or [tar
ball](https://github.com/JeffreyRacine/npiv/tarball/main), decompress
and run `R CMD INSTALL` on it.

Alternatively, you can install the development version but before doing
so Windows users have to first install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), while OS X
users have to first install
[Xcode](https://apps.apple.com/us/app/xcode/id497799835) and the command
line tools (in OS X 10.9 or higher, once you have Xcode installed, open
a terminal and run xcode-select –install). Note also that versions of
e.g. Rtools are paired with versions of R so ensure you have the latest
version of R installed prior to commencing this process.

After installing
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)/[Xcode](https://apps.apple.com/us/app/xcode/id497799835)
and **devtools** (via `install.packages("devtools")`), install the
development package using the following command:

``` r
library(devtools); install_github('JeffreyRacine/npiv')
```

## Examples

### Nonparametric instrumental variables estimation

We present a simple example to Engel curve estimation using data from
the British Family Expenditure Survey. We first load the data and sort
on log household expenditure for plotting purposes. We generate a grid
of log expenditure from 4.5 to 6.5 for plotting.

``` r
data("Engel95", package = "npiv")
Engel95 <- Engel95[order(Engel95$logexp),] 
attach(Engel95)
logexp.eval <- seq(4.5,6.5,length=100)
```

To nonparametrically regress `Y` on `X` using `W` as an instrument, use
`npiv(Y ~ X | W)` or `npiv(Y, X, W)`. By default, `npiv` uses a
data-driven choice of sieve dimension based on the method of [Chen,
Christensen, and Kankanala
(2024)](https://doi.org/10.1093/restud/rdae025). We include
`X.eval = logexp.eval` since we want to plot over a smaller region than
the support of log expenditure.

``` r
food_engel <- npiv(food, logexp, logwages, X.eval = logexp.eval)
```

We plot the estimated function and 95% uniform confidence bands using
the command `plot`. By default, the confidence bands are generated using
the method of [Chen, Christensen, and Kankanala
(2024)](https://doi.org/10.1093/restud/rdae025). We include
`showdata = TRUE` to overlay the data points.

``` r
plot(food_engel, showdata = TRUE)
```

![](man/figures/README-unnamed-chunk-5-1.png)<!-- -->

Confidence bands for the derivative of the structural function are
plotted by including the argument `type = "deriv"`.

``` r
plot(food_engel, type = "deriv")
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

### Nonparametric regression

We can estimate a conditional mean function by nonparametric regression
with data-driven choice of sieve dimension simply by passing the
regressor as the instrument. The following example estimates the Engel
curve for food by regression of `food` on `logexp`:

``` r
food_engel_reg_uniform <- npiv(food, logexp, logexp, X.eval = logexp.eval)
plot(food_engel_reg_uniform, showdata = TRUE)
```

![](man/figures/README-unnamed-chunk-7-1.png)<!-- -->

The plot is wiggly, indicating that the algorithm has selected a fairly
large sieve dimension. This can sometimes happen when the data are far
from uniform, because the default method uses splines with knots placed
uniformly over the support of the data. We can restore good performance
by including `knots = "quantiles"` to use splines with knots placed at
the quantiles of the data.

``` r
food_engel_reg_quantiles <- npiv(food, logexp, logexp, X.eval = logexp.eval, knots = "quantiles")
plot(food_engel_reg_quantiles, showdata = TRUE)
```

![](man/figures/README-unnamed-chunk-8-1.png)<!-- -->
