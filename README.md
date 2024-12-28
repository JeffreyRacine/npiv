# npiv

This is the R package `npiv` (Nonparametric Instrumental Variables Estimation and Inference) written by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored and maintained by Timothy Christensen (timothy.christensen@yale.edu).

## Description

This package implements methods introduced in [Chen, Christensen, and Kankanala (2024)](https://doi.org/10.1093/restud/rdae025) for estimating and constructing uniform confidence bands for nonparametric structural functions using instrumental variables, including data-driven choice of tuning parameters. It also provides functionality to construct uniform confidence bands using the method of [Chen and Christensen (2018)](https://doi.org/10.3982/QE722). All methods in this package apply to nonparametric regression as a special case.

## Installation

You can install by either downloading the [zip ball](https://github.com/JeffreyRacine/npiv/zipball/main) or [tar ball](https://github.com/JeffreyRacine/npiv/tarball/main), decompress and run `R CMD INSTALL` on it.

Alternatively, you can install the development version but before doing so Windows users have to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), while OS X users have to first install [Xcode](https://apps.apple.com/us/app/xcode/id497799835) and the command line tools (in OS X 10.9 or higher, once you have Xcode installed, open a terminal and run xcode-select --install). Note also that versions of e.g. Rtools are paired with versions of R so ensure you have the latest version of R installed prior to commencing this process.

After installing [Rtools](https://cran.r-project.org/bin/windows/Rtools/)/[Xcode](https://apps.apple.com/us/app/xcode/id497799835) and **devtools** (via install.packages("devtools")), install the development package using the following command:

```r
library(devtools); install_github('JeffreyRacine/npiv')
```
