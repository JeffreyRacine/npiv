# npiv

This is the R package `npiv' (Nonparametric Instrumental Variables Estimation) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored by Timothy M. Christensen (timothy.christensen@nyu.edu)

## Installation

You can install by either downloading the [zip
 ball](https://github.com/JeffreyRacine/npiv/zipball/main)
 or [tar
 ball](https://github.com/JeffreyRacine/npiv/tarball/main),
 decompress and run `R CMD INSTALL` on it.

Alternatively, you can install the development version but before
doing so Windows users have to first install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), while OS X
users have to first install
[Xcode](https://apps.apple.com/us/app/xcode/id497799835) and the
command line tools (in OS X 10.9 or higher, once you have Xcode
installed, open a terminal and run xcode-select --install). Note also
that versions of e.g. Rtools are paired with versions of R so ensure
you have the latest version of R installed prior to commencing this
process.

After installing
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)/[Xcode](https://apps.apple.com/us/app/xcode/id497799835)
and **devtools** (via install.packages("devtools")), install the
development package using the following command:

```r
library(devtools); install_github('JeffreyRacine/npiv')
```

