# network-model-R-package

[![Build Status](https://www.travis-ci.org/AIforGoodSimulator/network-model-R-package.svg?branch=master)](https://www.travis-ci.org/AIforGoodSimulator/network-model-R-package)

This package has been made to make network epidemic models of COVID19 spread in refugee camps.

To install run:
```R
devtools::install_github('AIforGoodSimulator/network-model-R-package')
```

EpiModel has recently received an update and we have yet to test
the compatibility of our code with the newest version. The package works
using the packages in renv.lock. To set that up run:
```R
renv::restore()
```
from within the repo root.

Read more about renv [here](https://rstudio.github.io/renv/)
