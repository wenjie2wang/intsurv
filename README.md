# intsurv

![CRAN_Status_Badge][cranVersion]
[![Build Status][travis_master]][travis]


The R package **intsurv** mainly provides function fitting extended Cox
proportional hazard model for uncertain survival data due to imperfect data
integration proposed by Wang (2017+).


## Development

The package is still under active development.

If it is able to pass the building check by Travis CI, you may consider
installing it with the help of **remotes** by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/intsurv")
```


## Reference

Wang, W., Chen, K., & Yan, J. (2017+).  Extended Cox Model by ECM Algorithm for
Uncertain Survival Records Due to Imperfect Data Integration. (working in
progress)


[cranVersion]: http://www.r-pkg.org/badges/version/intsurv
[travis]: https://travis-ci.org/wenjie2wang/intsurv
[travis_master]: https://travis-ci.org/wenjie2wang/intsurv.svg?branch=master
