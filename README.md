# intsurv

![CRAN_Status_Badge][cranVersion]
[![Build Status][travis_master]][travis]


The R package **intsurv** mainly provides functions fitting integrative Cox
proportional hazard model proposed by Wang (2017+) for uncertain survival data
due to imperfect data integration.


## Development

The package is still under active development.

If it is able to pass the building check by Travis CI, you may consider
installing it with the help of **remotes** by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/intsurv")
```


## Reference

Wang, W., Aseltine, R., Chen, K., & Yan, J. (2017+).  Integrative Survival
Analysis with Uncertain Event Records from Imperfect Data Integration. (working
in progress)


[cranVersion]: http://www.r-pkg.org/badges/version/intsurv
[travis]: https://travis-ci.org/wenjie2wang/intsurv
[travis_master]: https://travis-ci.org/wenjie2wang/intsurv.svg?branch=master
