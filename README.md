# intsurv

![CRAN_Status_Badge][cranVersion]
[![Build Status][travis_master]][travis]


The R package **intsurv** contains implementations of

- integrative Cox model with uncertain event times (Wang et al., 2019)
- Cox cure rate model with uncertain event status (Wang et al., 2019+)

and other survival analysis routines, including

- regular Cox cure rate model
- regularized Cox cure rate model with elastic net penalty
- weighted concordance index


## Installation

You may install the latest released version on CRAN by

```R
install.packages("intsurv")
```


## Get Started

Examples are provided for the main functions for model-fitting in the package.
One may get started from those examples and the function documentation.

```R
library(intsurv)
?iCoxph # integrative Cox model
?cox_cure # Cox cure rate model
?cox_cure_net # regularized Cox cure rate model
```


## Development

If the version under development is able to pass the building check by Travis
CI, you may consider installing it with the help of **remotes** by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/intsurv")
```


## Reference

Wang, W., Aseltine, R., Chen, K., & Yan, J. (2019).  Integrative Survival
Analysis with Uncertain Event Times in Application to a Suicide Risk
Study. *Annals of Applied Statistics*. (in press)

Wang, W., Chen, K., Luo, C., & Yan, J. (2019+). Cox Cure Model with Uncertain
Event Status with application to a Suicide Risk Study. *Working in Progress*.


[cranVersion]: http://www.r-pkg.org/badges/version/intsurv
[travis]: https://travis-ci.org/wenjie2wang/intsurv
[travis_master]: https://travis-ci.org/wenjie2wang/intsurv.svg?branch=master
