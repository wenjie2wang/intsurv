# intsurv

[![CRAN Version][cran-version]][cran]
[![Dependence][tinyverse-dep]][cran]
[![Build Status][gha-icon]][gha-url]
[![Code Coverage][codecov-main]][codecov]

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

If the version under development is able to pass the automated package checks,
one may consider installing it with the help of **remotes** by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/intsurv")
```


## References

- Wang, W., Aseltine, R. H., Chen, K., & Yan, J. (2020). Integrative Survival
  Analysis with Uncertain Event Times in Application to A Suicide Risk
  Study. *Annals of Applied Statistics*, 14(1), 51-73.
- Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen,
  K. (2020). Suicide Risk Modeling with Uncertain Diagnostic Records. *arXiv
  preprint arXiv:2009.02597*.


[cran]: https://cran.r-project.org/package=intsurv
[cran-version]: https://www.r-pkg.org/badges/version/intsurv
[tinyverse-dep]: https://tinyverse.netlify.com/badge/intsurv
[gha-icon]: https://github.com/wenjie2wang/intsurv/workflows/R-CMD-check/badge.svg
[gha-url]: https://github.com/wenjie2wang/intsurv/actions
[codecov]: https://codecov.io/gh/wenjie2wang/intsurv
[codecov-main]: https://codecov.io/gh/wenjie2wang/intsurv/branch/main/graph/badge.svg
