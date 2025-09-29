# intsurv

[![CRAN Version][cran-version]][cran]
[![Build Status][gha-icon]][gha-url]
[![Code Coverage][codecov-main]][codecov]

The {intsurv} R package contains implementations of

- integrative Cox model with uncertain event times (Wang et al., 2020)
- regularized Cox cure rate model with uncertain event status (Wang et al., 2023)

and other survival analysis routines, including

- Cox cure rate model
- regularized Cox cure rate model with elastic net penalty

## Installation

Install from CRAN by running:

```R
install.packages("intsurv")
```

## Usage

For examples and detailed documentation, check:

- CRAN reference manual: https://cran.r-project.org/web/packages/intsurv/refman/intsurv.html
- Package site: https://wwenjie.org/intsurv/

## References

- Wang, W., Aseltine, R. H., Chen, K., & Yan, J. (2020). Integrative Survival
  Analysis with Uncertain Event Times in Application to A Suicide Risk
  Study. *Annals of Applied Statistics*, 14(1), 51--73.
- Wang, W., Luo, C., Aseltine, R. H., Wang, F., Yan, J., & Chen,
  K. (2023). Survival Modeling of Suicide Risk with Rare and Uncertain
  Diagnoses. Statistics in Biosciences, 17(1), 35--61.

[cran]: https://cran.r-project.org/package=intsurv
[cran-version]: https://www.r-pkg.org/badges/version/intsurv
[gha-icon]: https://github.com/wenjie2wang/intsurv/actions/workflows/R-CMD-check.yaml/badge.svg
[gha-url]: https://github.com/wenjie2wang/intsurv/actions/workflows/R-CMD-check.yaml
[codecov]: https://app.codecov.io/gh/wenjie2wang/intsurv
[codecov-main]: https://codecov.io/gh/wenjie2wang/intsurv/branch/main/graph/badge.svg
