library(tinytest)
Rcpp::sourceCpp("test-CrossValidation.cpp")

set.seed(123)
nobs <- 300
nfolds <- 5
strata <- c(rep(0, nobs / 3), rep(1, nobs / 3), rep(2, nobs / 3))
strata <- sample(strata)
static_train_index <- c(99, 100)

train_size0 <- nobs / nfolds * (nfolds - 1)
train_size1 <- train_size0 - length(static_train_index)
train_size2 <- train_size0 + length(static_train_index)

res <- rt_cv(nobs, nfolds, static_train_index, strata)

## on static train index
expect_true(all(
    sapply(res$train_index, function(a) static_train_index %in% a)
))
expect_false(any(
    sapply(res$valid_index, function(a) static_train_index %in% a)
))

train_len0 <- sapply(res$train_index, function(a)
    length(setdiff(a, static_train_index)))
expect_true(all(train_len0 >= train_size1))
expect_true(all(train_len0 <= train_size2))

## on strata
strata_size0 <- train_size0 / 3
strata_max_size <- strata_size0 + length(static_train_index)
strata_min_size <- strata_size0 - length(static_train_index)
strata_size_mat <- sapply(res$train_index,
                          function(a) range(table(strata[a + 1])))
expect_true(all(
    strata_size_mat[1, ] >= strata_min_size
))
expect_true(all(
    strata_size_mat[2, ] <= strata_max_size
))
