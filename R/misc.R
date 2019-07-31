##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
##
## This file is part of the R package intsurv.
##
## The R package intsurv is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package intsurv is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

## wrap messages and keep proper line length
wrapMessages <- function(..., strwrap.args = list())
{
    x <- paste(...)
    wrap_x <- do.call(strwrap, c(list(x = x), strwrap.args))
    paste(wrap_x, collapse = "\n")
}

## is `x` object of class `foo`?
is_Survi <- function(x)
{
    inherits(x, "Survi")
}
is_iCoxph <- function(x)
{
    inherits(x, "iCoxph")
}
is_iCoxph.control <- function(x)
{
    inherits(x, "iCoxph.control")
}
is_iCoxph.start <- function(x)
{
    inherits(x, "iCoxph.start")
}

## remove NA's from vector `x`
rmNA <- function(x)
{
    x[! is.na(x)]
}

## computing L2-norm of vector x
L2norm <- function(x) {
    sqrt(sum(x ^ 2))
}
L2norm2 <- function(x) {
    sum(x ^ 2)
}

se_interQ <- function(x) {
    diff(stats::quantile(x, probs = c(0.25, 0.75))) /
        (stats::qnorm(0.75) - stats::qnorm(0.25))
}

## throw warnings if `...` is specified by mistake
warn_dots <- function(...) {
    dotsList <- list(...)
    .fun_name <- as.character(sys.call(- 1L)[[1L]])
    if (length(dotsList) > 0) {
        list_names <- names(dotsList)
        if (is.null(list_names)) {
            warning(
                sprintf("Some argument(s) went into `...` of %s()",
                        .fun_name),
                call. = FALSE
            )
        } else {
            list_names <- list_names[list_names != ""]
            if (length(list_names) > 2) {
                all_names <- paste(sprintf("'%s'", list_names), collapse = ", ")
                all_names <- gsub("(.+), (.+)$", "\\1, and \\2", all_names)
            } else {
                all_names <- paste(sprintf("'%s'", list_names),
                                   collapse = " and ")
            }
            warning(
                sprintf("The argument %s went into `...` of %s().",
                        all_names, .fun_name),
                call. = FALSE
            )
        }
    }
    invisible(NULL)
}
