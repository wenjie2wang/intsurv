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

## throw warnings if `...` is specified by mistake
warn_dots <- function(..., .fun_name = NULL) {
    dotsList <- list(...)
    if (length(dotsList) > 0) {
        warning(
            "Some arguments went into `...`",
            if (! is.null(.fun_name)) sprintf(" of %s()", .fun_name),
            ", which is however not used currently.",
            call. = FALSE
        )
    }
    invisible(NULL)
}
