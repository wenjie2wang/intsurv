##
## intsurv: Integrative Survival Models
## Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
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

## is `x` object of class `Survi`?
is.Survi <- function(x)
{
    inherits(x, "Survi")
}

## remove NA's from vector `x`
rmNA <- function(x)
{
    x[! is.na(x)]
}
