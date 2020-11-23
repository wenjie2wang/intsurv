//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2020  Wenjie Wang <wang@wwenjie.org>
//
// This file is part of the R package intsurv.
//
// The R package intsurv is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package intsurv is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <intsurv.h>

// [[Rcpp::export]]
Rcpp::NumericVector aggregateSum(const arma::vec& x,
                                 const arma::vec& indices,
                                 const bool simplify = true,
                                 const bool cumulative = false,
                                 const bool reversely = false)
{
    arma::vec res {
        Intsurv::aggregate_sum(x, indices, simplify, cumulative, reversely)
    };
    return Rcpp::NumericVector(res.begin(), res.end());
}

// [[Rcpp::export]]
Rcpp::NumericVector revcumsum(const arma::vec& x)
{
    arma::vec res { Intsurv::cum_sum(x, true) };
    return Rcpp::NumericVector(res.begin(), res.end());
}
