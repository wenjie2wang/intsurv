//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
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

#include "coxph_reg.h"
#include "utils.h"

// fitting Cox model by monotonic quadratic approximation algorithm
// that allows non-integer "event" and tied events
// reference: Bohning and Lindsay (1988) SIAM
// [[Rcpp::export]]
arma::vec rcpp_coxph(const arma::vec& time,
                     const arma::vec& event,
                     const arma::mat& x,
                     const arma::vec& offset = 0,
                     const arma::vec& start = 0,
                     const unsigned int max_iter = 1000,
                     const double rel_tol = 1e-6)
{
    Intsurv::CoxphReg object { Intsurv::CoxphReg(time, event, x) };
    object.set_offset(offset);
    object.fit(start, max_iter, rel_tol);
    return object.coef;
}


// fitting regularized Cox model with coordinate-majorizatio-descent algorithm
// that allows non-integer "event" and tied events
// [[Rcpp::export]]
arma::vec rcpp_reg_coxph(const arma::vec& time,
                         const arma::vec& event,
                         const arma::mat& x,
                         const double lambda = 0,
                         arma::vec penalty_factor = 0,
                         const arma::vec& offset = 0,
                         const arma::vec& start = 0,
                         const unsigned int max_iter = 1000,
                         const double rel_tol = 1e-6)
{
    Intsurv::CoxphReg object { Intsurv::CoxphReg(time, event, x) };
    object.set_offset(offset);
    object.regularized_fit(lambda, penalty_factor, start, max_iter, rel_tol);
    return object.coef;
}
