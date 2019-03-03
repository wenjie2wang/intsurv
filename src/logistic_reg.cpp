//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
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

#include "logistic_reg.h"
#include "l-bfgs.h"


// A hopefully fast subroutine fitting logistic model
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_logistic(const arma::mat& x, const arma::vec y)
{
    Intsurv::LogisticReg object {Intsurv::LogisticReg(x, y)};
    arma::vec beta {arma::zeros(x.n_cols)};
    Intsurv::control_line_search control_ls {
        Intsurv::control_line_search()
    };
    Intsurv::control_lbfgs control {Intsurv::control_lbfgs()};
    // TODO: may tweak the parameters based on inputs
    control.line_search = control_ls;
    control.epsilon = 1e-3;
    control.m = 6;
    Intsurv::lbfgs(beta, object, control);
    return Rcpp::NumericVector(beta.begin(), beta.end());
}
