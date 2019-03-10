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

#include "logistic_reg.h"

// fitting regular logistic model by monotonic quadratic approximation algorithm
// non-integer y vector is allowed
// reference: Bohning and Lindsay (1988) SIAM
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_logistic(const arma::mat& x,
                                  const arma::vec& y,
                                  const bool intercept = true,
                                  const arma::vec start = 0,
                                  const unsigned int max_iter = 1000,
                                  const double rel_tol = 1e-6)
{
    Intsurv::LogisticReg object { Intsurv::LogisticReg(x, y, intercept) };
    object.fit(start, max_iter, rel_tol);
    arma::vec beta {object.coef};
    return Rcpp::NumericVector(beta.begin(), beta.end());
}


// regularized logistic model by coordinate-majorization-descent algorithm
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_reg_logistic(const arma::mat& x,
                                      const arma::vec& y,
                                      const double lambda = 0,
                                      const arma::vec& penalty_factor = 0,
                                      const bool intercept = true,
                                      const arma::vec start = 0,
                                      const unsigned int max_iter = 1000,
                                      const double rel_tol = 1e-6)
{
    Intsurv::LogisticReg object { Intsurv::LogisticReg(x, y, intercept) };
    object.regularized_fit(lambda, penalty_factor,
                           start, max_iter, rel_tol);
    arma::vec beta { object.coef };
    return Rcpp::NumericVector(beta.begin(), beta.end());
}
