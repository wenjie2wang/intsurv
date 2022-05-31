//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2022  Wenjie Wang <wang@wwenjie.org>
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
#include <intsurv/CrossValidation.h>
#include <intsurv/utils.h>

// [[Rcpp::export]]
Rcpp::List rcpp_gen_cv_index(const unsigned int nobs,
                             const unsigned int nfolds,
                             const arma::uvec& strata,
                             const arma::uvec& static_train_index)
{
    intsurv::CrossValidation cv_obj {
        nobs, nfolds, static_train_index, strata
    };
    Rcpp::List train_list, valid_list;
    for (size_t i {0}; i < nfolds; ++i) {
        arma::uvec tmp1 { cv_obj.train_index_.at(i) + 1 };
        train_list.push_back(intsurv::arma2rvec(tmp1));
        arma::uvec tmp2 { cv_obj.test_index_.at(i) + 1 };
        valid_list.push_back(intsurv::arma2rvec(tmp2));
    }
    return Rcpp::List::create(
        Rcpp::Named("train_index") = train_list,
        Rcpp::Named("valid_index") = valid_list
        );
}
