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
#include <intsurv.h>

// Nelson-Aalen estimator for right censor data
// [[Rcpp::export]]
Rcpp::List rcpp_mcf_right(const arma::vec& time,
                          const arma::vec& event)
{
    Intsurv::NelsonAalen na_obj { time, event };
    return Rcpp::List::create(
        Rcpp::Named("time") = Intsurv::arma2rvec(na_obj.uni_event_time),
        Rcpp::Named("inst_rate") = Intsurv::arma2rvec(na_obj.inst_rate),
        Rcpp::Named("cum_rate") = Intsurv::arma2rvec(na_obj.cum_rate)
        );
}
