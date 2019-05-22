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


// fit regular Cox model that allows non-integer "event" and tied events
// [[Rcpp::export]]
Rcpp::List rcpp_coxph(const arma::vec& time,
                      const arma::vec& event,
                      const arma::mat& x,
                      const arma::vec& start = 0,
                      const unsigned int max_iter = 300,
                      const double rel_tol = 1e-5)
{
    Intsurv::CoxphReg object { Intsurv::CoxphReg(time, event, x) };
    object.fit(start, max_iter, rel_tol);
    object.compute_haz_surv_time();
    object.compute_censor_haz_surv_time();
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef),
        Rcpp::Named("negLogL") = object.negLogL,
        Rcpp::Named("nObs") = object.nObs,
        Rcpp::Named("unique_time") = Intsurv::arma2rvec(object.unique_time),
        Rcpp::Named("h0") = Intsurv::arma2rvec(object.h0_est),
        Rcpp::Named("H0") = Intsurv::arma2rvec(object.H0_est),
        Rcpp::Named("S0") = Intsurv::arma2rvec(object.S0_est),
        Rcpp::Named("hc") = Intsurv::arma2rvec(object.hc_est),
        Rcpp::Named("Hc") = Intsurv::arma2rvec(object.Hc_est),
        Rcpp::Named("Sc") = Intsurv::arma2rvec(object.Sc_est)
        );
}


// fit regularized Cox model with coordinate-majorizatio-descent algorithm
// that allows non-integer "event" and tied events
// for particular lambda's
// [[Rcpp::export]]
Rcpp::List rcpp_reg_coxph1(const arma::vec& time,
                           const arma::vec& event,
                           const arma::mat& x,
                           const double l1_lambda = 0,
                           const double l2_lambda = 0,
                           arma::vec l1_penalty_factor = 0,
                           const arma::vec& start = 0,
                           const unsigned int max_iter = 1000,
                           const double rel_tol = 1e-6)
{
    Intsurv::CoxphReg object { Intsurv::CoxphReg(time, event, x) };
    object.regularized_fit(l1_lambda, l2_lambda, l1_penalty_factor,
                           start, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef,
        Rcpp::Named("en_coef") = Intsurv::arma2rvec(object.en_coef),
        Rcpp::Named("negLogL") = object.negLogL,
        Rcpp::Named("nObs") = object.nObs,
        Rcpp::Named("l1_lambda_max") = object.l1_lambda_max,
        Rcpp::Named("l1_lambda") = object.l1_lambda,
        Rcpp::Named("l2_lambda") = object.l2_lambda,
        Rcpp::Named("l1_penalty_factor") = object.l1_penalty_factor
        );
}


// for a sequence of lambda's and a given alpha
// [[Rcpp::export]]
Rcpp::List rcpp_reg_coxph2(const arma::vec& time,
                           const arma::vec& event,
                           const arma::mat& x,
                           arma::vec lambda = 0,
                           const double alpha = 1,
                           const unsigned int& nlambda = 1,
                           double lambda_min_ratio = 1e-4,
                           arma::vec l1_penalty_factor = 0,
                           const unsigned int max_iter = 300,
                           const double rel_tol = 1e-6)
{
    Intsurv::CoxphReg object { Intsurv::CoxphReg(time, event, x) };
    object.regularized_fit(lambda, alpha, nlambda, lambda_min_ratio,
                           l1_penalty_factor, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef_mat,
        Rcpp::Named("en_coef") = object.en_coef_mat,
        Rcpp::Named("negLogL") = Intsurv::arma2rvec(object.negLogL_vec),
        Rcpp::Named("nObs") = object.nObs,
        Rcpp::Named("lambda_max") = object.l1_lambda_max,
        Rcpp::Named("lambda") = Intsurv::arma2rvec(object.lambda_vec),
        Rcpp::Named("alpha") = object.alpha,
        Rcpp::Named("l1_penalty_factor") = object.l1_penalty_factor
        );
}
