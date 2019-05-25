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

#include "coxph_cure.h"
#include "utils.h"


// fit regular Cox cure rate model by EM algorithm
// [[Rcpp::export]]
Rcpp::List coxph_cure(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 1000,
    const double& em_rel_tol = 1e-4,
    const unsigned int& cox_mstep_max_iter = 200,
    const double& cox_mstep_rel_tol = 1e-4,
    const unsigned int& cure_mstep_max_iter = 200,
    const double& cure_mstep_rel_tol = 1e-6,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool verbose = false
    )
{
    Intsurv::CoxphCure obj {
        Intsurv::CoxphCure(time, event, cox_x, cure_x, cure_intercept,
                           cox_standardize, cure_standardize)
    };
    obj.fit(cox_start, cure_start, em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol, verbose);
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("unique_time") = Intsurv::arma2rvec(obj.unique_time),
        Rcpp::Named("h0_est") = Intsurv::arma2rvec(obj.h0_est),
        Rcpp::Named("H0_est") = Intsurv::arma2rvec(obj.H0_est),
        Rcpp::Named("S0_est") = Intsurv::arma2rvec(obj.S0_est),
        Rcpp::Named("negLogL") = obj.negLogL,
        Rcpp::Named("nObs") = obj.nObs,
        Rcpp::Named("coef_df") = obj.coef_df,
        Rcpp::Named("num_iter") = obj.num_iter
        );
}


// fit regularized Cox cure rate model by EM algorithm,
// where the M-step utilized CMD algoritm
// [[Rcpp::export]]
Rcpp::List coxph_cure_reg(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const double& cox_l1_lambda = 0,
    const double& cox_l2_lambda = 0,
    const arma::vec& cox_l1_penalty_factor = 0,
    const double& cure_l1_lambda = 0,
    const double& cure_l2_lambda = 0,
    const arma::vec& cure_l1_penalty_factor = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 500,
    const double& em_rel_tol = 1e-5,
    const unsigned int& cox_mstep_max_iter = 200,
    const double& cox_mstep_rel_tol = 1e-4,
    const unsigned int& cure_mstep_max_iter = 200,
    const double& cure_mstep_rel_tol = 1e-4,
    const bool verbose = false
    )
{
    Intsurv::CoxphCure obj {
        Intsurv::CoxphCure(time, event, cox_x, cure_x, cure_intercept)
    };
    obj.regularized_fit(
        cox_l1_lambda, cox_l2_lambda,
        cure_l1_lambda, cure_l2_lambda,
        cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start, em_max_iter, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol, verbose
        );
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("cox_en_coef") = Intsurv::arma2rvec(obj.cox_en_coef),
        Rcpp::Named("cure_en_coef") = Intsurv::arma2rvec(obj.cure_en_coef),

        Rcpp::Named("unique_time") = Intsurv::arma2rvec(obj.unique_time),
        Rcpp::Named("h0_est") = Intsurv::arma2rvec(obj.h0_est),
        Rcpp::Named("H0_est") = Intsurv::arma2rvec(obj.H0_est),
        Rcpp::Named("S0_est") = Intsurv::arma2rvec(obj.S0_est),

        Rcpp::Named("negLogL") = obj.negLogL,
        Rcpp::Named("nObs") = obj.nObs,
        Rcpp::Named("coef_df") = obj.coef_df,
        Rcpp::Named("num_iter") = obj.num_iter,

        Rcpp::Named("cox_l1_lambda_max") = obj.cox_l1_lambda_max,
        Rcpp::Named("cox_l1_lambda") = obj.cox_l1_lambda,
        Rcpp::Named("cox_l2_lambda") = obj.cox_l2_lambda,
        Rcpp::Named("cox_l1_penalty_factor") =
        Intsurv::arma2rvec(obj.cox_l1_penalty_factor),

        Rcpp::Named("cure_l1_lambda_max") = obj.cure_l1_lambda_max,
        Rcpp::Named("cure_l1_lambda") = obj.cure_l1_lambda,
        Rcpp::Named("cure_l2_lambda") = obj.cure_l2_lambda,
        Rcpp::Named("cure_l1_penalty_factor") =
        Intsurv::arma2rvec(obj.cure_l1_penalty_factor)
        );
}
