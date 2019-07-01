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

#include "coxph_cure_uncer.h"

// Cox cure model with uncertain events without regularization
// [[Rcpp::export]]
Rcpp::List coxph_cure_uncer(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 300,
    const double& em_rel_tol = 1e-5,
    const unsigned int& cox_mstep_max_iter = 100,
    const double& cox_mstep_rel_tol = 1e-5,
    const unsigned int& cure_mstep_max_iter = 100,
    const double& cure_mstep_rel_tol = 1e-5,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool& spline_start = false,
    const unsigned int& iSpline_num_knots = 3,
    const unsigned int& iSpline_degree = 2,
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    Intsurv::CoxphCureUncer obj {
        Intsurv::CoxphCureUncer(time, event, cox_x, cure_x, cure_intercept,
                                cox_standardize, cure_standardize)
    };
    obj.fit(cox_start, cure_start,
            em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol,
            spline_start, iSpline_num_knots, iSpline_degree,
            tail_completion, tail_tau,
            pmin, early_stop, verbose
        );
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est),
            Rcpp::Named("hc") = Intsurv::arma2rvec(obj.hc_est),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(obj.Hc_est),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(obj.Sc_est)
            ),
        Rcpp::Named("prediction") = Rcpp::List::create(
            Rcpp::Named("cox_xBeta") = Intsurv::arma2rvec(obj.cox_xBeta),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xBeta),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob),
            Rcpp::Named("estep_cured") = Intsurv::arma2rvec(obj.estep_cured),
            Rcpp::Named("estep_event") = Intsurv::arma2rvec(obj.estep_event),
            Rcpp::Named("estep_censor") = Intsurv::arma2rvec(obj.estep_censor)
            ),
        Rcpp::Named("goodness") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.nObs,
            Rcpp::Named("coef_df") = obj.coef_df,
            Rcpp::Named("negLogL") = obj.negLogL,
            Rcpp::Named("c_index") = obj.c_index,
            Rcpp::Named("bic1") = obj.bic1,
            Rcpp::Named("bic2") = obj.bic2
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.num_iter
            )
        );
}


// fit regularized Cox cure rate model with uncertain events
// by EM algorithm, where the M-step utilized CMD algoritm
// [[Rcpp::export]]
Rcpp::List coxph_cure_uncer_reg(
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
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool& spline_start = false,
    const unsigned int& iSpline_num_knots = 3,
    const unsigned int& iSpline_degree = 2,
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    Intsurv::CoxphCureUncer obj {
        Intsurv::CoxphCureUncer(time, event, cox_x, cure_x, cure_intercept,
                                cox_standardize, cure_standardize)
    };
    obj.regularized_fit(
        cox_l1_lambda, cox_l2_lambda,
        cure_l1_lambda, cure_l2_lambda,
        cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start, em_max_iter, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol,
        spline_start, iSpline_num_knots, iSpline_degree,
        tail_completion, tail_tau,
        pmin, early_stop, verbose
        );
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("cox_en_coef") = Intsurv::arma2rvec(obj.cox_en_coef),
        Rcpp::Named("cure_en_coef") = Intsurv::arma2rvec(obj.cure_en_coef),

        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est),
            Rcpp::Named("hc") = Intsurv::arma2rvec(obj.hc_est),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(obj.Hc_est),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(obj.Sc_est)
            ),
        Rcpp::Named("prediction") = Rcpp::List::create(
            Rcpp::Named("cox_xBeta") = Intsurv::arma2rvec(obj.cox_xBeta),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xBeta),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob),
            Rcpp::Named("estep_cured") = Intsurv::arma2rvec(obj.estep_cured),
            Rcpp::Named("estep_event") = Intsurv::arma2rvec(obj.estep_event),
            Rcpp::Named("estep_censor") = Intsurv::arma2rvec(obj.estep_censor)
            ),
        Rcpp::Named("goodness") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.nObs,
            Rcpp::Named("coef_df") = obj.coef_df,
            Rcpp::Named("negLogL") = obj.negLogL,
            Rcpp::Named("c_index") = obj.c_index,
            Rcpp::Named("bic1") = obj.bic1,
            Rcpp::Named("bic2") = obj.bic2
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
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
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.num_iter
            )
        );
}