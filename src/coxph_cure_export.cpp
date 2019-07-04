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

#include "assessment.h"
#include "coxph_cure.h"


// fit regular Cox cure rate model by EM algorithm
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& bootstrap = 0,
    const unsigned int& em_max_iter = 1000,
    const double& em_rel_tol = 1e-4,
    const unsigned int& cox_mstep_max_iter = 200,
    const double& cox_mstep_rel_tol = 1e-4,
    const unsigned int& cure_mstep_max_iter = 200,
    const double& cure_mstep_rel_tol = 1e-6,
    const bool& firth = false,
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    // define object
    Intsurv::CoxphCure obj {
        Intsurv::CoxphCure(time, event, cox_x, cure_x, cure_intercept)
    };
    // model-fitting
    obj.fit(cox_start, cure_start,
            em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol,
            firth, tail_completion, tail_tau,
            pmin, early_stop, verbose);
    // initialize bootstrap estimates
    arma::mat boot_cox_coef_mat, boot_cure_coef_mat;
    if (bootstrap > 0) {
        boot_cox_coef_mat = arma::zeros(obj.cox_coef.n_elem, bootstrap);
        boot_cure_coef_mat = arma::zeros(obj.cure_coef.n_elem, bootstrap);
        arma::uvec case1_ind { arma::find(event > 0) };
        arma::uvec case2_ind { arma::find(event < 1) };
        for (size_t i {0}; i < bootstrap; ++i) {
            // generate a bootstrap sample
            arma::uvec boot_ind {
                Intsurv::vec_union(
                    Intsurv::bootstrap_sample(case1_ind),
                    Intsurv::bootstrap_sample(case2_ind)
                    )
            };
            Intsurv::CoxphCure boot_obj {
                Intsurv::CoxphCure(time.elem(boot_ind),
                                   event.elem(boot_ind),
                                   cox_x.rows(boot_ind),
                                   cure_x.rows(boot_ind),
                                   cure_intercept)
            };
            // fit the bootstarp sample
            boot_obj.fit(cox_start, cure_start,
                         em_max_iter, em_rel_tol,
                         cox_mstep_max_iter, cox_mstep_rel_tol,
                         cure_mstep_max_iter, cure_mstep_rel_tol,
                         firth, tail_completion, tail_tau,
                         pmin, early_stop, 0);
            boot_cox_coef_mat.col(i) = boot_obj.cox_coef;
            boot_cure_coef_mat.col(i) = boot_obj.cure_coef;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est)
            ),
        Rcpp::Named("prediction") = Rcpp::List::create(
            Rcpp::Named("cox_xBeta") = Intsurv::arma2rvec(obj.cox_xBeta),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xBeta),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob),
            Rcpp::Named("estep_cured") =
            Intsurv::arma2rvec(obj.estep_cured),
            Rcpp::Named("estep_susceptible") =
            Intsurv::arma2rvec(obj.estep_susceptible)
            ),
        Rcpp::Named("goodness") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.nObs,
            Rcpp::Named("coef_df") = obj.coef_df,
            Rcpp::Named("negLogL") = obj.negLogL,
            Rcpp::Named("c_index") = obj.c_index,
            Rcpp::Named("bic1") = obj.bic1,
            Rcpp::Named("bic2") = obj.bic2
            ),
        Rcpp::Named("bootstrap") = Rcpp::List::create(
            Rcpp::Named("B") = bootstrap,
            Rcpp::Named("cox_coef_mat") = boot_cox_coef_mat.t(),
            Rcpp::Named("cure_coef_mat") = boot_cure_coef_mat.t()
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.num_iter
            )
        );
}


// fit regularized Cox cure rate model by EM algorithm,
// where the M-step utilized CMD algoritm
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure_reg(
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
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    Intsurv::CoxphCure obj {
        Intsurv::CoxphCure(time, event, cox_x, cure_x, cure_intercept,
                           cox_standardize, cure_standardize)
    };
    obj.regularized_fit(
        cox_l1_lambda, cox_l2_lambda,
        cure_l1_lambda, cure_l2_lambda,
        cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start,
        em_max_iter, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol,
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
            Rcpp::Named("h0_est") = Intsurv::arma2rvec(obj.h0_est),
            Rcpp::Named("H0_est") = Intsurv::arma2rvec(obj.H0_est),
            Rcpp::Named("S0_est") = Intsurv::arma2rvec(obj.S0_est)
            ),
        Rcpp::Named("prediction") = Rcpp::List::create(
            Rcpp::Named("cox_xBeta") = Intsurv::arma2rvec(obj.cox_xBeta),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xBeta),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob),
            Rcpp::Named("estep_cured") =
            Intsurv::arma2rvec(obj.estep_cured),
            Rcpp::Named("estep_susceptible") =
            Intsurv::arma2rvec(obj.estep_susceptible)
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
