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
#include "logistic_reg.h"
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
    const double& em_rel_tol = 1e-3,
    const unsigned int& cox_mstep_max_iter = 30,
    const double& cox_mstep_rel_tol = 1e-3,
    const unsigned int& cure_mstep_max_iter = 30,
    const double& cure_mstep_rel_tol = 1e-6
    )
{
    // initialize
    Intsurv::CoxphReg cox_object { Intsurv::CoxphReg(time, event, cox_x) };
    arma::uvec cox_sort_ind { cox_object.get_sort_index() };

    arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
    arma::vec s_event { event.elem(cox_sort_ind) };
    arma::uvec case1_ind { arma::find(s_event > 0) };
    arma::uvec case2_ind { arma::find(s_event < 1) };

    Intsurv::LogisticReg cure_object {
        Intsurv::LogisticReg(cure_xx, s_event, cure_intercept)
    };
    if (cure_intercept) {
        // add intercept to the design matrix for cure rate model
        cure_xx = arma::join_horiz(arma::ones(cure_xx.n_rows), cure_xx);
    }
    arma::vec cox_beta { arma::zeros(cox_x.n_cols) };
    arma::vec cure_beta { arma::zeros(cure_xx.n_cols) };

    if (cox_start.n_elem == cox_x.n_cols) {
        cox_beta = cox_start;
    } else {
        arma::uvec tmp_idx { arma::find(event > 0) };
        Intsurv::CoxphReg tmp_object {
            Intsurv::CoxphReg(time.elem(tmp_idx),
                              event.elem(tmp_idx),
                              cox_x.rows(tmp_idx))
        };
        tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
        cox_beta = tmp_object.coef;
    }
    if (cure_start.n_elem == cure_xx.n_cols) {
        cure_beta = cure_start;
    } else {
        cure_object.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);
        cure_beta = cure_object.coef;
    }

    arma::vec p_vec { cure_object.predict(cure_beta) };
    cox_object.compute_haz_surv_time(cox_beta);

    size_t i {0};
    arma::vec estep_v { s_event }, log_v {0};
    double numer_j {0}, denom_j {0}, tol1 {0}, tol2 {0};
    double obs_ell {0};

    // main loop of EM algorithm
    while (true) {

        // E-step: compute v vector
        for (size_t j: case2_ind) {
            numer_j = p_vec(j) * cox_object.S_time(j);
            denom_j = numer_j + 1 - p_vec(j);
            estep_v(j) = numer_j / denom_j;
        }
        log_v = arma::log(estep_v);

        // M-step for the survival layer
        cox_object.set_offset(log_v);
        cox_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);

        // M-step for the Cure layer
        cure_object.update_y(estep_v);
        cure_object.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);

        // check convergence
        tol1 = Intsurv::rel_l2_norm(cox_object.coef, cox_beta);
        tol2 = Intsurv::rel_l2_norm(cure_object.coef, cure_beta);
        cox_beta = cox_object.coef;
        cure_beta = cure_object.coef;

        // update to last estimates
        cox_object.compute_haz_surv_time();
        p_vec = cure_object.predict(cure_beta);

        if ((tol1 < em_rel_tol && tol2 < em_rel_tol) ||
            i > em_max_iter) {
            // compute observed data likelihood
            // for case 1
            for (size_t j: case1_ind) {
                obs_ell += std::log(p_vec(j)) +
                    std::log(cox_object.h_time(j)) +
                    std::log(cox_object.S_time(j));
            }
            // for case 2
            for (size_t j: case2_ind) {
                obs_ell += std::log(
                    p_vec(j) * cox_object.S_time(j) + (1 - p_vec(j))
                    );
            }
            break;
        }
        // update iter
        ++i;
    }
    return Rcpp::List::create(
        Rcpp::Named("cox_beta") = cox_beta,
        Rcpp::Named("cure_beta") = cure_beta,
        Rcpp::Named("obs_ell") = obs_ell
        );
}


// fit Cox cure rate model by EM algorithm with adaptive lasso penalty
// [[Rcpp::export]]
Rcpp::List coxph_cure_reg(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_lambda = 0,
    const unsigned int& cox_nlambda = 1,
    double cox_lambda_min_ratio = 1e-4,
    const arma::vec& cox_penalty_factor = 0,
    const arma::vec& cure_lambda = 0,
    const unsigned int& cure_nlambda = 1,
    double cure_lambda_min_ratio = 1e-4,
    const arma::vec& cure_penalty_factor = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 200,
    const double& em_rel_tol = 1e-4,
    const unsigned int& cox_mstep_max_iter = 10,
    const double& cox_mstep_rel_tol = 1e-2,
    const unsigned int& cure_mstep_max_iter = 10,
    const double& cure_mstep_rel_tol = 1e-2
    )
{
    // initialize
    Intsurv::CoxphReg cox_object { Intsurv::CoxphReg(time, event, cox_x) };
    arma::uvec cox_sort_ind { cox_object.get_sort_index() };

    arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
    arma::vec s_event { event.elem(cox_sort_ind) };
    arma::uvec case1_ind { arma::find(s_event > 0) };
    arma::uvec case2_ind { arma::find(s_event < 1) };

    Intsurv::LogisticReg cure_object {
        Intsurv::LogisticReg(cure_xx, s_event, cure_intercept)
    };
    if (cure_intercept) {
        // add intercept to the design matrix for cure rate model
        cure_xx = arma::join_horiz(arma::ones(cure_xx.n_rows), cure_xx);
    }

    // penalty factor for Cox model
    arma::vec cox_penalty { arma::ones(cox_x.n_cols) };
    if (cox_penalty_factor.n_elem == cox_x.n_cols) {
        // re-scale so that sum(factor) = number of predictors
        cox_penalty = cox_penalty_factor * cox_x.n_cols /
            arma::sum(cox_penalty_factor);
    }

    // penalty factor for Cure model
    arma::vec cure_penalty { arma::ones(cure_x.n_cols) };
    if (cure_penalty_factor.n_elem == cure_x.n_cols) {
        // re-scale so that sum(factor) = number of predictors
        cure_penalty = cure_penalty_factor * cure_x.n_cols /
            arma::sum(cure_penalty_factor);
    }

    // all zeros coef
    arma::vec cox_beta { arma::zeros(cox_x.n_cols) };
    arma::vec cure_beta { arma::zeros(cure_xx.n_cols) };
    // the maximum (large enough) lambda that results in all-zero estimates
    arma::vec cox_grad_zero { arma::abs(cox_object.gradient(cox_beta)) };
    double cox_lambda_max {
        arma::max(cox_grad_zero / cox_penalty) / cox_x.n_rows
    };
    double log_cox_lambda_max { std::log(cox_lambda_max) };
    arma::vec cure_grad_zero { arma::abs(cure_object.gradient(cure_beta)) };
    double cure_lambda_max {
        arma::max(cure_grad_zero.tail(cure_penalty.n_elem) / cure_penalty) /
        cure_x.n_rows
    };
    double log_cure_lambda_max { std::log(cure_lambda_max) };

    // construct lambda sequence
    if (cox_x.n_cols > cox_x.n_rows) {
        cox_lambda_min_ratio = std::sqrt(cox_lambda_min_ratio);
    }
    if (cure_x.n_cols > cure_x.n_rows) {
        cure_lambda_min_ratio = std::sqrt(cure_lambda_min_ratio);
    }
    arma::vec cox_lambda_seq {
        arma::zeros(std::max(cox_nlambda, cox_lambda.n_elem))
    };
    arma::vec cure_lambda_seq {
        arma::zeros(std::max(cox_nlambda, cox_lambda.n_elem))
    };
    if (cox_lambda.n_elem == 1 && cox_nlambda > 1) {
        cox_lambda_seq = arma::exp(
            arma::linspace(log_cox_lambda_max,
                           log_cox_lambda_max +
                           std::log(cox_lambda_min_ratio),
                           cox_nlambda)
            );
    } else {
        cox_lambda_seq = cox_lambda;
    }
    if (cure_lambda.n_elem == 1 && cure_nlambda > 1) {
        cure_lambda_seq = arma::exp(
            arma::linspace(log_cure_lambda_max,
                           log_cure_lambda_max +
                           std::log(cure_lambda_min_ratio),
                           cure_nlambda)
            );
    } else {
        cure_lambda_seq = cure_lambda;
    }

    // update start estimates
    if (cox_start.n_elem == cox_x.n_cols) {
        cox_beta = cox_start;
    }
    if (cure_start.n_elem == cure_xx.n_cols) {
        cure_beta = cure_start;
    }

    arma::vec p_vec { cure_object.predict(cure_beta) };
    cox_object.compute_haz_surv_time(cox_beta);

    size_t i {0}, mat_idx {0};
    arma::vec estep_v { s_event }, log_v {0};
    double numer_j {0}, denom_j {0}, tol1 {0}, tol2 {0};
    double obs_ell {0};

    // set up the estimate matrix
    unsigned int est_ncol { cox_lambda_seq.n_elem * cure_lambda_seq.n_elem};
    arma::vec obs_ell_vec { arma::zeros(est_ncol) };
    arma::mat cox_beta_mat { arma::zeros(cox_x.n_cols, est_ncol) };
    arma::mat cure_beta_mat { arma::zeros(cure_xx.n_cols, est_ncol) };

    for (size_t l1 {0}; l1 < cox_lambda_seq.n_elem; ++l1) {
        for (size_t l2 {0}; l2 < cure_lambda_seq.n_elem; ++l2) {
            i = 0;
            obs_ell = 0;
            mat_idx = l1 * cure_lambda_seq.n_elem + l2;

            // main loop of EM algorithm
            while (true) {

                // E-step: compute v vector
                for (size_t j: case2_ind) {
                    numer_j = p_vec(j) * cox_object.S_time(j);
                    denom_j = numer_j + 1 - p_vec(j);
                    estep_v(j) = numer_j / denom_j;
                }
                log_v = arma::log(estep_v);

                // M-step for the survival layer
                cox_object.set_offset(log_v);
                cox_object.regularized_fit(
                    cox_lambda_seq(l1), cox_penalty,
                    cox_beta, cox_mstep_max_iter,
                    cox_mstep_rel_tol
                    );

                // M-step for the Cure layer
                cure_object.update_y(estep_v);
                cure_object.regularized_fit(
                    cure_lambda_seq(l2), cure_penalty,
                    cure_beta, cure_mstep_max_iter,
                    cure_mstep_rel_tol
                    );

                // check convergence
                tol1 = Intsurv::rel_l2_norm(cox_object.coef, cox_beta);
                tol2 = Intsurv::rel_l2_norm(cure_object.coef, cure_beta);
                cox_beta = cox_object.coef;
                cure_beta = cure_object.coef;

                // update to last estimates
                cox_object.compute_haz_surv_time();
                p_vec = cure_object.predict(cure_beta);

                if ((tol1 < em_rel_tol && tol2 < em_rel_tol) ||
                    i > em_max_iter) {
                    // compute observed data likelihood
                    // for case 1
                    for (size_t j: case1_ind) {
                        obs_ell += std::log(p_vec(j)) +
                            std::log(cox_object.h_time(j)) +
                            std::log(cox_object.S_time(j));
                    }
                    // for case 2
                    for (size_t j: case2_ind) {
                        obs_ell += std::log(p_vec(j) * cox_object.S_time(j) +
                                            (1 - p_vec(j)));
                    }
                    break;
                }
                // update iter
                ++i;
            }
            // break to here
            cox_beta_mat.col(mat_idx) = cox_beta;
            cure_beta_mat.col(mat_idx) = cure_beta;
            obs_ell_vec.row(mat_idx) = obs_ell;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("cox_beta_mat") = cox_beta_mat,
        Rcpp::Named("cure_beta_mat") = cure_beta_mat,
        Rcpp::Named("obs_ell") = obs_ell_vec
        );
}
