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

#include "coxph_reg.h"
#include "utils.h"

// fitting Cox model by monotonic quadratic approximation algorithm
// that allows non-integer "event" and tied events
// reference: Bohning and Lindsay (1988) SIAM
// [[Rcpp::export]]
arma::vec rcpp_coxph(const arma::vec& time,
                     const arma::vec& event,
                     const arma::mat& x,
                     const arma::vec& start = 0,
                     const unsigned int max_iter = 1000,
                     const double rel_tol = 1e-6)
{
    Intsurv::RcppCoxph object { Intsurv::RcppCoxph(time, event, x) };
    arma::mat b_mat { object.bl_cov_lowerbound_mat() };
    arma::mat inv_b_mat { arma::inv_sympd(b_mat) };
    arma::vec beta0 { arma::zeros(x.n_cols) };
    if (start.n_elem == x.n_cols) {
        beta0 = start;
    }
    arma::vec beta { beta0 }, h_vec { beta0 }, grad_vec { beta0 };
    size_t i {0};
    double b_new {0}, alpha {0};
    while (i < max_iter) {
        grad_vec = object.gradient(beta0);
        h_vec = - inv_b_mat * grad_vec;
        b_new = object.bl_step_lowerbound(x, h_vec);
        alpha = - arma::as_scalar(Intsurv::crossprod(h_vec, grad_vec)) / b_new;
        beta = beta0 + alpha * h_vec;
        if (Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
            break;
        }
        // update beta
        beta0 = beta;
        i++;
    }
    return beta;
}


// run one cycle of coordinate descent over a given active set
void reg_coxph_update(arma::vec& beta,
                      arma::uvec& is_active,
                      const Intsurv::RcppCoxph& object,
                      const arma::vec& penalty,
                      const bool update_active = false)
{
    // compute D_k
    arma::vec d_vec { object.cmd_lowerbound() };
    double dlj { 0 };
    double n_sample { static_cast<double>(object.sample_size()) };
    for (size_t j {0}; j < beta.n_elem; ++j) {
        if (is_active[j]) {
            dlj = object.gradient(beta, j) / n_sample;
            // update beta
            beta[j] = Intsurv::soft_threshold(
                d_vec[j] * beta[j] - dlj, penalty[j]) / d_vec[j];
            if (update_active) {
                // check if it has been shrinkaged to zero
                if (Intsurv::isAlmostEqual(beta[j], 0)) {
                    is_active[j] = 0;
                } else {
                    is_active[j] = 1;
                }
            }
        }
    }
}


// fitting regularized Cox model with coordinate-majorizatio-descent algorithm
// that allows non-integer "event" and tied events
// [[Rcpp::export]]
arma::vec rcpp_reg_coxph(const arma::vec& time,
                         const arma::vec& event,
                         const arma::mat& x,
                         const double lambda = 0,
                         arma::vec penalty_factor = 0,
                         const arma::vec& start = 0,
                         const unsigned int max_iter = 1000,
                         const double rel_tol = 1e-6)
{
    Intsurv::RcppCoxph object { Intsurv::RcppCoxph(time, event, x) };
    // declarations
    arma::vec beta0 { arma::zeros(x.n_cols) };
    if (start.n_elem == x.n_cols) {
        beta0 = start;
    }
    arma::vec beta { beta0 };
    arma::uvec is_active { arma::ones<arma::uvec>(x.n_cols) };
    arma::uvec is_active_stored { is_active };

    if (penalty_factor.n_elem == x.n_cols) {
        // re-scale so that sum(factor) = number of predictors
        penalty_factor = penalty_factor * x.n_cols / arma::sum(penalty_factor);
    } else {
        penalty_factor = arma::ones(x.n_cols);
    }
    penalty_factor *= lambda;

    // the lower bound for second derivative in cmd
    arma::vec d_vec { object.cmd_lowerbound() };

    size_t i {0};
    // use active-set if p > n ("helps when p >> n")
    if (x.n_cols > x.n_rows) {
        size_t ii {0};
        while (i < max_iter) {
            // cycles over the active set
            while (ii < max_iter) {
                reg_coxph_update(beta, is_active, object,
                                 penalty_factor, true);
                if (Intsurv::isAlmostEqual(Intsurv::l2_norm(beta), 0) ||
                    Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
                    break;
                }
                beta0 = beta;
                ii++;
            }
            is_active_stored = is_active;
            // run a full cycle over the converged beta
            is_active = arma::ones<arma::uvec>(x.n_cols);
            reg_coxph_update(beta, is_active, object,
                             penalty_factor, true);
            // check two active sets coincide
            if (arma::sum(arma::abs(is_active - is_active_stored)) > 0) {
                // if different, repeat this process
                ii = 0;
                i++;
            } else {
                break;
            }
        }
    } else {
        // regular coordinate descent
        while (i < max_iter) {
            reg_coxph_update(beta, is_active, object,
                             penalty_factor, false);
            if (Intsurv::isAlmostEqual(Intsurv::l2_norm(beta), 0) ||
                Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
                break;
            }
            beta0 = beta;
            i++;
        }
    }
    return beta;
}


// [[Rcpp::export]]
Rcpp::NumericVector grad_fun(const arma::vec& beta,
                             const arma::vec& time,
                             const arma::vec& event,
                             const arma::mat& x)
{
    Intsurv::RcppCoxph object { Intsurv::RcppCoxph(time, event, x) };
    arma::vec grad { object.gradient(beta) };
    return Rcpp::NumericVector(grad.begin(), grad.end());
}

// [[Rcpp::export]]
double obj_fun(const arma::vec& beta,
               const arma::vec& time,
               const arma::vec& event,
               const arma::mat& x)
{
    Intsurv::RcppCoxph object { Intsurv::RcppCoxph(time, event, x) };
    return object.objective(beta);
}

// [[Rcpp::export]]
double grad_fun_k(const arma::vec& beta,
                  const unsigned int k,
                  const arma::vec& time,
                  const arma::vec& event,
                  const arma::mat& x)
{
    Intsurv::RcppCoxph object { Intsurv::RcppCoxph(time, event, x) };
    return object.gradient(beta, k);
}
