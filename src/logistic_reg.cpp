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

// fitting regular logistic model by monotonic quadratic approximation algorithm
// non-integer y vector is allowed
// reference: Bohning and Lindsay (1988) SIAM
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_logistic(const arma::mat& x,
                                  const arma::vec& y,
                                  const arma::vec start = 0,
                                  const unsigned int max_iter = 1000,
                                  const double rel_tol = 1e-6)
{
    Intsurv::LogisticReg object { Intsurv::LogisticReg(x, y) };
    // declarations
    arma::vec beta0 { arma::zeros(x.n_cols) };
    if (start.n_elem == x.n_cols) {
        beta0 = start;
    }
    arma::vec beta { beta0 };
    arma::vec eta { arma::zeros(y.n_elem) };
    arma::vec y_hat { eta };
    arma::mat iter_mat { 4 * arma::inv_sympd(x.t() * x) * x.t() };
    size_t i {0};
    while (i < max_iter) {
        eta = x * beta0;
        y_hat = object.linkinv(eta);
        beta = beta0 + iter_mat * (y - y_hat);
        if (Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
            break;
        }
        // update beta
        beta0 = beta;
        i++;
    }
    if (i == max_iter) {
        Rcpp::Rcout << "Reached the maximum iteration number." << std::endl;
    }
    return Rcpp::NumericVector(beta.begin(), beta.end());
}


// run one cycle of coordinate descent over a given active set
void reg_logistic_update(arma::vec& beta,
                         arma::uvec& is_active,
                         const Intsurv::LogisticReg& object,
                         const arma::rowvec& b_vec,
                         const arma::vec& penalty,
                         const bool update_active = false)
{
    double dlj { 0 };
    unsigned int n { object.sample_size() };
    double n_y { static_cast<double>(n) };
    for (size_t j {0}; j < beta.n_elem; ++j) {
        if (is_active[j]) {
            dlj = object.gradient(beta, j) / n_y;
            // update beta
            beta[j] = Intsurv::soft_threshold(
                b_vec[j] * beta[j] - dlj, penalty[j]) / b_vec[j];
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



// regularized logistic model by coordinate-majorization-descent algorithm
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_reg_logistic(const arma::mat& x,
                                      const arma::vec& y,
                                      const double lambda = 0,
                                      arma::vec penalty_factor = 0,
                                      const arma::vec start = 0,
                                      const unsigned int max_iter = 1000,
                                      const double rel_tol = 1e-6)
{
    Intsurv::LogisticReg object { Intsurv::LogisticReg(x, y) };
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
    arma::rowvec b_vec { arma::sum(arma::square(x), 0) / (4 * y.n_elem) };
    size_t i {0};

    // use active-set if p > n ("helps when p >> n")
    if (x.n_cols > x.n_rows) {
        size_t ii {0};
        while (i < max_iter) {
            // cycles over the active set
            while (ii < max_iter) {
                reg_logistic_update(beta, is_active, object,
                                    b_vec, penalty_factor, true);
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
            reg_logistic_update(beta, is_active, object,
                                b_vec, penalty_factor, true);
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
            reg_logistic_update(beta, is_active, object,
                                b_vec, penalty_factor, false);
            if (Intsurv::isAlmostEqual(Intsurv::l2_norm(beta), 0) ||
                Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
                break;
            }
            beta0 = beta;
            i++;
        }
    }
    // throw message/warning if the maximum iteration number is reached
    if (i == max_iter) {
        Rcpp::Rcout << "Reached the maximum iteration number." << std::endl;
    }
    return Rcpp::NumericVector(beta.begin(), beta.end());
}
