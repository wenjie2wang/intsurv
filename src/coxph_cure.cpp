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
Rcpp::NumericVector coxph_cure(const arma::vec& time,
                               const arma::vec& event,
                               const arma::mat& cox_x,
                               const arma::mat& cure_x,
                               const bool cure_intercept = true,
                               const arma::vec& cox_start = 0,
                               const arma::vec& cure_start = 0,
                               const unsigned int& max_iter = 1000,
                               const double& rel_tol = 1e-6)
{
    // initialize
    Intsurv::CoxphReg cox_object { Intsurv::CoxphReg(time, event, cox_x) };
    arma::uvec cox_sort_ind { cox_object.get_sort_index() };

    arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
    arma::vec s_event { event.elem(cox_sort_ind) };
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
        tmp_object.fit(cox_beta, max_iter, rel_tol);
        cox_beta = tmp_object.coef;
    }
    if (cure_start.n_elem == cure_xx.n_cols) {
        cure_beta = cure_start;
    } else {
        cure_object.fit(cure_beta, max_iter, rel_tol);
        cure_beta = cure_object.coef;
    }

    arma::vec p_vec { cure_object.predict(cure_beta) };
    cox_object.set_offset(arma::log(s_event));
    cox_object.compute_haz_surv_time(cox_beta);

    arma::uvec miss_ind { arma::find(s_event < 1) };
    size_t i {0};
    arma::vec estep_v { s_event }, log_v {0};
    double numer_j {0}, denom_j {0}, tol1 {0}, tol2 {0};

    // main loop of EM algorithm
    while (i < max_iter) {

        // E-step: compute v vector
        for (size_t j: miss_ind) {
            numer_j = p_vec(j) * cox_object.S_time(j);
            denom_j = numer_j + 1 - p_vec(j);
            estep_v(j) = numer_j / denom_j;
        }
        log_v = arma::log(estep_v);

        // M-step for the survival layer
        cox_object.set_offset(log_v);
        cox_object.fit(cox_beta, max_iter, rel_tol);

        // M-step for the Cure layer
        cure_object.update_y(estep_v);
        cure_object.fit(cure_beta, max_iter, rel_tol);

        // check convergence
        tol1 = Intsurv::rel_l2_norm(cox_object.coef, cox_beta);
        tol2 = Intsurv::rel_l2_norm(cure_object.coef, cure_beta);
        cox_beta = cox_object.coef;
        cure_beta = cure_object.coef;
        if (tol1 < rel_tol && tol2 < rel_tol) {
            break;
        }

        // update to last estimates
        cox_object.compute_haz_surv_time();
        p_vec = cure_object.predict(cure_beta);

        // update iter
        ++i;
    }
    arma::vec out { arma::join_vert(cox_beta, cure_beta) };
    return Rcpp::NumericVector(out.begin(), out.end());
}
