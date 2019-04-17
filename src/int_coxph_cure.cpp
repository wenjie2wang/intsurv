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
#include "nonparametric.h"
#include "utils.h"

// #include "../working/debug_arma.h"

// integrative Cox cure rate model by EM algorithm
// [[Rcpp::export]]
Rcpp::List int_coxph_cure(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const double& prob_event_start = 0.5,
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

    arma::vec s_event { cox_object.get_event() };
    arma::vec s_time { cox_object.get_time() };
    arma::mat s_cox_x { cox_object.get_x() };

    arma::uvec case1_ind { arma::find(s_event == 1) };
    arma::uvec case2_ind { arma::find(s_event == 0) };
    arma::uvec case12_ind { Intsurv::vec_union(case1_ind, case2_ind) };
    arma::uvec case13_ind { arma::find(s_event > 0) };
    arma::uvec case23_ind { arma::find(s_event < 1) };
    arma::uvec case3_ind { Intsurv::vec_intersection(case13_ind, case23_ind) };

    arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
    Intsurv::LogisticReg cure_object {
        Intsurv::LogisticReg(cure_xx, s_event, cure_intercept)
    };
    if (cure_intercept) {
        // add intercept to the design matrix for cure rate model
        cure_xx = arma::join_horiz(arma::ones(cure_xx.n_rows), cure_xx);
    }

    arma::vec cox_beta { arma::zeros(cox_x.n_cols) };
    if (cox_start.n_elem == cox_x.n_cols) {
        cox_beta = cox_start;
    } else {
        Intsurv::CoxphReg tmp_object {
            Intsurv::CoxphReg(s_time.elem(case1_ind),
                              s_event.elem(case1_ind),
                              s_cox_x.rows(case1_ind))
        };
        tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
        cox_beta = tmp_object.coef;
    }
    arma::vec cox_exp_x_beta = arma::exp(s_cox_x * cox_beta);

    arma::vec cure_beta { arma::zeros(cure_xx.n_cols) };
    if (cure_start.n_elem == cure_xx.n_cols) {
        cure_beta = cure_start;
    } else {
        cure_object.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);
        cure_beta = cure_object.coef;
    }
    arma::vec p_vec { cure_object.predict(cure_beta) };

    // initialize baseline hazard functions
    // for events
    Intsurv::NelsonAalen nelen_event {
        Intsurv::NelsonAalen(s_time.elem(case12_ind),
                             s_event.elem(case12_ind))
    };
    cox_object.h0_time = nelen_event.step_inst_rate(s_time);
    cox_object.h_time = cox_object.h0_time % cox_exp_x_beta;
    cox_object.H0_time = nelen_event.step_cum_rate(s_time);
    cox_object.H_time = cox_object.H0_time % cox_exp_x_beta;
    cox_object.S0_time = arma::exp(- cox_object.H0_time);
    cox_object.S_time = arma::exp(- cox_object.H_time);
    // for censoring
    Intsurv::NelsonAalen nelen_censor {
        Intsurv::NelsonAalen(s_time.elem(case12_ind),
                             1 - s_event.elem(case12_ind))
    };
    cox_object.hc_time = nelen_censor.step_inst_rate(s_time);
    cox_object.Hc_time = nelen_censor.step_cum_rate(s_time);
    cox_object.Sc_time = arma::exp(- cox_object.Hc_time);

    size_t i {0};
    double numer_j {0};
    double tol1 {0}, tol2 {0};
    double prob_event { prob_event_start }, prob_event_new { prob_event };
    double m1 {0}, m2 {0}, m3 {0};
    double w1 {0}, w2 {0}, w1_sum {0}, w2_sum {0};
    double m12_common {0};
    arma::vec estep_m { s_event }, log_m { 0 };
    double obs_ell {0};

    // main loop of EM algorithm
    while (true) {

        obs_ell = 0;

        // E-step: compute the v vector for case 2
        for (size_t j: case2_ind) {
            numer_j = p_vec(j) * cox_object.S_time(j);
            estep_m(j) = 1 / ((1 - p_vec(j)) / numer_j + 1);
        }

        // reset for updating prob_event
        w1_sum = 0;
        w2_sum = 0;

        // Rcpp::Rcout << "\nh0(t)\n"
        //             << cox_object.h0_time << std::endl;
        // Rcpp::Rcout << "\nS0(t)\n"
        //             << cox_object.S0_time << std::endl;
        // Rcpp::Rcout << "\nh_c(t)\n"
        //             << cox_object.hc_time << std::endl;
        // Rcpp::Rcout << "\nS_c(t)\n"
        //             << cox_object.Sc_time << std::endl;

        // Rcpp::Rcout << "\nS_c(t)\n"
        //             << cox_object.Sc_time.elem(case3_ind) << std::endl;
        // Rcpp::Rcout << "\nh_c(t)\n"
        //             << cox_object.hc_time.elem(case3_ind) << std::endl;
        // Rcpp::Rcout << "\np_j\n" << p_vec.elem(case3_ind) << std::endl;

        arma::mat m123_mat { arma::zeros(case3_ind.n_elem, 3) };
        arma::mat w123_mat { arma::zeros(case3_ind.n_elem, 3) };
        size_t ii { 0 };

        // E-step: compute the w vector for case 3
        for (size_t j: case3_ind) {
            // the Sc_time(j) is cancelled out
            m12_common = p_vec(j) * cox_object.S_time(j);
            m1 = prob_event * cox_object.h_time(j) * m12_common;
            m2 = (1 - prob_event) * cox_object.hc_time(j) * m12_common;
            m3 = (1 - p_vec(j)) * cox_object.hc_time(j);

            w1 = 1 / ((m2 + m3) / m1 + 1);
            w1_sum += w1;
            w2 = 1 / ((m1 + m3) / m2 + 1);
            w2_sum += w2;
            estep_m(j) = w1 + w2;
            s_event(j) = w1;

            m123_mat(ii, 0) = cox_object.h_time(j) * cox_object.S_time(j);
            m123_mat(ii, 1) = cox_object.hc_time(j) * cox_object.Sc_time(j);
            m123_mat(ii, 2) = m3;

            w123_mat(ii, 0) = w1;
            w123_mat(ii, 1) = w2;
            w123_mat(ii, 2) = 1 / ((m1 + m2) / m3 + 1);;

            ++ii;
        }

        // Rcpp::Rcout << "\nm123_mat\n" << m123_mat << std::endl;
        // Rcpp::Rcout << "\ncol_sum_m123\n" << arma::sum(m123_mat) << std::endl;

        // Rcpp::Rcout << "\nw123_mat\n" << w123_mat << std::endl;
        // Rcpp::Rcout << "\ncol_sum_w123\n" << arma::sum(w123_mat) << std::endl;

        // Rcpp::Rcout << "\nM_j\n" << estep_m.elem(case3_ind) << std::endl;
        // Rcpp::Rcout << "\nevent\n" << s_event.elem(case3_ind) << std::endl;

        // M-step for the survival layer
        log_m = arma::log(estep_m);
        cox_object.set_offset(log_m);
        cox_object.update_event_weight(s_event);
        cox_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);

        // M-step for the cure rate layer
        cure_object.update_y(estep_m);
        cure_object.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);

        // Rcpp::Rcout << w1_sum + w2_sum << std::endl;

        // M-step for prob_event
        prob_event_new = 1 / (w2_sum / w1_sum + 1);

        Rcpp::Rcout << i << " times iteration: prob_event = " <<
            prob_event_new << std::endl;

        // check convergence
        tol1 = Intsurv::rel_l2_norm(cox_object.coef, cox_beta);
        tol2 = Intsurv::rel_l2_norm(cure_object.coef, cure_beta);
        // tol3 = std::abs(prob_event - prob_event_new);

        // update to last estimates
        cox_beta = cox_object.coef;
        cure_beta = cure_object.coef;
        prob_event = prob_event_new;
        cox_object.compute_haz_surv_time();
        cox_object.compute_censor_haz_surv_time();
        p_vec = cure_object.predict(cure_beta);

        // compute the observed data log-likelihood
        // for case 1
        for (size_t j: case1_ind) {
            obs_ell += std::log(p_vec(j)) +
                std::log(cox_object.h_time(j)) +
                std::log(cox_object.S_time(j)) +
                std::log(cox_object.Sc_time(j));
        }
        // for case 2
        for (size_t j: case2_ind) {
            obs_ell += std::log(p_vec(j) * cox_object.S_time(j) +
                                1 - p_vec(j)) +
                std::log(cox_object.Sc_time(j)) +
                std::log(cox_object.hc_time(j));
        }
        // for case 3
        for (size_t j: case3_ind) {
            m1 = p_vec(j) * prob_event * cox_object.h_time(j) *
                cox_object.S_time(j) * cox_object.Sc_time(j);
            m2 = p_vec(j) * (1 - prob_event) * cox_object.S_time(j) *
                cox_object.hc_time(j) * cox_object.Sc_time(j);
            m3 = (1 - p_vec(j)) * cox_object.hc_time(j) *
                cox_object.Sc_time(j);
            obs_ell += std::log(m1 + m2 + m3);
        }

        Rcpp::Rcout << "obs_ell:" << obs_ell << std::endl;

        if ((tol1 < em_rel_tol && tol2 < em_rel_tol) ||
            i > em_max_iter) {
            // break here
            break;
        }
        // update iter
        ++i;
    }
    return Rcpp::List::create(
        Rcpp::Named("cox_beta") = cox_beta,
        Rcpp::Named("cure_beta") = cure_beta,
        Rcpp::Named("prob_event") = prob_event,
        Rcpp::Named("obs_ell") = obs_ell
        );
}
