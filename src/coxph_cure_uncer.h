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

#ifndef COXPH_CURE_UNCER_H
#define COXPH_CURE_UNCER_H

#include <RcppArmadillo.h>
#include "coxph_reg.h"
#include "logistic_reg.h"
#include "nonparametric.h"
#include "utils.h"
#include "splines.h"


namespace Intsurv {

    class CoxphCureUncer {
    private:
        CoxphReg cox_obj;
        LogisticReg cure_obj;
        unsigned int cox_p;
        unsigned int cure_p;
        arma::uvec case1_ind;
        arma::uvec case2_ind;
        arma::uvec cer_ind;     // index of certain event indicators
        arma::uvec case3_ind;

    public:
        arma::vec cox_coef;
        arma::vec cure_coef;
        double negLogL;
        unsigned int nObs;        // number of observations
        unsigned int num_iter;    // number of iterations

        // the "big enough" L1 lambda => zero coef
        double cox_l1_lambda_max;
        double cure_l1_lambda_max;

        // regularized by particular lambdas
        arma::vec en_cox_coef;  // elastic net estimates
        arma::vec en_cure_coef; // elastic net estimates
        double cox_l1_lambda;
        double cox_l2_lambda;
        arma::vec cox_l1_penalty_factor;
        double cure_l1_lambda;
        double cure_l2_lambda;
        arma::vec cure_l1_penalty_factor;

        // default constructor
        CoxphCureUncer() {}

        // constructors
        CoxphCureUncer(
            const arma::vec& time,
            const arma::vec& event,
            const arma::mat& cox_x,
            const arma::mat& cure_x,
            const bool& cure_intercept = true,
            const bool& cox_standardize = true,
            const bool& cure_standardize = true
            )
        {
            // replace NA or NaN event indicator with 0.5
            // (or any number between 0 and 1)
            arma::vec event0na { event };
            const double const4na { 0.5 };
            event0na.replace(arma::datum::nan, const4na);

            // create the CoxphReg object
            this->cox_obj = CoxphReg(time, event0na, cox_x, cox_standardize);

            // pre-process x and y
            this->cox_p = cox_x.n_cols;
            this->cure_p = cure_x.n_cols +
                static_cast<unsigned int>(cure_intercept);
            this->nObs = cox_x.n_rows;
            arma::uvec cox_sort_ind { cox_obj.get_sort_index() };
            arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
            arma::vec s_event { event0na.elem(cox_sort_ind) };
            this->case1_ind = arma::find(s_event > const4na);
            this->case2_ind = arma::find(s_event < const4na);
            this->cer_ind = Intsurv::vec_union(case1_ind, case2_ind);
            this->case3_ind = arma::find(s_event == const4na);

            // create the LogisticReg object
            this->cure_obj = LogisticReg(cure_xx, s_event, cure_intercept,
                                         cure_standardize);

        }


        // function members
        // fit the Cox cure mode with uncertain events by EM algorithm
        inline void fit(
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const unsigned int& em_max_iter,
            const double& em_rel_tol,
            const unsigned int& cox_mstep_max_iter,
            const double& cox_mstep_rel_tol,
            const unsigned int& cure_mstep_max_iter,
            const double& cure_mstep_rel_tol,
            const bool& spline_start,
            const unsigned int& iSpline_num_knots,
            const unsigned int& iSpline_degree
            );

        // fit regularized Cox cure model with adaptive elastic net penalty
        // for perticular lambda's
        inline void regularized_fit(
            const double& cox_l1_lambda,
            const double& cox_l2_lambda,
            const double& cure_l1_lambda,
            const double& cure_l2_lambda,
            const arma::vec& cox_l1_penalty_factor,
            const arma::vec& cure_l1_penalty_factor,
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const unsigned int& em_max_iter,
            const double& em_rel_tol,
            const unsigned int& cox_mstep_max_iter,
            const double& cox_mstep_rel_tol,
            const unsigned int& cure_mstep_max_iter,
            const double& cure_mstep_rel_tol,
            const bool& spline_start,
            const unsigned int& iSpline_num_knots,
            const unsigned int& iSpline_degree
            );


    };                          // end of class definition


    // fit the Cox cure mode by EM algorithm
    inline void CoxphCureUncer::fit(
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int& em_max_iter = 500,
        const double& em_rel_tol = 1e-5,
        const unsigned int& cox_mstep_max_iter = 50,
        const double& cox_mstep_rel_tol = 1e-3,
        const unsigned int& cure_mstep_max_iter = 50,
        const double& cure_mstep_rel_tol = 1e-3,
        const bool& spline_start = false,
        const unsigned int& iSpline_num_knots = 3,
        const unsigned int& iSpline_degree = 2
        )
    {
        // get pre-processed design matrix, time, and event
        arma::mat cox_x { cox_obj.get_x() };
        arma::vec time { cox_obj.get_time() };
        arma::vec event { cox_obj.get_event() };

        // initialize cox_beta
        arma::vec cox_beta { arma::zeros(this->cox_p) };
        if (cox_start.n_elem == this->cox_p) {
            cox_beta = cox_start;
        } else {
            Intsurv::CoxphReg tmp_object {
                Intsurv::CoxphReg(time.elem(case1_ind),
                                  event.elem(case1_ind),
                                  cox_x.rows(case1_ind))
            };
            tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
            cox_beta = tmp_object.coef;
        }
        arma::vec cox_exp_x_beta = arma::exp(cox_obj.get_x() * cox_beta);

        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(this->cure_p) };
        if (cure_start.n_elem == this->cure_p) {
            cure_beta = cure_start;
        } else {
            cure_obj.fit(cure_beta, cure_mstep_max_iter,
                         cure_mstep_rel_tol);
            cure_beta = cure_obj.coef;
        }
        arma::vec p_vec { cure_obj.predict(cure_beta) };

        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(this->cer_ind),
                        event.elem(this->cer_ind))
        };
        cox_obj.h0_time = nelen_event.step_inst_rate(time);
        cox_obj.H0_time = nelen_event.step_cum_rate(time);
        // for censoring
        Intsurv::NelsonAalen nelen_censor {
            Intsurv::NelsonAalen(time.elem(this->cer_ind),
                                 1 - event.elem(this->cer_ind))
        };
        cox_obj.hc_time = nelen_censor.step_inst_rate(time);
        cox_obj.Hc_time = nelen_censor.step_cum_rate(time);
        // update related function values
        cox_obj.h_time = cox_obj.h0_time % cox_exp_x_beta;
        cox_obj.H_time = cox_obj.H0_time % cox_exp_x_beta;
        cox_obj.S0_time = arma::exp(- cox_obj.H0_time);
        cox_obj.S_time = arma::exp(- cox_obj.H_time);
        cox_obj.Sc_time = arma::exp(- cox_obj.Hc_time);

        // further smooth hazard and survival estimates by splines
        if (spline_start) {
            // set up spline bases
            arma::vec internal_knots_event {
                get_internal_knots(time.elem(case1_ind), iSpline_num_knots)
            };
            arma::vec internal_knots_censor {
                get_internal_knots(time.elem(case2_ind), iSpline_num_knots)
            };
            arma::vec boundary_knots { get_boundary_knots(time) };
            arma::mat iSpline_mat_event {
                iSpline(time, iSpline_degree,
                        internal_knots_event, boundary_knots)
            };
            arma::mat iSpline_mat_censor {
                iSpline(time, iSpline_degree,
                        internal_knots_censor, boundary_knots)
            };
            // fit non-negative least square
            arma::vec H0_coef {
                nnls(cox_obj.H0_time, iSpline_mat_event)
            };
            arma::vec Hc_coef {
                nnls(cox_obj.Hc_time, iSpline_mat_censor)
            };
            arma::mat mSpline_mat_event {
                mSpline(time, iSpline_degree,
                        internal_knots_event, boundary_knots)
            };
            arma::mat mSpline_mat_censor {
                mSpline(time, iSpline_degree,
                        internal_knots_censor, boundary_knots)
            };
            // update H0, Hc, h0, and hc, etc.
            cox_obj.h0_time = mSpline_mat_event * H0_coef;
            cox_obj.h_time = cox_obj.h0_time % cox_exp_x_beta;
            cox_obj.H0_time = iSpline_mat_event * H0_coef;
            cox_obj.H_time = cox_obj.H0_time % cox_exp_x_beta;
            cox_obj.S0_time = arma::exp(- cox_obj.H0_time);
            cox_obj.S_time = arma::exp(- cox_obj.H_time);
            cox_obj.hc_time = mSpline_mat_censor * Hc_coef;
            cox_obj.Hc_time = iSpline_mat_censor * Hc_coef;
            cox_obj.Sc_time = arma::exp(- cox_obj.Hc_time);
        }

        // intialization for the main loop
        size_t i {0};
        double obs_ell {0};

        // main loop of EM algorithm
        while (true) {
            arma::vec estep_m { event };
            double numer_j {0}, tol1 {0}, tol2 {0};
            double m1 {0}, m2 {0}, m3 {0}, w1 {0}, w2 {0};

            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind) {
                numer_j = p_vec(j) * cox_obj.S_time(j);
                estep_m(j) = 1 / ((1 - p_vec(j)) / numer_j + 1);
            }

            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind) {
                m1 = p_vec(j) *
                    cox_obj.h_time(j) *
                    cox_obj.S_time(j) *
                    cox_obj.Sc_time(j);
                m2 = p_vec(j) *
                    cox_obj.hc_time(j) *
                    cox_obj.Sc_time(j) *
                    cox_obj.S_time(j);
                m3 = (1 - p_vec(j)) *
                    cox_obj.hc_time(j) *
                    cox_obj.Sc_time(j);
                w1 = 1 / ((m2 + m3) / m1 + 1);
                w2 = 1 / ((m1 + m3) / m2 + 1);
                estep_m(j) = w1 + w2;
                event(j) = w1;
            }

            // M-step for the survival layer
            cox_obj.set_offset(arma::log(estep_m));
            cox_obj.update_event_weight(event);
            cox_obj.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);

            // M-step for the Cure layer
            cure_obj.update_y(estep_m);
            cure_obj.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);

            // check convergence
            tol1 = Intsurv::rel_l2_norm(cox_obj.coef, cox_beta);
            tol2 = Intsurv::rel_l2_norm(cure_obj.coef, cure_beta);

            // update to last estimates
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;

            // early exit if has any `nan`
            if (cox_beta.has_nan() || cure_beta.has_nan()) {
                obs_ell = - arma::datum::inf;
                throw std::range_error(
                    "The negative log-likelihood function went to infinite."
                    );
                break;
            }
            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // prepare the E-step for next iteration
            cox_obj.compute_haz_surv_time();
            cox_obj.compute_censor_haz_surv_time();
            p_vec = cure_obj.predict(cure_beta);

            // check convergence
            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i > em_max_iter) {
                // compute observed data likelihood
                // for case 1
                for (size_t j: case1_ind) {
                    obs_ell += std::log(p_vec(j)) +
                        std::log(cox_obj.h_time(j)) +
                        std::log(cox_obj.S_time(j)) +
                        std::log(cox_obj.Sc_time(j));
                }
                // for case 2
                for (size_t j: case2_ind) {
                    obs_ell +=
                        std::log(p_vec(j) * cox_obj.S_time(j) + 1 - p_vec(j)) +
                        std::log(cox_obj.Sc_time(j)) +
                        std::log(cox_obj.hc_time(j));
                }
                // for case 3
                for (size_t j: case3_ind) {
                    m1 = p_vec(j) *
                        cox_obj.h_time(j) * cox_obj.S_time(j) *
                        cox_obj.Sc_time(j);
                    m2 = p_vec(j) *
                        cox_obj.hc_time(j) * cox_obj.Sc_time(j) *
                        cox_obj.S_time(j);
                    m3 = (1 - p_vec(j)) *
                        cox_obj.hc_time(j) * cox_obj.Sc_time(j);
                    obs_ell += std::log(m1 + m2 + m3);
                }
                break;
            }
            // update iter
            ++i;
        } // end of the EM algorithm
        // prepare outputs
        this->cox_coef = cox_beta;
        this->cure_coef = cure_beta;
        this->negLogL = - obs_ell;
        this->num_iter = i;
    }


}


#endif
