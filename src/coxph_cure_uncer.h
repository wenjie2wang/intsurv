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
#include "assessment.h"
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
        // number of predictors in survival layer
        unsigned int cox_p;
        // number of predictors in cure layer with possible intercept
        unsigned int cure_p;
        // number of predictors in cure layer without possible intercept
        unsigned int cure_p0;
        arma::uvec case1_ind;
        arma::uvec case2_ind;
        arma::uvec cer_ind;     // index of rows with correct event indicators
        arma::uvec case3_ind;
        unsigned int max_event_time_ind; // index of the maximum event time

    public:
        arma::vec cox_coef;
        arma::vec cure_coef;
        unsigned int coef_df;     // degree of freedom of coef estimates
        double negLogL;           // negative log-likelihood
        unsigned int nObs;        // number of observations
        unsigned int num_iter;    // number of iterations
        // BIC: log(num_obs) * coef_df + 2 * negLogL
        double bic1;
        // BIC: log(num_certain_event) * coef_df + 2 * negLogL
        double bic2;
        double c_index;           // weighted C-index

        // for each subject and in the original order of X
        arma::vec cox_xBeta;        // score from the survival layer
        arma::vec cure_xBeta;       // score from the cure layer
        arma::vec susceptible_prob; // probability of being susceptible
        // values in the last E-step
        arma::vec estep_cured;
        arma::vec estep_event;
        arma::vec estep_censor;

        // hazard and survival function estimates at unique time
        arma::vec unique_time;
        arma::vec h0_est;
        arma::vec H0_est;
        arma::vec S0_est;
        arma::vec hc_est;
        arma::vec Hc_est;
        arma::vec Sc_est;

        // the "big enough" L1 lambda => zero coef
        double cox_l1_lambda_max;
        double cure_l1_lambda_max;

        // regularized by particular lambdas
        arma::vec cox_en_coef;  // elastic net estimates
        arma::vec cure_en_coef; // elastic net estimates
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
            this->cure_p0 = cure_x.n_cols;
            this->cure_p = this->cure_p0 +
                static_cast<unsigned int>(cure_intercept);
            this->nObs = cox_x.n_rows;
            arma::uvec cox_sort_ind { cox_obj.get_sort_index() };
            arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
            arma::vec s_event { event0na.elem(cox_sort_ind) };
            this->case1_ind = arma::find(s_event > const4na);
            this->case2_ind = arma::find(s_event < const4na);
            this->cer_ind = vec_union(case1_ind, case2_ind);
            this->case3_ind = arma::find(s_event == const4na);
            this->max_event_time_ind = arma::max(this->case1_ind);
            // create the LogisticReg object
            this->cure_obj = LogisticReg(cure_xx, s_event, cure_intercept,
                                         cure_standardize);

        }


        // function members
        // helper functions
        inline unsigned int get_cox_p() {
            return this->cox_p;
        }
        inline unsigned int get_cure_p() {
            return this->cure_p;
        }

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
            const unsigned int& iSpline_degree,
            const unsigned int& tail_completion,
            double tail_tau,
            const double& pmin,
            const unsigned int& early_stop,
            const unsigned int& verbose
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
            const unsigned int& iSpline_degree,
            const unsigned int& tail_completion,
            double tail_tau,
            const double& pmin,
            const unsigned int& early_stop,
            const unsigned int& verbose
            );

        // compute BIC
        inline void compute_bic1() {
            this->bic1 = std::log(nObs) * coef_df + 2 * negLogL;
        }
        inline void compute_bic2() {
            this->bic2 = std::log(case1_ind.n_elem) *
                coef_df + 2 * negLogL;
        }

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
        const unsigned int& iSpline_degree = 2,
        const unsigned int& tail_completion = 1,
        double tail_tau = -1,
        const double& pmin = 1e-5,
        const unsigned int& early_stop = 0,
        const unsigned int& verbose = 0
        )
    {
        // get pre-processed design matrix, time, and event
        const arma::mat& cox_x { cox_obj.get_x() };
        const arma::vec& time { cox_obj.get_time() };
        arma::vec event { cox_obj.get_event() };

        // initialize cox_beta
        arma::vec cox_beta { arma::zeros(this->cox_p) };
        if (cox_start.n_elem == this->cox_p) {
            cox_beta = cox_start;
        } else {
            CoxphReg tmp_object {
                CoxphReg(time.elem(case1_ind),
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

        // initialize coef estimates
        cure_obj.coef = cure_beta;
        cox_obj.coef = cox_beta;

        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(this->cer_ind),
                        event.elem(this->cer_ind))
        };
        cox_obj.h0_time = nelen_event.step_inst_rate(time);
        cox_obj.H0_time = nelen_event.step_cum_rate(time);
        // for censoring
        NelsonAalen nelen_censor {
            NelsonAalen(time.elem(this->cer_ind),
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
        arma::vec p_vec { arma::zeros(nObs) };
        arma::vec estep_m { event };
        size_t i {0};
        double obs_ell {0}, obs_ell_old { - arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;

        // prepare for tail completion
        const arma::uvec case23_ind { vec_union(case2_ind, case3_ind) };
        double max_event_time { time(this->max_event_time_ind) };
        if (tail_tau < 0)
            tail_tau = arma::datum::inf;

        // allow users to stop here
        Rcpp::checkUserInterrupt();

        // main loop of EM algorithm
        while (true) {

            // update to the latest estimates
            p_vec = cure_obj.predict(cure_obj.coef, pmin);
            cox_obj.compute_haz_surv_time();
            cox_obj.compute_censor_haz_surv_time();

            // prepare for exponential tail completion method
            double s0_tau {0}, etail_lambda {0};
            if (tail_completion == 2) {
                s0_tau = cox_obj.S0_time(max_event_time_ind);
                etail_lambda = - std::log(s0_tau / max_event_time);
            }
            // tail completion for the conditional survival function
            // for case 2 and case 3
            for (size_t j: case23_ind) {
                // tail completion for the conditional survival function
                switch(tail_completion) {
                    case 0:
                        // tail completion after the given tail_tau
                        // by default, it means no tail completion
                        if (time(j) > tail_tau) {
                            cox_obj.S_time(j) = 0;
                            cox_obj.S0_time(j) = 0;
                        }
                        break;
                    case 1:
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            cox_obj.S_time(j) = 0;
                            cox_obj.S0_time(j) = 0;
                        }
                        break;
                    case 2:
                        // exponential tail by Peng (2003)
                        if (time(j) > max_event_time) {
                            cox_obj.S0_time(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            cox_obj.S_time(j) = std::pow(
                                cox_obj.S0_time(j),
                                std::exp(cox_obj.xBeta(j))
                                );
                        }
                        break;
                    default:    // do nothing, otherwise
                        break;
                }
            }

            // compute observed data log-likelihood
            obs_ell = 0;
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
                double m12_common {
                    p_vec(j) * cox_obj.S0_time(j) * cox_obj.Sc_time(j)
                };
                double m1 { cox_obj.h_time(j) * m12_common };
                double m2 { cox_obj.hc_time(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * cox_obj.hc_time(j) * cox_obj.Sc_time(j)
                };
                obs_ell += std::log(m1 + m2 + m3);
            }

            // if verbose
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(50, '=')
                            << "\niteration: " << i
                            << "\n  Cox coef: "
                            << arma2rvec(cox_obj.coef)
                            << "\n    relative diff: " << tol1
                            << "\n  cure coef: "
                            << arma2rvec(cure_obj.coef)
                            << "\n    relative diff: " << tol2
                            << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << std::endl;
            }

            bool early_exit { false };
            // early exit if has any `nan`
            if (cox_obj.coef.has_nan() || cure_obj.coef.has_nan()) {
                obs_ell = - arma::datum::inf;
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function went to infinite."
                            << std::endl;
                early_exit = true;
            }
            // early exit if the regularized objective function increased, which
            // is technically impossible and thus can serve as a warning
            if (obs_ell < obs_ell_old) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "The observed data log-likelihood decreased."
                                << std::endl;
                }
                early_exit = early_exit || early_stop;
            }
            // return the estimates from last step
            if (early_exit) {
                if (verbose) {
                    Rcpp::Rcout << "Ended the EM algorithm after iteration "
                                << i
                                << " with estimates from last step."
                                << std::endl;
                }
                // take the estimates from the last step
                cox_obj.coef = cox_beta;
                cure_obj.coef = cure_beta;
                // update hazard and survival function estimates
                cox_obj.compute_haz_surv_time();
                cox_obj.S0_time = s0_wi_tail;
                cox_obj.S_time = s_wi_tail;
                cox_obj.compute_censor_haz_surv_time();
                cox_obj.est_haz_surv();
                // break here
                break;
            }

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i >= em_max_iter) {
                if (verbose) {
                    if (i < em_max_iter) {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached convergence after " << i
                                    << " iterations\n" << std::endl;
                    } else {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached the max iteration number."
                                    << std::endl;
                    }
                }
                // compute hazard and survival function estimates
                cox_obj.est_haz_surv();
                // get out of the loop here
                break;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // record estimates from last step
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;
            obs_ell_old = obs_ell;
            s0_wi_tail = cox_obj.S0_time;
            s_wi_tail = cox_obj.S_time;

            // update iter for the next iteration
            ++i;

            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind) {
                double numer_j { p_vec(j) * cox_obj.S_time(j)};
                estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
            }

            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind) {
                double m12_common {
                    p_vec(j) * cox_obj.S_time(j) * cox_obj.Sc_time(j)
                };
                double m1 { cox_obj.h_time(j) * m12_common };
                double m2 { cox_obj.hc_time(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * cox_obj.hc_time(j) * cox_obj.Sc_time(j)
                };
                double m { m1 + m2 + m3 };
                double w1 { m1 / m };
                double w2 { m2 / m };
                // some special care for subjects in case 3
                // since event(j) cannot be either 0 or 1!
                w1 = ge_le(w1, pmin, 1 - pmin);
                w2 = ge_le(w2, pmin, 1 - pmin);
                estep_m(j) = w1 + w2;
                event(j) = w1;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning the M-step for the survival layer:"
                            << std::endl;
            }
            cox_obj.set_offset(arma::log(estep_m));
            cox_obj.update_event_weight(event);
            cox_obj.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol,
                        early_stop, verbose > 2);
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nThe M-step for the survival layer was done."
                            << std::endl;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the Cure layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning M-step for the cure layer:";
            }
            cure_obj.update_y(estep_m);
            cure_obj.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol,
                         pmin, early_stop, verbose > 2);
            if (verbose > 1) {
                Rcpp::Rcout << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done."
                            << std::endl;
            }

            // update tolerance
            tol1 = rel_l2_norm(cox_obj.coef, cox_beta);
            tol2 = rel_l2_norm(cure_obj.coef, cure_beta);

        } // end of the EM algorithm

        // reset cox_obj and cure_obj in case of further usage
        cox_obj.reset_offset();
        cure_obj.update_y(cox_obj.get_event());

        // prepare outputs
        this->cox_coef = cox_beta;
        this->cure_coef = cure_beta;
        this->unique_time = cox_obj.unique_time;
        this->h0_est = cox_obj.h0_est;
        this->H0_est = cox_obj.H0_est;
        this->S0_est = cox_obj.S0_est;
        this->hc_est = cox_obj.hc_est;
        this->Hc_est = cox_obj.Hc_est;
        this->Sc_est = cox_obj.Sc_est;
        this->negLogL = - obs_ell;
        this->coef_df = cox_obj.coef_df + cure_obj.coef_df;
        this->num_iter = i;
        this->compute_bic1();
        this->compute_bic2();

        // prepare scores and prob in their original order
        arma::uvec rev_ord { cox_obj.get_rev_sort_index() };
        this->cox_xBeta = cox_obj.xBeta.elem(rev_ord);
        this->cure_xBeta = cure_obj.xBeta.elem(rev_ord);
        this->susceptible_prob = cure_obj.prob_vec.elem(rev_ord);

        // compute posterior probabilities from E-step
        for (size_t j: case2_ind) {
            double numer_j { p_vec(j) *  cox_obj.S_time(j)};
            estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        for (size_t j: case3_ind) {
            double m12_common {
                p_vec(j) * cox_obj.S_time(j) * cox_obj.Sc_time(j)
            };
            double m1 { cox_obj.h_time(j) * m12_common };
            double m2 { cox_obj.hc_time(j) * m12_common };
            double m3 {
                (1 - p_vec(j)) * cox_obj.hc_time(j) * cox_obj.Sc_time(j)
            };
            double m { m1 + m2 + m3 };
            double w1 { m1 / m };
            double w2 { m2 / m };
            w1 = ge_le(w1, pmin, 1 - pmin);
            w2 = ge_le(w2, pmin, 1 - pmin);
            estep_m(j) = w1 + w2;
            event(j) = w1;
        }
        this->estep_cured = 1 - estep_m.elem(rev_ord);
        this->estep_event = event(rev_ord);
        this->estep_censor = 1 - this->estep_cured - this->estep_event;

        // compute weighted c-index over the certain records
        this->c_index =
            Intsurv::Concordance(cox_obj.get_time().elem(cer_ind),
                                 cox_obj.get_event().elem(cer_ind),
                                 cox_obj.xBeta.elem(cer_ind),
                                 cure_obj.prob_vec.elem(cer_ind)).index;
    }


    // fit regularized Cox cure model with uncertain events
    // with adaptive elastic net penalty for perticular lambda's
    // lambda_1 * lasso * factors + lambda_2 * ridge
    inline void CoxphCureUncer::regularized_fit(
        const double& cox_l1_lambda = 0,
        const double& cox_l2_lambda = 0,
        const double& cure_l1_lambda = 0,
        const double& cure_l2_lambda = 0,
        const arma::vec& cox_l1_penalty_factor = 0,
        const arma::vec& cure_l1_penalty_factor = 0,
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int& em_max_iter = 500,
        const double& em_rel_tol = 1e-5,
        const unsigned int& cox_mstep_max_iter = 200,
        const double& cox_mstep_rel_tol = 1e-5,
        const unsigned int& cure_mstep_max_iter = 200,
        const double& cure_mstep_rel_tol = 1e-5,
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
        // get pre-processed design matrix, time, and event
        const arma::mat& cox_x { cox_obj.get_x() };
        const arma::vec& time { cox_obj.get_time() };
        arma::vec event { cox_obj.get_event() };

        // L1 penalty factor for Cox model
        arma::vec cox_l1_penalty { arma::ones(this->cox_p) };
        if (cox_l1_penalty_factor.n_elem == cox_p) {
            // re-scale so that sum(factor) = number of predictors
            cox_l1_penalty = cox_l1_penalty_factor * cox_p /
                arma::sum(cox_l1_penalty_factor);
        }
        this->cox_l1_penalty_factor = cox_l1_penalty;
        // L1 penalty factor for cure model
        arma::vec cure_l1_penalty { arma::ones(this->cure_p0) };
        if (cure_l1_penalty_factor.n_elem == this->cure_p0) {
            // re-scale so that sum(factor) = number of predictors
            cure_l1_penalty = cure_l1_penalty_factor * cure_p0 /
                arma::sum(cure_l1_penalty_factor);
        }
        this->cure_l1_penalty_factor = cure_l1_penalty;

        // initialized with all zeros coef
        arma::vec cox_beta { arma::zeros(this->cox_p) };
        arma::vec cure_beta { arma::zeros(this->cure_p) };

        // compute the large enough lambdas that result in all-zero estimates
        arma::vec cox_grad_zero { arma::abs(cox_obj.gradient(cox_beta)) };
        this->cox_l1_lambda_max =
            arma::max(cox_grad_zero / cox_l1_penalty) / this->nObs;
        arma::vec cure_grad_zero {
            arma::abs(cure_obj.gradient(cure_beta, pmin))
        };
        this->cure_l1_lambda_max =
            arma::max(cure_grad_zero.tail(cure_l1_penalty.n_elem) /
                      cure_l1_penalty) / this->nObs;

        // early stop if we want lambda_max by setting em_max_iter = 0
        if (em_max_iter == 0) {
            return;
        }

        // initialize cox_beta
        if (cox_start.n_elem == this->cox_p) {
            cox_beta = cox_start;
        } else {
            cox_beta = arma::zeros(this->cox_p);
        }
        arma::vec cox_exp_x_beta = arma::exp(cox_obj.get_x() * cox_beta);
        // initialize cure_beta
        if (cure_start.n_elem == this->cure_p) {
            cure_beta = cure_start;
        } else {
            cure_beta = arma::zeros(this->cure_p);
        }

        // initialize coef estimates
        cure_obj.coef = cure_beta;
        cox_obj.coef = cox_beta;

        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(this->cer_ind),
                        event.elem(this->cer_ind))
        };
        cox_obj.h0_time = nelen_event.step_inst_rate(time);
        cox_obj.H0_time = nelen_event.step_cum_rate(time);
        // for censoring
        NelsonAalen nelen_censor {
            NelsonAalen(time.elem(this->cer_ind),
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
        arma::vec p_vec { arma::zeros(nObs) };
        arma::vec estep_m { event };
        size_t i {0};
        double obs_ell {0};
        double reg_obj {0}, reg_obj_old { arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;

        // prepare for tail completion
        const arma::uvec case23_ind { vec_union(case2_ind, case3_ind) };
        double max_event_time { time(this->max_event_time_ind) };
        if (tail_tau < 0)
            tail_tau = arma::datum::inf;

        // allow users to stop here
        Rcpp::checkUserInterrupt();

        // main loop of EM algorithm
        while (true) {

            // update to the latest estimates
            p_vec = cure_obj.predict(cure_obj.coef, pmin);
            cox_obj.compute_haz_surv_time();
            cox_obj.compute_censor_haz_surv_time();

            // prepare for exponential tail completion method
            double s0_tau {0}, etail_lambda {0};
            if (tail_completion == 2) {
                s0_tau = cox_obj.S0_time(max_event_time_ind);
                etail_lambda = - std::log(s0_tau / max_event_time);
            }
            // tail completion for the conditional survival function
            // for case 2 and case 3
            for (size_t j: case23_ind) {
                switch(tail_completion) {
                    case 0:
                        // tail completion after the given tail_tau
                        // by default, it means no tail completion
                        if (time(j) > tail_tau) {
                            cox_obj.S_time(j) = 0;
                            cox_obj.S0_time(j) = 0;
                        }
                        break;
                    case 1:
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            cox_obj.S_time(j) = 0;
                            cox_obj.S0_time(j) = 0;
                        }
                        break;
                    case 2:
                        // exponential tail by Peng (2003)
                        if (time(j) > max_event_time) {
                            cox_obj.S0_time(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            cox_obj.S_time(j) = std::pow(
                                cox_obj.S0_time(j),
                                std::exp(cox_obj.xBeta(j))
                                );
                        }
                        break;
                    default:    // do nothing, otherwise
                        break;
                }
            }

            // compute observed data log-likelihood
            obs_ell = 0;
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
                double m1 {
                    p_vec(j) * cox_obj.h_time(j) *
                    cox_obj.S_time(j) * cox_obj.Sc_time(j)
                };
                double m2 {
                    p_vec(j) * cox_obj.hc_time(j) *
                    cox_obj.Sc_time(j) * cox_obj.S_time(j)
                };
                double m3 {
                    (1 - p_vec(j)) * cox_obj.hc_time(j) * cox_obj.Sc_time(j)
                };
                obs_ell += std::log(m1 + m2 + m3);
            }
            // compuete the regularized objective function
            double reg_cox {
                cox_l1_lambda * l1_norm(cox_obj.coef % cox_l1_penalty) +
                cox_l2_lambda * sum_of_square(cox_obj.coef)
            };
            double reg_cure {
                cure_l1_lambda *
                l1_norm(cure_obj.coef.tail(cure_p0) % cure_l1_penalty) +
                cure_l2_lambda *
                sum_of_square(cure_obj.coef.tail(cure_p0))
            };
            reg_obj = - obs_ell / this->nObs + reg_cox + reg_cure;

            // if verbose
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(50, '=')
                            << "\niteration: " << i
                            << "\n  Cox coef: "
                            << arma2rvec(cox_obj.coef)
                            << "\n    relative diff: " << tol1
                            << "\n  cure coef: "
                            << arma2rvec(cure_obj.coef)
                            << "\n    relative diff: " << tol2
                            << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << "\n  regularized objective function: "
                            << reg_obj
                            << "\n    penalty on Cox model: "
                            << reg_cox
                            << "\n    penalty on cure layer: "
                            << reg_cure
                            << std::endl;
            }

            bool early_exit { false };
            // early exit if has any `nan`
            if (cox_obj.coef.has_nan() || cure_obj.coef.has_nan()) {
                obs_ell = - arma::datum::inf;
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function went to infinite."
                            << std::endl;
                early_exit = true;
            }
            // early exit if the regularized objective function increased, which
            // is technically impossible and thus can serve as a warning
            if (reg_obj > reg_obj_old) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "The observed data log-likelihood decreased."
                                << std::endl;
                }
                early_exit = early_exit || early_stop;
            }
            // return the estimates from last step
            if (early_exit) {
                if (verbose) {
                    Rcpp::Rcout << "Ended the EM algorithm after iteration "
                                << i
                                << " with estimates from last step."
                                << std::endl;
                }
                // take the estimates from the last step
                cox_obj.coef = cox_beta;
                cure_obj.coef = cure_beta;
                // update hazard and survival function estimates
                cox_obj.compute_haz_surv_time();
                cox_obj.S0_time = s0_wi_tail;
                cox_obj.S_time = s_wi_tail;
                cox_obj.compute_censor_haz_surv_time();
                cox_obj.est_haz_surv();
                // break here
                break;
            }

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i >= em_max_iter) {
                if (verbose) {
                    if (i < em_max_iter) {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached convergence after " << i
                                    << " iterations\n" << std::endl;
                    } else {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached the max iteration number."
                                    << std::endl;
                    }
                }
                // compute hazard and survival function estimates
                cox_obj.est_haz_surv();
                // get out of the loop here
                break;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // record estimates from last step
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;
            reg_obj_old = reg_obj;
            s0_wi_tail = cox_obj.S0_time;
            s_wi_tail = cox_obj.S_time;

            // update iter for the next iteration
            ++i;

            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind) {
                double numer_j { p_vec(j) * cox_obj.S_time(j)};
                estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
            }

            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind) {
                double m12_common {
                    p_vec(j) * cox_obj.S0_time(j) * cox_obj.Sc_time(j)
                };
                double m1 { cox_obj.h_time(j) * m12_common };
                double m2 { cox_obj.hc_time(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * cox_obj.hc_time(j) * cox_obj.Sc_time(j)
                };
                double m { m1 + m2 + m3 };
                double w1 { m1 / m };
                double w2 { m2 / m };
                // some special care for subjects in case 3
                // since event(j) cannot be either 0 or 1!
                w1 = ge_le(w1, pmin, 1 - pmin);
                w2 = ge_le(w2, pmin, 1 - pmin);
                estep_m(j) = w1 + w2;
                event(j) = w1;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning the M-step for the survival layer:"
                            << std::endl;
            }
            cox_obj.set_offset(arma::log(estep_m));
            cox_obj.update_event_weight(event);
            cox_obj.regularized_fit(
                cox_l1_lambda, cox_l2_lambda, cox_l1_penalty_factor,
                cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol,
                early_stop, verbose > 2
                );
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nThe M-step for the survival layer was done."
                            << std::endl;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the Cure layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning M-step for the cure layer:";
            }
            cure_obj.update_y(estep_m);
            cure_obj.regularized_fit(
                cure_l1_lambda, cure_l2_lambda, cure_l1_penalty_factor,
                cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol,
                pmin, early_stop, verbose > 2
                );
            if (verbose > 1) {
                Rcpp::Rcout << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done."
                            << std::endl;
            }

            // update tolerance
            tol1 = rel_l2_norm(cox_obj.coef, cox_beta);
            tol2 = rel_l2_norm(cure_obj.coef, cure_beta);

        } // end of the EM algorithm

        // reset cox_obj and cure_obj in case of further usage
        cox_obj.reset_offset();
        cure_obj.update_y(cox_obj.get_event());

        // prepare outputs
        this->cox_coef = cox_obj.coef;
        this->cure_coef = cure_obj.coef;
        this->cox_en_coef = cox_obj.en_coef;
        this->cure_en_coef = cure_obj.en_coef;

        this->unique_time = cox_obj.unique_time;
        this->h0_est = cox_obj.h0_est;
        this->H0_est = cox_obj.H0_est;
        this->S0_est = cox_obj.S0_est;
        this->hc_est = cox_obj.hc_est;
        this->Hc_est = cox_obj.Hc_est;
        this->Sc_est = cox_obj.Sc_est;

        this->negLogL = - obs_ell;
        this->coef_df = cox_obj.coef_df + cure_obj.coef_df;
        this->cox_l1_lambda = cox_l1_lambda;
        this->cox_l2_lambda = cox_l2_lambda;
        this->cure_l1_lambda = cure_l1_lambda;
        this->cure_l2_lambda = cure_l2_lambda;
        this->num_iter = i;
        this->compute_bic1();
        this->compute_bic2();

        // prepare scores and prob in their original order
        arma::uvec rev_ord { cox_obj.get_rev_sort_index() };
        this->cox_xBeta = cox_obj.xBeta.elem(rev_ord);
        this->cure_xBeta = cure_obj.xBeta.elem(rev_ord);
        this->susceptible_prob = cure_obj.prob_vec.elem(rev_ord);

        // compute posterior probabilities from E-step
        for (size_t j: case2_ind) {
            double numer_j { p_vec(j) *  cox_obj.S_time(j)};
            estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        for (size_t j: case3_ind) {
            double m12_common {
                p_vec(j) * cox_obj.S0_time(j) * cox_obj.Sc_time(j)
            };
            double m1 { cox_obj.h_time(j) * m12_common };
            double m2 { cox_obj.hc_time(j) * m12_common };
            double m3 {
                (1 - p_vec(j)) * cox_obj.hc_time(j) * cox_obj.Sc_time(j)
            };
            double m { m1 + m2 + m3 };
            double w1 { m1 / m };
            double w2 { m2 / m };
            estep_m(j) = w1 + w2;
            event(j) = w1;
        }
        this->estep_cured = 1 - estep_m.elem(rev_ord);
        this->estep_event = event(rev_ord);
        this->estep_censor = 1 - this->estep_cured - this->estep_event;

        // compute weighted c-index over the certain records
        this->c_index =
            Intsurv::Concordance(cox_obj.get_time().elem(cer_ind),
                                 cox_obj.get_event().elem(cer_ind),
                                 cox_obj.xBeta.elem(cer_ind),
                                 cure_obj.prob_vec.elem(cer_ind)).index;
    }


}


#endif
