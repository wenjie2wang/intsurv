//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2021  Wenjie Wang <wang@wwenjie.org>
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

#ifndef INTSURV_COXPH_CURE_MCAR_H
#define INTSURV_COXPH_CURE_MCAR_H

#include <RcppArmadillo.h>
#include "assessment.h"
#include "coxph_reg.h"
#include "logistic_reg.h"
#include "nonparametric.h"
#include "utils.h"


namespace Intsurv {

    class CoxphCureMcar {
    public:
        CoxphReg cox_obj_;
        LogisticReg cure_obj_;
        // number of predictors in survival layer
        unsigned int cox_p_;
        // number of predictors in cure layer with possible intercept
        unsigned int cure_p_;
        // number of predictors in cure layer without possible intercept
        unsigned int cure_p0_;
        arma::uvec case1_ind_;
        arma::uvec case2_ind_;
        arma::uvec cer_ind_;    // index of rows with correct event indicators
        arma::uvec case3_ind_;
        unsigned int max_event_time_ind_; // index of the maximum event time
        double pmin_;
        bool cox_standardize_;
        bool cure_standardize_;

        // outputs
        arma::vec cox_coef_;
        arma::vec cure_coef_;
        unsigned int coef_df_;  // degree of freedom of coef estimates
        double neg_ll_;         // negative log-likelihood
        unsigned int n_obs_;    // number of observations
        double dn_obs_;         // double version of n_obs_
        unsigned int n_event_;  // number of certain events
        unsigned int n_iter_;   // number of iterations
        // BIC: log(num_obs) * coef_df_ + 2 * neg_ll_
        double bic1_;
        // BIC: log(num_certain_event) * coef_df_ + 2 * neg_ll_
        double bic2_;
        double aic_;            // AIC: 2 * coef_df_ + 2 * neg_ll_
        double c_index_;        // weighted C-index

        // for each subject and in the original order of X
        arma::vec cox_xbeta_;        // score from the survival layer
        arma::vec cure_xbeta_;       // score from the cure layer
        arma::vec susceptible_prob_; // probability of being susceptible
        // values in the last E-step
        arma::vec estep_cured_;
        arma::vec estep_event_;
        arma::vec estep_censor_;

        // hazard and survival function estimates at unique time
        arma::vec unique_time_;
        arma::vec h0_est_;
        arma::vec H0_est_;
        arma::vec S0_est_;
        arma::vec hc_est_;
        arma::vec Hc_est_;
        arma::vec Sc_est_;

        // the "big enough" L1 lambda => zero coef
        double cox_l1_lambda_max_;
        double cure_l1_lambda_max_;

        // regularized by particular lambdas
        // arma::vec cox_en_coef_;  // elastic net estimates
        // arma::vec cure_en_coef_; // elastic net estimates
        double cox_l1_lambda_;
        double cox_l2_lambda_;
        arma::vec cox_l1_penalty_factor_;
        double cure_l1_lambda_;
        double cure_l2_lambda_;
        arma::vec cure_l1_penalty_factor_;

        // default constructor
        CoxphCureMcar() {}

        // constructors
        CoxphCureMcar(
            const arma::vec& time,
            const arma::vec& event,
            const arma::mat& cox_x,
            const arma::mat& cure_x,
            const bool cure_intercept = true,
            const bool cox_standardize = true,
            const bool cure_standardize = true,
            const arma::vec& cox_offset = 0,
            const arma::vec& cure_offset = 0
            )
        {
            // replace NA or NaN event indicator with 0.5
            // (or any number between 0 and 1)
            arma::vec event0na { event };
            const double const4na { 0.5 };
            event0na.replace(arma::datum::nan, const4na);

            // create the CoxphReg object
            cox_obj_ = CoxphReg(time, event0na, cox_x, cox_standardize);
            cox_obj_.set_offset(cox_offset, false);
            // pre-process x and y
            cox_p_ = cox_x.n_cols;
            cure_p0_ = cure_x.n_cols;
            cure_p_ = cure_p0_ +
                static_cast<unsigned int>(cure_intercept);
            n_obs_ = cox_x.n_rows;
            dn_obs_ = static_cast<double>(n_obs_);
            arma::uvec cox_sort_ind { cox_obj_.ord_ };
            arma::vec s_event { event0na.elem(cox_sort_ind) };
            // initialize offset terms
            arma::vec s_cure_offset;
            if (cure_offset.n_elem == 1 || cure_offset.empty()) {
                s_cure_offset = arma::zeros(n_obs_);
            } else if (cure_offset.n_elem == n_obs_) {
                s_cure_offset = cure_offset.elem(cox_sort_ind);
            } else {
                throw std::length_error(
                    "The length of offset must match sample size.");
            }
            case1_ind_ = arma::find(s_event > const4na);
            case2_ind_ = arma::find(s_event < const4na);
            n_event_ = case1_ind_.n_elem;
            cer_ind_ = vec_union(case1_ind_, case2_ind_);
            case3_ind_ = arma::find(s_event == const4na);
            max_event_time_ind_ = arma::max(case1_ind_);
            // create the LogisticReg object
            cure_obj_ = LogisticReg(cure_x.rows(cox_sort_ind),
                                    s_event, cure_intercept,
                                    cure_standardize);
            cure_obj_.set_offset(s_cure_offset);
            // avoid standardization after each iteration
            cox_standardize_ = cox_standardize;
            cox_obj_.standardize_ = false;
            cure_standardize_ = cure_standardize;
            cure_obj_.standardize_ = false;
        }


        // function members
        // fit the Cox cure model with uncertain events by EM algorithm
        inline void fit(
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const unsigned int em_max_iter,
            const double em_rel_tol,
            const unsigned int cox_mstep_max_iter,
            const double cox_mstep_rel_tol,
            const unsigned int cure_mstep_max_iter,
            const double cure_mstep_rel_tol,
            const unsigned int tail_completion,
            double tail_tau,
            const double pmin,
            const unsigned int early_stop,
            const unsigned int verbose
            );

        // fit regularized Cox cure model with adaptive elastic net penalty
        // for perticular lambda's
        inline void regularized_fit(
            const double cox_l1_lambda,
            const double cox_l2_lambda,
            const double cure_l1_lambda,
            const double cure_l2_lambda,
            const arma::vec& cox_l1_penalty_factor,
            const arma::vec& cure_l1_penalty_factor,
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const unsigned int em_max_iter,
            const double em_rel_tol,
            const unsigned int cox_mstep_max_iter,
            const double cox_mstep_rel_tol,
            const unsigned int cure_mstep_max_iter,
            const double cure_mstep_rel_tol,
            const unsigned int tail_completion,
            double tail_tau,
            const double pmin,
            const unsigned int early_stop,
            const unsigned int verbose
            );

        // function to compute the observe data log-likelihood function
        // for given fitted model and estimates
        inline double obs_log_likelihood() const;

        // for given fitted model and a new set of data
        // all the inputs will be sorted inside of the function
        // so their copies are asked here
        inline double obs_log_likelihood(
            arma::vec new_time,
            arma::vec new_event,
            arma::mat new_cox_x,
            arma::mat new_cure_x,
            arma::vec new_cox_offset,
            arma::vec new_cure_offset
            ) const;

        // compute BIC
        inline void compute_bic1() {
            bic1_ = std::log(n_obs_) * coef_df_ + 2 * neg_ll_;
        }
        inline void compute_bic2() {
            bic2_ = std::log(case1_ind_.n_elem) *
                coef_df_ + 2 * neg_ll_;
        }
        inline void compute_aic() {
            aic_ = 2 * (coef_df_ + neg_ll_);
        }

    };                          // end of class definition


    // fit the Cox cure mode by EM algorithm
    inline void CoxphCureMcar::fit(
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int em_max_iter = 500,
        const double em_rel_tol = 1e-5,
        const unsigned int cox_mstep_max_iter = 50,
        const double cox_mstep_rel_tol = 1e-3,
        const unsigned int cure_mstep_max_iter = 50,
        const double cure_mstep_rel_tol = 1e-3,
        const unsigned int tail_completion = 1,
        double tail_tau = -1,
        const double pmin = 1e-5,
        const unsigned int early_stop = 0,
        const unsigned int verbose = 0
        )
    {
        pmin_ = pmin;
        // get pre-processed design matrix, time, and event
        const arma::vec& time { cox_obj_.time_ };
        arma::vec event { cox_obj_.event_ };

        // initialize cox_beta
        arma::vec cox_beta { arma::zeros(cox_p_) };
        if (cox_start.n_elem == cox_p_) {
            cox_beta = cox_start;
        } else {
            CoxphReg tmp_object {
                time.elem(case1_ind_),
                event.elem(case1_ind_),
                cox_obj_.x_.rows(case1_ind_)
            };
            tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
            cox_beta = tmp_object.coef_;
        }
        arma::vec cox_exp_x_beta = arma::exp(cox_obj_.x_ * cox_beta);

        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(cure_p_) };
        if (cure_start.n_elem == cure_p_) {
            cure_beta = cure_start;
        } else {
            cure_obj_.fit(cure_beta, cure_mstep_max_iter,
                          cure_mstep_rel_tol);
            cure_beta = cure_obj_.coef_;
        }

        // initialize coef estimates
        cure_obj_.coef_ = cure_beta;
        cox_obj_.coef_ = cox_beta;

        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(cer_ind_),
                        event.elem(cer_ind_))
        };
        cox_obj_.h0_time_ = nelen_event.step_inst_rate(time);
        cox_obj_.H0_time_ = nelen_event.step_cum_rate(time);
        // for censoring
        NelsonAalen nelen_censor {
            NelsonAalen(time.elem(cer_ind_),
                        1 - event.elem(cer_ind_))
        };
        cox_obj_.hc_time_ = nelen_censor.step_inst_rate(time);
        cox_obj_.Hc_time_ = nelen_censor.step_cum_rate(time);
        // update related function values
        cox_obj_.h_time_ = cox_obj_.h0_time_ % cox_exp_x_beta;
        cox_obj_.H_time_ = cox_obj_.H0_time_ % cox_exp_x_beta;
        cox_obj_.S0_time_ = arma::exp(- cox_obj_.H0_time_);
        cox_obj_.S_time_ = arma::exp(- cox_obj_.H_time_);
        cox_obj_.Sc_time_ = arma::exp(- cox_obj_.Hc_time_);

        // intialization for the main loop
        arma::vec p_vec { arma::zeros(n_obs_) };
        arma::vec estep_m { event };
        size_t i {0};
        double obs_ell {0}, obs_ell_old { - arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;

        // prepare for tail completion
        const arma::uvec case23_ind { vec_union(case2_ind_, case3_ind_) };
        double max_event_time { time(max_event_time_ind_) };
        if (tail_tau < 0)
            tail_tau = arma::datum::inf;

        // allow users to stop here
        Rcpp::checkUserInterrupt();

        // main loop of EM algorithm
        while (true) {

            // update to the latest estimates
            p_vec = cure_obj_.predict(cure_obj_.coef_);
            cox_obj_.compute_haz_surv_time();
            cox_obj_.compute_censor_haz_surv_time();

            // prepare for exponential tail completion method
            double s0_tau {0}, etail_lambda {0};
            if (tail_completion == 2) {
                s0_tau = cox_obj_.S0_time_(max_event_time_ind_);
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
                            cox_obj_.S_time_(j) = 0;
                            cox_obj_.S0_time_(j) = 0;
                        }
                        break;
                    case 1:
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            cox_obj_.S_time_(j) = 0;
                            cox_obj_.S0_time_(j) = 0;
                        }
                        break;
                    case 2:
                        // exponential tail by Peng (2003)
                        if (time(j) > max_event_time) {
                            cox_obj_.S0_time_(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            cox_obj_.S_time_(j) = std::pow(
                                cox_obj_.S0_time_(j),
                                std::exp(cox_obj_.xbeta_(j))
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
            for (size_t j: case1_ind_) {
                obs_ell += std::log(p_vec(j)) +
                    std::log(cox_obj_.h_time_(j)) -
                    cox_obj_.H_time_(j) - cox_obj_.Hc_time_(j);
            }
            // for case 2
            for (size_t j: case2_ind_) {
                obs_ell +=
                    std::log(p_vec(j) * cox_obj_.S_time_(j) + 1 - p_vec(j)) -
                    cox_obj_.Hc_time_(j) + std::log(cox_obj_.hc_time_(j));
            }
            // for case 3
            for (size_t j: case3_ind_) {
                double m12_common {
                    p_vec(j) * cox_obj_.S0_time_(j) * cox_obj_.Sc_time_(j)
                };
                double m1 { cox_obj_.h_time_(j) * m12_common };
                double m2 { cox_obj_.hc_time_(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
                };
                obs_ell += std::log(m1 + m2 + m3);
            }

            // if verbose
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(50, '=')
                            << "\niteration: " << i
                            << "\n  Cox coef: "
                            << arma2rvec(cox_obj_.coef_)
                            << "\n    relative diff: " << tol1
                            << "\n  cure coef: "
                            << arma2rvec(cure_obj_.coef_)
                            << "\n    relative diff: " << tol2
                            << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << "\n";
            }

            bool early_exit { false };
            // early exit if has any `nan`
            if (cox_obj_.coef_.has_nan() || cure_obj_.coef_.has_nan()) {
                obs_ell = - arma::datum::inf;
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function went to infinite."
                            << "\n";
                early_exit = true;
            }
            // early exit if the regularized objective function increased, which
            // is technically impossible and thus can serve as a warning
            if (obs_ell < obs_ell_old) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "The observed data log-likelihood decreased."
                                << "\n";
                }
                early_exit = early_exit || early_stop;
            }
            // return the estimates from last step
            if (early_exit) {
                if (verbose) {
                    Rcpp::Rcout << "Ended the EM algorithm after iteration "
                                << i
                                << " with estimates from last step."
                                << "\n";
                }
                // take the estimates from the last step
                cox_obj_.coef_ = cox_beta;
                cure_obj_.coef_ = cure_beta;
                // update hazard and survival function estimates
                cox_obj_.compute_haz_surv_time();
                cox_obj_.S0_time_ = s0_wi_tail;
                cox_obj_.S_time_ = s_wi_tail;
                cox_obj_.compute_censor_haz_surv_time();
                cox_obj_.est_haz_surv();
                // use old obs likelihood
                obs_ell = obs_ell_old;
                // break here
                break;
            }

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i >= em_max_iter) {
                if (verbose) {
                    if (i < em_max_iter) {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached convergence after " << i
                                    << " iterations\n" << "\n";
                    } else {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached the max iteration number."
                                    << "\n";
                    }
                }
                // compute hazard and survival function estimates
                cox_obj_.est_haz_surv();
                // get out of the loop here
                break;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // record estimates from last step
            cox_beta = cox_obj_.coef_;
            cure_beta = cure_obj_.coef_;
            obs_ell_old = obs_ell;
            s0_wi_tail = cox_obj_.S0_time_;
            s_wi_tail = cox_obj_.S_time_;

            // update iter for the next iteration
            ++i;

            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind_) {
                double numer_j { p_vec(j) * cox_obj_.S_time_(j)};
                estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
            }

            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind_) {
                double m12_common {
                    p_vec(j) * cox_obj_.S_time_(j) * cox_obj_.Sc_time_(j)
                };
                double m1 { cox_obj_.h_time_(j) * m12_common };
                double m2 { cox_obj_.hc_time_(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
                };
                double m { m1 + m2 + m3 };
                double w1 { m1 / m };
                double w2 { m2 / m };
                // some special care for subjects in case 3
                // since event(j) cannot be either 0 or 1!
                set_pmin_bound(w1, pmin);
                set_pmin_bound(w2, pmin);
                estep_m(j) = w1 + w2;
                event(j) = w1;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning the M-step for the survival layer:"
                            << "\n";
            }
            cox_obj_.set_offset_haz(arma::log(estep_m));
            cox_obj_.update_event_weight(event);
            cox_obj_.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol,
                         early_stop == 1, verbose > 2);
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nThe M-step for the survival layer was done."
                            << "\n";
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the Cure layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning M-step for the cure layer:";
            }
            cure_obj_.update_y(estep_m);
            cure_obj_.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol,
                          pmin, early_stop == 1, verbose > 2);
            if (verbose > 1) {
                Rcpp::Rcout << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done."
                            << "\n";
            }

            // update tolerance
            tol1 = rel_l2_norm(cox_obj_.coef_, cox_beta);
            tol2 = rel_l2_norm(cure_obj_.coef_, cure_beta);

        } // end of the EM algorithm

        // standardize
        if (cox_standardize_) {
            cox_obj_.standardize_ = true;
            cox_obj_.rescale_estimates();
            cox_obj_.est_haz_surv();
        }
        if (cure_standardize_) {
            cure_obj_.standardize_ = true;
            cure_obj_.rescale_coef();
        }

        // reset cox_obj_ and cure_obj_ in case of further usage
        cox_obj_.reset_offset_haz();
        // cox_obj_.set_offset(offset0);
        cure_obj_.update_y(cox_obj_.event_);

        // prepare outputs
        cox_coef_ = cox_obj_.coef_;
        cure_coef_ = cure_obj_.coef_;
        unique_time_ = cox_obj_.unique_time_;
        h0_est_ = cox_obj_.h0_est_;
        H0_est_ = cox_obj_.H0_est_;
        S0_est_ = cox_obj_.S0_est_;
        hc_est_ = cox_obj_.hc_est_;
        Hc_est_ = cox_obj_.Hc_est_;
        Sc_est_ = cox_obj_.Sc_est_;
        neg_ll_ = - obs_ell;
        coef_df_ = cox_obj_.coef_df_ + cure_obj_.coef_df_;
        n_iter_ = i;
        compute_bic1();
        compute_bic2();
        compute_aic();

        // prepare scores and prob in their original order
        const arma::uvec& rev_ord { cox_obj_.rev_ord_ };
        cox_xbeta_ = cox_obj_.xbeta_.elem(rev_ord);
        cure_xbeta_ = cure_obj_.xbeta_.elem(rev_ord);

        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj_.prob_vec_ };
        p_vec_event.elem(case1_ind_).ones();
        susceptible_prob_ = cure_obj_.prob_vec_.elem(rev_ord);

        // compute posterior probabilities from E-step
        for (size_t j: case2_ind_) {
            double numer_j { p_vec(j) *  cox_obj_.S_time_(j)};
            estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        for (size_t j: case3_ind_) {
            double m12_common {
                p_vec(j) * cox_obj_.S_time_(j) * cox_obj_.Sc_time_(j)
            };
            double m1 { cox_obj_.h_time_(j) * m12_common };
            double m2 { cox_obj_.hc_time_(j) * m12_common };
            double m3 {
                (1 - p_vec(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
            };
            double m { m1 + m2 + m3 };
            double w1 { m1 / m };
            double w2 { m2 / m };
            set_pmin_bound(w1, pmin);
            set_pmin_bound(w2, pmin);
            estep_m(j) = w1 + w2;
            event(j) = w1;
        }
        estep_cured_ = 1 - estep_m.elem(rev_ord);
        estep_event_ = event(rev_ord);
        estep_censor_ = 1 - estep_cured_ - estep_event_;

        // compute weighted c-index over the certain records
        c_index_ = Intsurv::Concordance(
            cox_obj_.time_.elem(cer_ind_),
            cox_obj_.event_.elem(cer_ind_),
            cox_obj_.xbeta_.elem(cer_ind_),
            p_vec_event.elem(cer_ind_)
            ).index_;
    }


    // fit regularized Cox cure model with uncertain events
    // with adaptive elastic net penalty for perticular lambda's
    // lambda_1 * lasso * factors + lambda_2 * ridge
    inline void CoxphCureMcar::regularized_fit(
        const double cox_l1_lambda = 0,
        const double cox_l2_lambda = 0,
        const double cure_l1_lambda = 0,
        const double cure_l2_lambda = 0,
        const arma::vec& cox_l1_penalty_factor = 0,
        const arma::vec& cure_l1_penalty_factor = 0,
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int em_max_iter = 500,
        const double em_rel_tol = 1e-5,
        const unsigned int cox_mstep_max_iter = 200,
        const double cox_mstep_rel_tol = 1e-5,
        const unsigned int cure_mstep_max_iter = 200,
        const double cure_mstep_rel_tol = 1e-5,
        const unsigned int tail_completion = 1,
        double tail_tau = -1,
        const double pmin = 1e-5,
        const unsigned int early_stop = 0,
        const unsigned int verbose = 0
        )
    {
        pmin_ = pmin;
        // get pre-processed design matrix, time, and event
        const arma::vec& time { cox_obj_.time_ };
        arma::vec event { cox_obj_.event_ };

        // L1 penalty factor for Cox model
        arma::vec cox_l1_penalty { arma::ones(cox_p_) };
        if (cox_l1_penalty_factor.n_elem == cox_p_) {
            // re-scale so that sum(factor) = number of predictors
            cox_l1_penalty = cox_l1_penalty_factor * cox_p_ /
                arma::sum(cox_l1_penalty_factor);
        }
        cox_l1_penalty_factor_ = cox_l1_penalty;
        // L1 penalty factor for cure model
        arma::vec cure_l1_penalty { arma::ones(cure_p0_) };
        if (cure_l1_penalty_factor.n_elem == cure_p0_) {
            // re-scale so that sum(factor) = number of predictors
            cure_l1_penalty = cure_l1_penalty_factor * cure_p0_ /
                arma::sum(cure_l1_penalty_factor);
        }
        cure_l1_penalty_factor_ = cure_l1_penalty;

        // initialized with all zeros coef
        arma::vec cox_beta { arma::zeros(cox_p_) };
        arma::vec cure_beta { arma::zeros(cure_p_) };

        // compute the large enough lambdas that result in all-zero estimates
        arma::vec cox_grad_zero { arma::abs(cox_obj_.gradient(cox_beta)) };
        arma::vec cure_grad_zero {
            arma::abs(cure_obj_.gradient(cure_beta))
        };
        cure_grad_zero = cure_grad_zero.tail(cure_l1_penalty.n_elem);
        // excluding variable with zero penalty factor
        arma::uvec cox_active_l1_penalty { arma::find(cox_l1_penalty > 0) };
        arma::uvec cure_active_l1_penalty { arma::find(cure_l1_penalty > 0) };
        cox_l1_lambda_max_ = arma::max(
            cox_grad_zero.elem(cox_active_l1_penalty) /
            cox_l1_penalty.elem(cox_active_l1_penalty)
            ) / n_obs_;
        cure_l1_lambda_max_ = arma::max(
            cure_grad_zero.elem(cure_active_l1_penalty) /
            cure_l1_penalty.elem(cure_active_l1_penalty)
            ) / n_obs_;

        // early stop if we want lambda_max by setting em_max_iter = 0
        if (em_max_iter == 0) {
            return;
        }

        // initialize cox_beta
        if (cox_start.n_elem == cox_p_) {
            cox_beta = cox_start;
        } else {
            cox_beta = arma::zeros(cox_p_);
        }
        arma::vec cox_exp_x_beta = arma::exp(cox_obj_.x_ * cox_beta);
        // initialize cure_beta
        if (cure_start.n_elem == cure_p_) {
            cure_beta = cure_start;
        } else {
            cure_beta = arma::zeros(cure_p_);
        }

        // initialize coef estimates
        cure_obj_.coef_ = cure_beta;
        cox_obj_.coef_ = cox_beta;
        cure_obj_.coef_df_ = compute_coef_df(cure_beta);
        cox_obj_.coef_df_ = compute_coef_df(cox_beta);

        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(cer_ind_),
                        event.elem(cer_ind_))
        };
        cox_obj_.h0_time_ = nelen_event.step_inst_rate(time);
        cox_obj_.H0_time_ = nelen_event.step_cum_rate(time);
        // for censoring
        NelsonAalen nelen_censor {
            NelsonAalen(time.elem(cer_ind_),
                        1 - event.elem(cer_ind_))
        };
        cox_obj_.hc_time_ = nelen_censor.step_inst_rate(time);
        cox_obj_.Hc_time_ = nelen_censor.step_cum_rate(time);
        // update related function values
        cox_obj_.h_time_ = cox_obj_.h0_time_ % cox_exp_x_beta;
        cox_obj_.H_time_ = cox_obj_.H0_time_ % cox_exp_x_beta;
        cox_obj_.S0_time_ = arma::exp(- cox_obj_.H0_time_);
        cox_obj_.S_time_ = arma::exp(- cox_obj_.H_time_);
        cox_obj_.Sc_time_ = arma::exp(- cox_obj_.Hc_time_);

        // intialization for the main loop
        arma::vec p_vec { arma::zeros(n_obs_) };
        arma::vec estep_m { event };
        size_t i {0};
        double obs_ell {0};
        double reg_obj {0};
        double reg_obj_old { arma::datum::inf }, obs_ell_old { reg_obj_old };
        double bic1_old { arma::datum::inf }, bic2_old { bic1_old };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;

        // prepare for tail completion
        const arma::uvec case23_ind { vec_union(case2_ind_, case3_ind_) };
        double max_event_time { time(max_event_time_ind_) };
        if (tail_tau < 0)
            tail_tau = arma::datum::inf;

        // allow users to stop here
        Rcpp::checkUserInterrupt();

        // main loop of EM algorithm
        while (true) {
            // update to the latest estimates
            p_vec = cure_obj_.predict(cure_obj_.coef_);
            cox_obj_.compute_haz_surv_time();
            cox_obj_.compute_censor_haz_surv_time();

            // prepare for exponential tail completion method
            double s0_tau {0}, etail_lambda {0};
            if (tail_completion == 2) {
                s0_tau = cox_obj_.S0_time_(max_event_time_ind_);
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
                            cox_obj_.S_time_(j) = 0;
                            cox_obj_.S0_time_(j) = 0;
                        }
                        break;
                    case 1:
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            cox_obj_.S_time_(j) = 0;
                            cox_obj_.S0_time_(j) = 0;
                        }
                        break;
                    case 2:
                        // exponential tail by Peng (2003)
                        if (time(j) > max_event_time) {
                            cox_obj_.S0_time_(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            cox_obj_.S_time_(j) = std::pow(
                                cox_obj_.S0_time_(j),
                                std::exp(cox_obj_.xbeta_(j))
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
            for (size_t j: case1_ind_) {
                obs_ell += std::log(p_vec(j)) +
                    std::log(cox_obj_.h_time_(j)) +
                    std::log(cox_obj_.S_time_(j)) +
                    std::log(cox_obj_.Sc_time_(j));
            }
            // for case 2
            for (size_t j: case2_ind_) {
                obs_ell +=
                    std::log(p_vec(j) * cox_obj_.S_time_(j) + 1 - p_vec(j)) +
                    std::log(cox_obj_.Sc_time_(j)) +
                    std::log(cox_obj_.hc_time_(j));
            }
            // for case 3
            for (size_t j: case3_ind_) {
                double m1 {
                    p_vec(j) * cox_obj_.h_time_(j) *
                    cox_obj_.S_time_(j) * cox_obj_.Sc_time_(j)
                };
                double m2 {
                    p_vec(j) * cox_obj_.hc_time_(j) *
                    cox_obj_.Sc_time_(j) * cox_obj_.S_time_(j)
                };
                double m3 {
                    (1 - p_vec(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
                };
                obs_ell += std::log(m1 + m2 + m3);
            }
            // compuete the regularized objective function
            double reg_cox {
                cox_l1_lambda_ * l1_norm(cox_obj_.coef_ % cox_l1_penalty) +
                cox_l2_lambda_ * sum_of_square(cox_obj_.coef_)
            };
            double reg_cure {
                cure_l1_lambda_ *
                l1_norm(cure_obj_.coef_.tail(cure_p0_) % cure_l1_penalty) +
                cure_l2_lambda_ *
                sum_of_square(cure_obj_.coef_.tail(cure_p0_))
            };
            reg_obj = - obs_ell / n_obs_ + reg_cox + reg_cure;

            // compute bic
            neg_ll_ = - obs_ell;
            coef_df_ = cox_obj_.coef_df_ + cure_obj_.coef_df_;
            compute_bic1();
            compute_bic2();
            compute_aic();

            // if verbose
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(50, '=')
                            << "\niteration: " << i
                            << "\n  Cox coef: "
                            << arma2rvec(cox_obj_.coef_)
                            << "\n    relative diff: " << tol1
                            << "\n  cure coef: "
                            << arma2rvec(cure_obj_.coef_)
                            << "\n    relative diff: " << tol2
                            << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << "\n  regularized objective function: "
                            << reg_obj
                            << "\n    penalty on Cox model: "
                            << reg_cox
                            << "\n    penalty on cure layer: "
                            << reg_cure
                            << "\n";
            }

            bool early_exit { false };
            // early exit if has any `nan`
            if (cox_obj_.coef_.has_nan() || cure_obj_.coef_.has_nan()) {
                obs_ell = - arma::datum::inf;
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function went to infinite."
                            << "\n";
                early_exit = true;
            }
            // early exit if the regularized objective function increased, which
            // is technically impossible and thus can serve as a warning
            if (reg_obj > reg_obj_old) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "The observed data log-likelihood decreased."
                                << "\n";
                }
                early_exit = early_exit || early_stop == 1;
            }
            // early exit if bic increased
            if (bic1_ > bic1_old && bic2_ > bic2_old) {
                if (verbose) {
                    Rcpp::Rcout << "The BIC increased."
                                << "\n";
                }
                early_exit = early_exit || early_stop == 2;
            }
            // return the estimates from last step
            if (early_exit) {
                if (verbose) {
                    Rcpp::Rcout << "Ended the EM algorithm after iteration "
                                << i
                                << " with estimates from last step."
                                << "\n";
                }
                // take the estimates from the last step
                cox_obj_.coef_ = cox_beta;
                cure_obj_.coef_ = cure_beta;
                // update coef_df_ and en_coef
                // cox_obj_.set_en_coef(cox_l2_lambda_);
                // cure_obj_.set_en_coef(cure_l2_lambda_);

                // update hazard and survival function estimates
                cox_obj_.compute_haz_surv_time();
                cox_obj_.S0_time_ = s0_wi_tail;
                cox_obj_.S_time_ = s_wi_tail;
                cox_obj_.compute_censor_haz_surv_time();
                cox_obj_.est_haz_surv();
                // use old obs likelihood
                obs_ell = obs_ell_old;
                // break here
                break;
            }

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i >= em_max_iter) {
                if (verbose) {
                    if (i < em_max_iter) {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached convergence after " << i
                                    << " iterations\n" << "\n";
                    } else {
                        Rcpp::Rcout << "\n" << std::string(50, '=') << "\n"
                                    << "reached the max iteration number."
                                    << "\n";
                    }
                }
                // compute hazard and survival function estimates
                cox_obj_.est_haz_surv();
                // get out of the loop here
                break;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // record estimates from last step
            cox_beta = cox_obj_.coef_;
            cure_beta = cure_obj_.coef_;
            obs_ell_old = obs_ell;
            reg_obj_old = reg_obj;
            s0_wi_tail = cox_obj_.S0_time_;
            s_wi_tail = cox_obj_.S_time_;

            // update iter for the next iteration
            ++i;

            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind_) {
                double numer_j { p_vec(j) * cox_obj_.S_time_(j)};
                estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
            }

            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind_) {
                double m12_common {
                    p_vec(j) * cox_obj_.S0_time_(j) * cox_obj_.Sc_time_(j)
                };
                double m1 { cox_obj_.h_time_(j) * m12_common };
                double m2 { cox_obj_.hc_time_(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
                };
                double m { m1 + m2 + m3 };
                double w1 { m1 / m };
                double w2 { m2 / m };
                // some special care for subjects in case 3
                // since event(j) cannot be either 0 or 1!
                set_pmin_bound(w1, pmin);
                set_pmin_bound(w2, pmin);
                estep_m(j) = w1 + w2;
                event(j) = w1;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning the M-step for the survival layer:"
                            << "\n";
            }
            cox_obj_.set_offset_haz(arma::log(estep_m));
            cox_obj_.update_event_weight(event);
            cox_obj_.regularized_fit(
                cox_l1_lambda_, cox_l2_lambda_, cox_l1_penalty_factor_,
                cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol,
                early_stop == 1, verbose > 2
                );
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nThe M-step for the survival layer was done."
                            << "\n";
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the Cure layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning M-step for the cure layer:";
            }
            cure_obj_.update_y(estep_m);
            cure_obj_.regularized_fit(
                cure_l1_lambda_, cure_l2_lambda_, cure_l1_penalty_factor_,
                cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol,
                pmin, early_stop == 1, verbose > 2
                );
            if (verbose > 1) {
                Rcpp::Rcout << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done."
                            << "\n";
            }

            // update tolerance
            tol1 = rel_l2_norm(cox_obj_.coef_, cox_beta);
            tol2 = rel_l2_norm(cure_obj_.coef_, cure_beta);

        } // end of the EM algorithm

        // standardize
        if (cox_standardize_) {
            cox_obj_.standardize_ = cox_standardize_;
            cox_obj_.rescale_estimates();
            cox_obj_.est_haz_surv();
        }
        if (cure_standardize_) {
            cure_obj_.standardize_ = cure_standardize_;
            cure_obj_.rescale_coef();
        }

        // reset cox_obj_ and cure_obj_ in case of further usage
        cox_obj_.reset_offset_haz();
        // cox_obj_.set_offset(offset0);
        cure_obj_.update_y(cox_obj_.event_);

        // prepare outputs
        cox_coef_ = cox_obj_.coef_;
        cure_coef_ = cure_obj_.coef_;
        // cox_en_coef_ = cox_obj_.en_coef_;
        // cure_en_coef_ = cure_obj_.en_coef_;

        unique_time_ = cox_obj_.unique_time_;
        h0_est_ = cox_obj_.h0_est_;
        H0_est_ = cox_obj_.H0_est_;
        S0_est_ = cox_obj_.S0_est_;
        hc_est_ = cox_obj_.hc_est_;
        Hc_est_ = cox_obj_.Hc_est_;
        Sc_est_ = cox_obj_.Sc_est_;

        neg_ll_ = - obs_ell;
        coef_df_ = cox_obj_.coef_df_ + cure_obj_.coef_df_;
        cox_l1_lambda_ = cox_l1_lambda;
        cox_l2_lambda_ = cox_l2_lambda;
        cure_l1_lambda_ = cure_l1_lambda;
        cure_l2_lambda_ = cure_l2_lambda;
        n_iter_ = i;
        compute_bic1();
        compute_bic2();
        compute_aic();

        // prepare scores and prob in their original order
        const arma::uvec& rev_ord { cox_obj_.ord_ };
        cox_xbeta_ = cox_obj_.xbeta_.elem(rev_ord);
        cure_xbeta_ = cure_obj_.xbeta_.elem(rev_ord);

        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj_.prob_vec_ };
        p_vec_event.elem(case1_ind_).ones();
        susceptible_prob_ = cure_obj_.prob_vec_.elem(rev_ord);

        // compute posterior probabilities from E-step
        for (size_t j: case2_ind_) {
            double numer_j { p_vec(j) *  cox_obj_.S_time_(j)};
            estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        for (size_t j: case3_ind_) {
            double m12_common {
                p_vec(j) * cox_obj_.S0_time_(j) * cox_obj_.Sc_time_(j)
            };
            double m1 { cox_obj_.h_time_(j) * m12_common };
            double m2 { cox_obj_.hc_time_(j) * m12_common };
            double m3 {
                (1 - p_vec(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
            };
            double m { m1 + m2 + m3 };
            double w1 { m1 / m };
            double w2 { m2 / m };
            estep_m(j) = w1 + w2;
            event(j) = w1;
        }
        estep_cured_ = 1 - estep_m.elem(rev_ord);
        estep_event_ = event(rev_ord);
        estep_censor_ = 1 - estep_cured_ - estep_event_;

        // compute weighted c-index over the certain records
        c_index_ = Intsurv::Concordance(
            cox_obj_.time_.elem(cer_ind_),
            cox_obj_.event_.elem(cer_ind_),
            cox_obj_.xbeta_.elem(cer_ind_),
            p_vec_event.elem(cer_ind_)
            ).index_;
    }

    // function to compute the observe data log-likelihood function
    // for given fitted model and estimates
    inline double CoxphCureMcar::obs_log_likelihood() const
    {
        double obs_ell { 0 };
        arma::vec sus_prob { cure_obj_.prob_vec_ };
        // for case 1
        for (size_t j: case1_ind_) {
            obs_ell += std::log(sus_prob(j)) +
                std::log(cox_obj_.h_time_(j)) -
                cox_obj_.H_time_(j) - cox_obj_.Hc_time_(j);
        }
        // for case 2
        for (size_t j: case2_ind_) {
            obs_ell +=
                std::log(sus_prob(j) * cox_obj_.S_time_(j) + 1 - sus_prob(j)) -
                cox_obj_.Hc_time_(j) + std::log(cox_obj_.hc_time_(j));
        }
        // for case 3
        for (size_t j: case3_ind_) {
            double m12_common {
                sus_prob(j) * cox_obj_.S0_time_(j) * cox_obj_.Sc_time_(j)
            };
            double m1 { cox_obj_.h_time_(j) * m12_common };
            double m2 { cox_obj_.hc_time_(j) * m12_common };
            double m3 {
                (1 - sus_prob(j)) * cox_obj_.hc_time_(j) * cox_obj_.Sc_time_(j)
            };
            obs_ell += std::log(m1 + m2 + m3);
        }
        return obs_ell;
    }

    // for given fitted model and a new set of data
    inline double CoxphCureMcar::obs_log_likelihood(
        // all the inputs will be sorted inside of the function
        // so their copies are asked here
        arma::vec new_time,
        arma::vec new_event,
        arma::mat new_cox_x,
        arma::mat new_cure_x,
        arma::vec new_cox_offset = 0,
        arma::vec new_cure_offset = 0
        ) const
    {
        // check if the number of covariates matchs the fitted model
        if (new_cox_x.n_cols != cox_p_) {
            throw std::range_error(
                "The number of columns ('new_cox_x') must match the model."
                );
        }
        if (new_cure_x.n_cols != cure_p0_) {
            throw std::range_error(
                "The number of columns ('new_cure_x') must match the model."
                );
        }
        // number of observations
        unsigned int new_n_obs { new_cox_x.n_rows };
        if (new_cure_x.n_rows != new_n_obs ||
            new_time.n_elem != new_n_obs ||
            new_event.n_elem != new_n_obs) {
            throw std::range_error(
                "The number of rows of the new data must be the same."
                );
        }
        const double const4na { 0.5 };
        new_event.replace(arma::datum::nan, const4na);
        double obs_ell { 0 };
        // sort based on time and event
        // time: ascending order
        // event: events first, then censoring at the same time point
        arma::uvec des_event_ind { arma::sort_index(new_event, "descend") };
        arma::uvec asc_time_ind {
            arma::stable_sort_index(new_time.elem(des_event_ind), "ascend")
        };
        arma::uvec ord { des_event_ind.elem(asc_time_ind) };
        // do actual sorting
        new_time = new_time.elem(ord);
        new_event = new_event.elem(ord);
        new_cox_x = new_cox_x.rows(ord);
        new_cure_x = new_cure_x.rows(ord);
        // process offset terms
        if (new_cox_offset.n_elem != new_cox_x.n_rows) {
            new_cox_offset = arma::zeros(new_cox_x.n_rows);
        } else {
            new_cox_offset = new_cox_offset(ord);
        }
        if (new_cure_offset.n_elem != new_cure_x.n_rows) {
            new_cure_offset = arma::zeros(new_cure_x.n_rows);
        } else {
            new_cure_offset = new_cure_offset(ord);
        }

        // add intercept if needed
        if (cure_p_ > cure_p0_) {
            new_cure_x = arma::join_horiz(
                arma::ones(new_cure_x.n_rows), new_cure_x
                );
        }
        arma::uvec new_case1_ind { arma::find(new_event > const4na) };
        arma::uvec new_case2_ind { arma::find(new_event < const4na) };
        arma::uvec new_case3_ind { arma::find(new_event == const4na) };

        // construct the baseline survival curve
        // tail completion has already been applied to S0_est_
        arma::vec S0_vec {
            arma::join_cols(arma::ones<arma::vec>(1), S0_est_)
        };
        arma::vec Sc0_vec {
            arma::join_cols(arma::ones<arma::vec>(1), Sc_est_)
        };
        // baseline estimates
        arma::vec S_vec {
            step_fun(new_time, unique_time_, S0_vec)
        };
        arma::vec Sc_vec {
            step_fun(new_time, unique_time_, Sc0_vec)
        };
        arma::vec H_vec { - arma::log(S_vec) };
        arma::vec Hc_vec { - arma::log(Sc_vec) };
        // only consider positive values
        arma::uvec which_h { arma::find(h0_est_ > 0) };
        arma::uvec which_hc { arma::find(hc_est_ > 0) };
        arma::vec h_vec {
            step_fun2(new_time,
                      unique_time_.elem(which_h),
                      h0_est_.elem(which_h))
        };
        arma::vec hc_vec {
            step_fun2(new_time,
                      unique_time_.elem(which_hc),
                      hc_est_.elem(which_hc))
        };
        // apply x * beta
        // compute parts for the new data
        arma::vec new_cox_xbeta {
            mat2vec(new_cox_x * cox_coef_) + new_cox_offset
        };
        arma::vec exp_cox_xbeta { arma::exp(new_cox_xbeta) };
        h_vec %= exp_cox_xbeta;
        H_vec %= exp_cox_xbeta;
        S_vec = arma::exp(- H_vec);
        arma::vec new_cure_xgamma {
            mat2vec(new_cure_x * cure_coef_) + new_cure_offset
        };
        arma::vec p_vec { 1 / (1 + arma::exp(- new_cure_xgamma)) };
        set_pmin_bound(p_vec, cure_obj_.pmin_);
        // for case 1
        for (size_t j: new_case1_ind) {
            obs_ell += std::log(p_vec(j)) +
                std::log(h_vec(j)) -
                cox_obj_.H_time_(j) - cox_obj_.Hc_time_(j);
        }
        // for case 2
        for (size_t j: new_case2_ind) {
            obs_ell +=
                std::log(p_vec(j) * S_vec(j) + 1 - p_vec(j)) -
                Hc_vec(j) + std::log(hc_vec(j));
        }
        // for case 3
        for (size_t j: new_case3_ind) {
            double m1 {
                p_vec(j) * h_vec(j) *
                S_vec(j) * Sc_vec(j)
            };
            double m2 {
                p_vec(j) * hc_vec(j) *
                Sc_vec(j) * S_vec(j)
            };
            double m3 {
                (1 - p_vec(j)) * hc_vec(j) * Sc_vec(j)
            };
            obs_ell += std::log(m1 + m2 + m3);
        }
        return obs_ell;
    }

}


#endif
