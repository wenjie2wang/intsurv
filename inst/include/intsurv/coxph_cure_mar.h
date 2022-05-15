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

#ifndef INTSURV_COXPH_CURE_MAR_H
#define INTSURV_COXPH_CURE_MAR_H

#include <RcppArmadillo.h>
#include "assessment.h"
#include "coxph_reg.h"
#include "logistic_reg.h"
#include "nonparametric.h"
#include "control.h"
#include "utils.h"

namespace Intsurv {

    class CoxphCureMar
    {
    protected:
        // caches
        unsigned int cure_p0_;
        double dn_obs_;          // double version of n_obs_
        bool surv_standardize0_; // original standardization option
        bool cure_standardize0_; // original standardization option
        arma::vec event0_;       // sorted original input event for recovery

    public:
        CoxphReg surv_obj_;
        LogisticReg cure_obj_;
        LogisticReg mar_obj_;

        arma::uvec case1_ind_;
        arma::uvec case2_ind_;
        arma::uvec cer_ind_;    // index of rows with correct event indicators
        arma::uvec case3_ind_;

        unsigned int surv_p_;
        unsigned int cure_p_;
        unsigned int max_event_time_ind_;

        // outputs
        arma::vec surv_coef_;
        arma::vec cure_coef_;
        arma::vec mar_coef_;
        arma::vec surv_xbeta_;  // score from the survival layer
        arma::vec cure_xbeta_;  // score from the cure layer
        arma::vec mar_xbeta_;  // score from the cure layer

        // for each subject and in the original order of X
        arma::vec susceptible_prob_; // probability of being susceptible
        arma::vec mar_prob_;         // probability of case 3
        // values in the last E-step
        arma::vec estep_cured_;
        arma::vec estep_event_;
        arma::vec estep_censor_;

        unsigned int surv_coef_df_;
        unsigned int cure_coef_df_;
        unsigned int mar_coef_df_;
        unsigned int coef_df_;

        double neg_ll_;         // negative log-likelihood
        unsigned int n_obs_;    // number of observations
        unsigned int n_event_;  // number of certain events
        // BIC: log(num_obs) * coef_df_ + 2 * neg_ll_
        double bic1_;
        // BIC: log(num_certain_event) * coef_df_ + 2 * neg_ll_
        double bic2_;
        double aic_;            // AIC: 2 * coef_df_ + 2 * neg_ll_
        double c_index_;        // weighted C-index

        // hazard and survival function estimates at unique time
        arma::vec unique_time_;
        arma::vec h0_est_;
        arma::vec H0_est_;
        arma::vec S0_est_;
        arma::vec hc_est_;
        arma::vec Hc_est_;
        arma::vec Sc_est_;

        // controls
        Control control_;
        unsigned int n_iter_;   // number of iterations

        // default constructor
        CoxphCureMar() {}

        // constructors
        CoxphCureMar(
            const arma::vec& time,
            const arma::vec& event,
            const arma::mat& surv_x,
            const arma::mat& cure_x,
            const arma::mat& mar_x,
            const Control& control,
            const Control& surv_control,
            const Control& cure_control,
            const Control& mar_control) :
            control_ (control)
        {
            // replace NA or NaN event indicator with 0.5
            // (or any number between 0 and 1)
            arma::vec event0na { event };
            event0_ = event;
            const double const4na { 0.5 };
            event0na.replace(arma::datum::nan, const4na);
            // create the CoxphReg object
            surv_obj_ = CoxphReg(time, event0na, surv_x, surv_control);
            // pre-process x and y
            surv_p_ = surv_x.n_cols;
            cure_p0_ = cure_x.n_cols;
            cure_p_ = cure_p0_ +
                static_cast<unsigned int>(cure_control.intercept_);
            n_obs_ = surv_x.n_rows;
            dn_obs_ = static_cast<double>(n_obs_);
            const arma::uvec& surv_sort_ind { surv_obj_.ord_ };
            event0_ = event0_.elem(surv_sort_ind);
            const arma::vec& s_event { surv_obj_.event_ };
            case1_ind_ = arma::find(s_event > const4na);
            case2_ind_ = arma::find(s_event < const4na);
            n_event_ = case1_ind_.n_elem;
            cer_ind_ = vec_union(case1_ind_, case2_ind_);
            case3_ind_ = arma::find(s_event == const4na);
            max_event_time_ind_ = arma::max(case1_ind_);
            // create a cure object
            cure_obj_ = LogisticReg(cure_x.rows(surv_sort_ind),
                                    s_event,
                                    cure_control);
            // create a mar object
            arma::vec mar_y { arma::zeros(n_obs_) };
            mar_y.elem(case3_ind_).ones();
            mar_obj_ = LogisticReg(mar_x.rows(surv_sort_ind),
                                   mar_y,
                                   mar_control);
            // avoid standardization after each iteration
            surv_standardize0_ = surv_obj_.control_.standardize_;
            surv_obj_.control_.set_standardize(false);
            cure_standardize0_ = cure_obj_.control_.standardize_;
            cure_obj_.control_.set_standardize(false);
        }

        // function members
        inline arma::vec get_event(const bool resort = false) const
        {
            if (resort) {
                return event0_.elem(surv_obj_.rev_ord_);
            }
            return event0_;
        }

        // fit mar part
        inline void mar_fit()
        {
            mar_obj_.fit();
            mar_coef_ = mar_obj_.coef_;
            mar_xbeta_ = mar_obj_.get_xbeta();
            mar_prob_ = mar_obj_.predict();
        }
        inline void mar_net_fit()
        {
            mar_obj_.net_fit();
            mar_coef_ = mar_obj_.coef_;
            mar_xbeta_ = mar_obj_.get_xbeta();
            mar_prob_ = mar_obj_.predict();
        }
        inline void mar_net_path()
        {
            mar_obj_.net_path();
            mar_coef_ = mar_obj_.coef_;
            mar_xbeta_ = mar_obj_.get_xbeta();
            mar_prob_ = mar_obj_.predict();
        }

        // fit the Cox cure model with uncertain events by EM algorithm
        inline void fit();

        // fit regularized Cox cure model with adaptive elastic net penalty
        // for perticular lambda's
        inline void net_fit();

        // function to compute the observe data log-likelihood function
        // for given fitted model and estimates
        inline double obs_log_likelihood() const;

        // for given fitted model and a new set of data
        // all the inputs will be sorted inside of the function
        // so their copies are asked here
        inline double obs_log_likelihood(
            arma::vec new_time,
            arma::vec new_event,
            arma::mat new_surv_x,
            arma::mat new_cure_x,
            arma::vec new_surv_offset,
            arma::vec new_cure_offset
            ) const;

        inline double obs_log_likelihood(const CoxphCureMar& object) const;

        // compute BIC
        inline void compute_bic1() {
            bic1_ = std::log(dn_obs_) * coef_df_ + 2 * neg_ll_;
        }
        inline void compute_bic2() {
            bic2_ = std::log(static_cast<double>(case1_ind_.n_elem)) *
                coef_df_ + 2 * neg_ll_;
        }
        inline void compute_aic() {
            aic_ = 2 * (coef_df_ + neg_ll_);
        }

    };                          // end of class definition


    // fit the Cox cure mode by EM algorithm
    inline void CoxphCureMar::fit()
    {
        // initialize surv_beta
        const arma::vec& time { surv_obj_.time_ };
        const arma::vec& event { surv_obj_.event_ };
        arma::vec surv_beta { arma::zeros(surv_p_) };
        if (surv_obj_.control_.start_.n_elem == surv_p_) {
            surv_obj_.set_start();
            surv_beta = surv_obj_.control_.start_;
        } else {
            Control tmp_control { surv_obj_.control_ };
            tmp_control.set_standardize(false)->set_verbose(0);
            CoxphReg tmp_object {
                time.elem(case1_ind_),
                event.elem(case1_ind_),
                surv_obj_.x_.rows(case1_ind_),
                tmp_control
            };
            tmp_object.fit();
            surv_beta = tmp_object.coef_;
        }
        surv_obj_.control_.start_ = surv_beta;
        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(cure_p_) };
        if (cure_obj_.control_.start_.n_elem == cure_p_) {
            cure_obj_.set_start();
            cure_beta = cure_obj_.control_.start_;
        } else {
            Control tmp_control { cure_obj_.control_ };
            cure_obj_.control_.set_standardize(false)->
                set_verbose(0);
            cure_obj_.fit();
            cure_beta = cure_obj_.coef_;
            cure_obj_.control_ = tmp_control;
        }
        cure_obj_.control_.start_ = cure_beta;
        surv_obj_.coef0_ = surv_obj_.coef_ = surv_beta;
        cure_obj_.coef0_ = cure_obj_.coef_ = cure_beta;
        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(cer_ind_),
                        event.elem(cer_ind_))
        };
        surv_obj_.h0_time_ = nelen_event.step_inst_rate(time);
        surv_obj_.H0_time_ = nelen_event.step_cum_rate(time);
        // for censoring
        NelsonAalen nelen_censor {
            NelsonAalen(time.elem(cer_ind_),
                        1 - event.elem(cer_ind_))
        };
        surv_obj_.hc_time_ = nelen_censor.step_inst_rate(time);
        surv_obj_.Hc_time_ = nelen_censor.step_cum_rate(time);
        // update related function values
        arma::vec surv_exp_x_beta {
            arma::exp(surv_obj_.x_ * surv_beta)
        };
        surv_obj_.h_time_ = surv_obj_.h0_time_ % surv_exp_x_beta;
        surv_obj_.H_time_ = surv_obj_.H0_time_ % surv_exp_x_beta;
        surv_obj_.S0_time_ = arma::exp(- surv_obj_.H0_time_);
        surv_obj_.S_time_ = arma::exp(- surv_obj_.H_time_);
        surv_obj_.Sc_time_ = arma::exp(- surv_obj_.Hc_time_);
        // intialization for the main loop
        arma::vec p_vec { arma::zeros(n_obs_) };
        arma::vec estep_m { event }, estep_event { event };
        double obs_ell {0}, obs_ell_old { - arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;
        // prepare for tail completion
        const arma::uvec case23_ind { vec_union(case2_ind_, case3_ind_) };
        double max_event_time { time(max_event_time_ind_) };
        if (control_.tail_tau_ < 0)
            control_.tail_tau_ = arma::datum::inf;
        // for exp tail completion only
        double s0_tau, etail_lambda;
        // allow users to stop here
        Rcpp::checkUserInterrupt();
        // main loop of EM algorithm
        for (size_t i {0}; i < control_.max_iter_; ++i) {
            n_iter_ = i + 1;
            // update to the latest estimates
            p_vec = cure_obj_.predict();
            surv_obj_.compute_haz_surv_time();
            surv_obj_.compute_censor_haz_surv_time();
            // record estimates from last step
            surv_beta = surv_obj_.coef_;
            cure_beta = cure_obj_.coef_;
            s0_wi_tail = surv_obj_.S0_time_;
            s_wi_tail = surv_obj_.S_time_;
            // set them as starting values
            surv_obj_.control_.set_start(surv_beta);
            cure_obj_.control_.set_start(cure_beta);
            // tail completion for the conditional survival function
            switch(control_.tail_completion_) {
                case 0:
                    for (size_t j: case2_ind_) {
                        // tail completion after the given tail_tau
                        // by default, it means no tail completion
                        if (time(j) > control_.tail_tau_) {
                            surv_obj_.S_time_(j) = 0.0;
                            surv_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 1:
                    for (size_t j: case2_ind_) {
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            surv_obj_.S_time_(j) = 0.0;
                            surv_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 2:
                    s0_tau = surv_obj_.S0_time_(max_event_time_ind_);
                    etail_lambda = - std::log(s0_tau / max_event_time);
                    for (size_t j: case2_ind_) {
                        // exponential tail by Peng (2003)
                        arma::vec surv_xbeta { surv_obj_.get_xbeta() };
                        if (time(j) > max_event_time) {
                            surv_obj_.S0_time_(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            surv_obj_.S_time_(j) = std::pow(
                                surv_obj_.S0_time_(j),
                                std::exp(surv_xbeta(j))
                                );
                        }
                    }
                    break;
                default:    // do nothing, otherwise
                    break;
            }
            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind_) {
                double numer_j { p_vec(j) * surv_obj_.S_time_(j)};
                estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
            }
            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind_) {
                double m12_common {
                    p_vec(j) * surv_obj_.S_time_(j) * surv_obj_.Sc_time_(j)
                };
                double m1 { surv_obj_.h_time_(j) * m12_common };
                double m2 { surv_obj_.hc_time_(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * surv_obj_.hc_time_(j) *
                    surv_obj_.Sc_time_(j)
                };
                double m { m1 + m2 + m3 };
                double w1 { m1 / m };
                double w2 { m2 / m };
                // some special care for subjects in case 3
                // since event(j) cannot be either 0 or 1!
                set_pmin_bound(w1, cure_obj_.control_.pmin_);
                set_pmin_bound(w2, cure_obj_.control_.pmin_);
                estep_m(j) = w1 + w2;
                estep_event(j) = w1;
            }
            Rcpp::checkUserInterrupt();
            // M-step for the survival layer
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "Running the M-step for the survival layer:"
                            << "\n";
            }
            surv_obj_.set_offset_haz(arma::log(estep_m));
            surv_obj_.update_event_weight(estep_event);
            surv_obj_.fit();
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "The M-step for the survival layer was done."
                            << "\n";
            }
            Rcpp::checkUserInterrupt();
            // M-step for the Cure layer
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "Running M-step for the cure layer:"
                            << "\n";
            }
            cure_obj_.update_y(estep_m);
            cure_obj_.fit();
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "The M-step for the cure layer was done."
                            << "\n";
            }
            // check if has any `nan`
            if (surv_obj_.coef_.has_nan() || cure_obj_.coef_.has_nan()) {
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function might go to infinite.\n"
                            << "Try to return etimates from last step.\n";
                // take the estimates from the last step
                surv_obj_.coef_ = surv_beta;
                cure_obj_.coef_ = cure_beta;
                surv_obj_.compute_haz_surv_time();
                surv_obj_.S0_time_ = s0_wi_tail;
                surv_obj_.S_time_ = s_wi_tail;
                // surv_obj_.compute_censor_haz_surv_time();
                surv_obj_.est_haz_surv();
                break;
            }
            // update tolerance
            tol1 = rel_l1_norm(surv_obj_.coef_, surv_beta);
            tol2 = rel_l1_norm(cure_obj_.coef_, cure_beta);
            if (control_.verbose_ > 0) {
                // compute observed data log-likelihood
                obs_ell = 0;
                // for case 1
                for (size_t j: case1_ind_) {
                    obs_ell += std::log(mar_prob_(j)) +
                        std::log(p_vec(j)) +
                        std::log(surv_obj_.h_time_(j)) -
                        surv_obj_.H_time_(j) - surv_obj_.Hc_time_(j);
                }
                // for case 2
                for (size_t j: case2_ind_) {
                    obs_ell += std::log(mar_prob_(j)) +
                        std::log(p_vec(j) * surv_obj_.S_time_(j) +
                                 1 - p_vec(j)) -
                        surv_obj_.Hc_time_(j) +
                        std::log(surv_obj_.hc_time_(j));
                }
                // for case 3
                for (size_t j: case3_ind_) {
                    double m12_common {
                        p_vec(j) * surv_obj_.S0_time_(j) * surv_obj_.Sc_time_(j)
                    };
                    double m1 { surv_obj_.h_time_(j) * m12_common };
                    double m2 { surv_obj_.hc_time_(j) * m12_common };
                    double m3 {
                        (1 - p_vec(j)) * surv_obj_.hc_time_(j) *
                        surv_obj_.Sc_time_(j)
                    };
                    obs_ell += std::log(1 - mar_prob_(j)) +
                        std::log(m1 + m2 + m3);
                }
                Rcpp::Rcout << "\n"
                            << std::string(50, '=')
                            << "\niteration: "
                            << n_iter_;
                if (control_.verbose_ > 2) {
                    Rcpp::Rcout << "\n  surv coef: "
                                << arma2rvec(surv_obj_.coef_)
                                << "\n    relative diff: "
                                << tol1
                                << "\n  cure coef: "
                                << arma2rvec(cure_obj_.coef_)
                                << "\n    relative diff: "
                                << tol2;
                }
                Rcpp::Rcout << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << "\n";
                // check if the observed data log-likelihood decreased, which
                // is technically impossible and thus can serve as a warning
                // likely to be a coding bug
                if (obs_ell < obs_ell_old) {
                    Rcpp::Rcout << "Warning: the observed data "
                                << "log-likelihood somehow decreased.\n";
                }
                obs_ell_old = obs_ell;
            }
            // check convergence
            if (tol1 < control_.epsilon_ && tol2 < control_.epsilon_) {
                if (control_.verbose_ > 0) {
                    Rcpp::Rcout << "\n"
                                << std::string(50, '=')
                                << "\n"
                                << "reached convergence after "
                                << n_iter_
                                << " iterations\n\n";
                }
                // compute hazard and survival function estimates
                surv_obj_.est_haz_surv();
                break;
            }
        } // end of the EM algorithm
        if (control_.verbose_ > 0 && n_iter_ == control_.max_iter_) {
            msg("Rearched maximum number of iterations");
        }
        // standardize
        if (surv_standardize0_) {
            surv_obj_.control_.set_standardize(surv_standardize0_);
            surv_obj_.rescale_estimates();
            surv_obj_.est_haz_surv();
        }
        if (cure_standardize0_) {
            cure_obj_.control_.set_standardize(cure_standardize0_);
            cure_obj_.rescale_coef();
        }
        // reset surv_obj_ and cure_obj_ in case of further usage
        surv_obj_.reset_offset_haz();
        // surv_obj_.set_offset(offset0);
        cure_obj_.update_y(surv_obj_.event_);
        // prepare outputs
        surv_coef_ = surv_obj_.coef_;
        cure_coef_ = cure_obj_.coef_;
        unique_time_ = surv_obj_.unique_time_;
        h0_est_ = surv_obj_.h0_est_;
        H0_est_ = surv_obj_.H0_est_;
        S0_est_ = surv_obj_.S0_est_;
        hc_est_ = surv_obj_.hc_est_;
        Hc_est_ = surv_obj_.Hc_est_;
        Sc_est_ = surv_obj_.Sc_est_;
        neg_ll_ = - obs_log_likelihood();
        surv_coef_df_ = compute_coef_df(surv_obj_.coef_);
        cure_coef_df_ = compute_coef_df(cure_obj_.coef_);
        mar_coef_df_ = compute_coef_df(mar_obj_.coef_);
        coef_df_ = surv_coef_df_ + cure_coef_df_ + mar_coef_df_;
        compute_bic1();
        compute_bic2();
        compute_aic();
        // prepare scores and prob in their original order
        const arma::uvec& rev_ord { surv_obj_.rev_ord_ };
        surv_xbeta_ = surv_obj_.get_xbeta().elem(rev_ord);
        cure_xbeta_ = cure_obj_.get_xbeta().elem(rev_ord);
        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj_.predict() };
        susceptible_prob_ = p_vec_event.elem(rev_ord);
        p_vec_event.elem(case1_ind_).ones();
        // compute posterior probabilities from E-step
        for (size_t j: case2_ind_) {
            double numer_j { p_vec(j) *  surv_obj_.S_time_(j)};
            estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        for (size_t j: case3_ind_) {
            double m12_common {
                p_vec(j) * surv_obj_.S_time_(j) * surv_obj_.Sc_time_(j)
            };
            double m1 { surv_obj_.h_time_(j) * m12_common };
            double m2 { surv_obj_.hc_time_(j) * m12_common };
            double m3 {
                (1 - p_vec(j)) * surv_obj_.hc_time_(j) * surv_obj_.Sc_time_(j)
            };
            double m { m1 + m2 + m3 };
            double w1 { m1 / m };
            double w2 { m2 / m };
            set_pmin_bound(w1, cure_obj_.control_.pmin_);
            set_pmin_bound(w2, cure_obj_.control_.pmin_);
            estep_m(j) = w1 + w2;
            estep_event(j) = w1;
        }
        estep_cured_ = 1 - estep_m.elem(rev_ord);
        estep_event_ = event(rev_ord);
        estep_censor_ = 1 - estep_cured_ - estep_event_;
        // compute weighted c-index over the certain records
        c_index_ = Concordance(
            surv_obj_.time_.elem(cer_ind_),
            surv_obj_.event_.elem(cer_ind_),
            surv_xbeta_.elem(cer_ind_),
            p_vec_event.elem(cer_ind_)
            ).index_;
    }


    // fit regularized Cox cure model with uncertain events
    // with adaptive elastic net penalty for perticular lambda's
    // lambda_1 * lasso * factors + lambda_2 * ridge
    inline void CoxphCureMar::net_fit()
    {
        // start
        surv_obj_.set_start();
        cure_obj_.set_start();
        // penalty factor
        surv_obj_.set_penalty_factor();
        cure_obj_.set_penalty_factor();
        // max l1 lambda
        surv_obj_.set_l1_lambda_max();
        cure_obj_.set_l1_lambda_max();
        // early stop: return lambda_max if max_iter = 0
        if (control_.max_iter_ == 0) {
            return;
        }
        // set the start estimates
        arma::vec surv_beta { surv_obj_.control_.start_ };
        arma::vec cure_beta { cure_obj_.control_.start_ };
        surv_obj_.coef0_ = surv_obj_.coef_ = surv_beta;
        cure_obj_.coef0_ = cure_obj_.coef_ = cure_beta;
        // get pre-processed design matrix, time, and event
        const arma::vec& time { surv_obj_.time_ };
        const arma::vec& event { surv_obj_.event_ };
        arma::vec surv_exp_x_beta = arma::exp(surv_obj_.x_ * surv_beta);
        // initialize baseline hazard functions for the E-step
        // for events
        NelsonAalen nelen_event {
            NelsonAalen(time.elem(cer_ind_),
                        event.elem(cer_ind_))
        };
        surv_obj_.h0_time_ = nelen_event.step_inst_rate(time);
        surv_obj_.H0_time_ = nelen_event.step_cum_rate(time);
        // for censoring
        NelsonAalen nelen_censor {
            NelsonAalen(time.elem(cer_ind_),
                        1 - event.elem(cer_ind_))
        };
        surv_obj_.hc_time_ = nelen_censor.step_inst_rate(time);
        surv_obj_.Hc_time_ = nelen_censor.step_cum_rate(time);
        // update related function values
        surv_obj_.h_time_ = surv_obj_.h0_time_ % surv_exp_x_beta;
        surv_obj_.H_time_ = surv_obj_.H0_time_ % surv_exp_x_beta;
        surv_obj_.S0_time_ = arma::exp(- surv_obj_.H0_time_);
        surv_obj_.S_time_ = arma::exp(- surv_obj_.H_time_);
        surv_obj_.Sc_time_ = arma::exp(- surv_obj_.Hc_time_);
        // intialization for the main loop
        arma::vec p_vec { arma::zeros(n_obs_) };
        arma::vec estep_m { event }, estep_event { event };
        double obs_ell {0};
        double reg_obj {0}, reg_obj_old { arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;
        // prepare for tail completion
        const arma::uvec case23_ind { vec_union(case2_ind_, case3_ind_) };
        double max_event_time { time(max_event_time_ind_) };
        if (control_.tail_tau_ < 0)
            control_.tail_tau_ = arma::datum::inf;
        // prepare for exponential tail completion method
        double s0_tau {0}, etail_lambda {0};
        // allow users to stop here
        Rcpp::checkUserInterrupt();
        n_iter_ = 0;
        // main loop of EM algorithm
        for (size_t i {0}; i < control_.max_iter_; ++i) {
            n_iter_ = i + 1;
            // update to the latest estimates
            p_vec = cure_obj_.predict();
            surv_obj_.compute_haz_surv_time();
            surv_obj_.compute_censor_haz_surv_time();
            // record estimates from last step
            surv_beta = surv_obj_.coef_;
            cure_beta = cure_obj_.coef_;
            s0_wi_tail = surv_obj_.S0_time_;
            s_wi_tail = surv_obj_.S_time_;
            // set them as starting values
            surv_obj_.control_.set_start(surv_beta);
            cure_obj_.control_.set_start(cure_beta);
            // tail completion for the conditional survival function
            switch(control_.tail_completion_) {
                case 0:
                    for (size_t j: case2_ind_) {
                        // tail completion after the given tail_tau
                        // by default, it means no tail completion
                        if (time(j) > control_.tail_tau_) {
                            surv_obj_.S_time_(j) = 0.0;
                            surv_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 1:
                    for (size_t j: case2_ind_) {
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            surv_obj_.S_time_(j) = 0.0;
                            surv_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 2:
                    s0_tau = surv_obj_.S0_time_(max_event_time_ind_);
                    etail_lambda = - std::log(s0_tau / max_event_time);
                    for (size_t j: case2_ind_) {
                        // exponential tail by Peng (2003)
                        arma::vec surv_xbeta { surv_obj_.get_xbeta() };
                        if (time(j) > max_event_time) {
                            surv_obj_.S0_time_(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            surv_obj_.S_time_(j) = std::pow(
                                surv_obj_.S0_time_(j),
                                std::exp(surv_xbeta(j))
                                );
                        }
                    }
                    break;
                default:    // do nothing, otherwise
                    break;
            }
            // E-step: compute the v vector for case 2
            for (size_t j: case2_ind_) {
                double numer_j { p_vec(j) * surv_obj_.S_time_(j)};
                estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
            }
            // E-step: compute the w vector for case 3
            for (size_t j: case3_ind_) {
                double m12_common {
                    p_vec(j) * surv_obj_.S0_time_(j) * surv_obj_.Sc_time_(j)
                };
                double m1 { surv_obj_.h_time_(j) * m12_common };
                double m2 { surv_obj_.hc_time_(j) * m12_common };
                double m3 {
                    (1 - p_vec(j)) * surv_obj_.hc_time_(j) *
                    surv_obj_.Sc_time_(j)
                };
                double m { m1 + m2 + m3 };
                double w1 { m1 / m };
                double w2 { m2 / m };
                // some special care for subjects in case 3
                // since event(j) cannot be either 0 or 1!
                set_pmin_bound(w1, cure_obj_.control_.pmin_);
                set_pmin_bound(w2, cure_obj_.control_.pmin_);
                estep_m(j) = w1 + w2;
                estep_event(j) = w1;
            }
            Rcpp::checkUserInterrupt();
            // M-step for the survival layer
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "Running the M-step for the survival layer:"
                            << "\n";
            }
            surv_obj_.set_offset_haz(arma::log(estep_m));
            surv_obj_.update_event_weight(estep_event);
            surv_obj_.net_fit();
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "The M-step for the survival layer was done."
                            << "\n";
            }
            Rcpp::checkUserInterrupt();
            // M-step for the Cure layer
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "Running M-step for the cure layer:"
                            << "\n";
            }
            cure_obj_.update_y(estep_m);
            cure_obj_.net_fit();
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "The M-step for the cure layer was done."
                            << "\n";
            }
            // check if has any `nan`
            if (surv_obj_.coef_.has_nan() || cure_obj_.coef_.has_nan()) {
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function might go to infinite.\n"
                            << "Try to return etimates from last step.\n";
                // take the estimates from the last step
                surv_obj_.coef_ = surv_beta;
                cure_obj_.coef_ = cure_beta;
                surv_obj_.compute_haz_surv_time();
                surv_obj_.S0_time_ = s0_wi_tail;
                surv_obj_.S_time_ = s_wi_tail;
                // surv_obj_.compute_censor_haz_surv_time();
                surv_obj_.est_haz_surv();
                break;
            }
            // update tolerance
            tol1 = rel_l1_norm(surv_obj_.coef_, surv_beta);
            tol2 = rel_l1_norm(cure_obj_.coef_, cure_beta);
            if (control_.verbose_ > 0) {
                // compuete the regularized objective function
                double reg_surv { surv_obj_.net_penalty() };
                double reg_cure { cure_obj_.net_penalty() };
                // compute observed data log-likelihood
                obs_ell = 0;
                // for case 1
                for (size_t j: case1_ind_) {
                    obs_ell += std::log(mar_prob_(j)) +
                        std::log(p_vec(j)) +
                        std::log(surv_obj_.h_time_(j)) -
                        surv_obj_.H_time_(j) - surv_obj_.Hc_time_(j);
                }
                // for case 2
                for (size_t j: case2_ind_) {
                    obs_ell += std::log(mar_prob_(j)) +
                        std::log(p_vec(j) * surv_obj_.S_time_(j) +
                                 1 - p_vec(j)) -
                        surv_obj_.Hc_time_(j) + std::log(surv_obj_.hc_time_(j));
                }
                // for case 3
                for (size_t j: case3_ind_) {
                    double m12_common {
                        p_vec(j) * surv_obj_.S0_time_(j) * surv_obj_.Sc_time_(j)
                    };
                    double m1 { surv_obj_.h_time_(j) * m12_common };
                    double m2 { surv_obj_.hc_time_(j) * m12_common };
                    double m3 {
                        (1 - p_vec(j)) * surv_obj_.hc_time_(j) *
                        surv_obj_.Sc_time_(j)
                    };
                    obs_ell += std::log(mar_prob_(j)) +
                        std::log(m1 + m2 + m3);
                }
                reg_obj = - obs_ell / dn_obs_ + reg_surv + reg_cure;
                Rcpp::Rcout << "\n"
                            << std::string(50, '=')
                            << "\niteration: "
                            << n_iter_;
                if (control_.verbose_ > 2) {
                    Rcpp::Rcout << "\n  surv coef: "
                                << arma2rvec(surv_obj_.coef_)
                                << "\n    relative diff: "
                                << tol1
                                << "\n  cure coef: "
                                << arma2rvec(cure_obj_.coef_)
                                << "\n    relative diff: "
                                << tol2;
                }
                Rcpp::Rcout << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << "\n  regularized objective function: "
                            << reg_obj
                            << "\n    penalty on surv part: "
                            << reg_surv
                            << "\n    penalty on cure part: "
                            << reg_cure
                            << "\n";
                if (reg_obj > reg_obj_old) {
                    Rcpp::Rcout << "Warning: the objective function"
                                << " somehow increased.\n";
                }
                reg_obj_old = reg_obj;
            }
            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();
            // check convergence
            if (tol1 < control_.epsilon_ && tol2 < control_.epsilon_) {
                if (control_.verbose_ > 0) {
                    Rcpp::Rcout << "\n"
                                << std::string(50, '=')
                                << "\n"
                                << "reached convergence after "
                                << n_iter_
                                << " iterations\n\n";
                }
                // compute hazard and survival function estimates
                surv_obj_.est_haz_surv();
                break;
            }
        } // end of the EM algorithm
        if (control_.verbose_ > 0 && n_iter_ == control_.max_iter_) {
            msg("Rearched maximum number of iterations");
        }
        // standardize
        if (surv_standardize0_) {
            surv_obj_.control_.set_standardize(surv_standardize0_);
            surv_obj_.rescale_estimates();
            surv_obj_.est_haz_surv();
        }
        if (cure_standardize0_) {
            cure_obj_.control_.set_standardize(cure_standardize0_);
            cure_obj_.rescale_coef();
        }
        // reset surv_obj_ and cure_obj_ in case of further usage
        surv_obj_.reset_offset_haz();
        // surv_obj_.set_offset(offset0);
        cure_obj_.update_y(surv_obj_.event_);
        // prepare outputs
        surv_coef_ = surv_obj_.coef_;
        cure_coef_ = cure_obj_.coef_;
        unique_time_ = surv_obj_.unique_time_;
        h0_est_ = surv_obj_.h0_est_;
        H0_est_ = surv_obj_.H0_est_;
        S0_est_ = surv_obj_.S0_est_;
        hc_est_ = surv_obj_.hc_est_;
        Hc_est_ = surv_obj_.Hc_est_;
        Sc_est_ = surv_obj_.Sc_est_;
        neg_ll_ = - obs_log_likelihood();
        surv_coef_df_ = compute_coef_df(surv_obj_.coef_);
        cure_coef_df_ = compute_coef_df(cure_obj_.coef_);
        mar_coef_df_ = compute_coef_df(mar_obj_.coef_);
        coef_df_ = surv_coef_df_ + cure_coef_df_ + mar_coef_df_;
        compute_bic1();
        compute_bic2();
        compute_aic();
        // prepare scores and prob in their original order
        const arma::uvec& rev_ord { surv_obj_.ord_ };
        surv_xbeta_ = surv_obj_.get_xbeta().elem(rev_ord);
        cure_xbeta_ = cure_obj_.get_xbeta().elem(rev_ord);
        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj_.predict() };
        susceptible_prob_ = p_vec_event.elem(rev_ord);
        p_vec_event.elem(case1_ind_).ones();
        // compute posterior probabilities from E-step
        for (size_t j: case2_ind_) {
            double numer_j { p_vec(j) *  surv_obj_.S_time_(j)};
            estep_m(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        for (size_t j: case3_ind_) {
            double m12_common {
                p_vec(j) * surv_obj_.S0_time_(j) * surv_obj_.Sc_time_(j)
            };
            double m1 { surv_obj_.h_time_(j) * m12_common };
            double m2 { surv_obj_.hc_time_(j) * m12_common };
            double m3 {
                (1 - p_vec(j)) * surv_obj_.hc_time_(j) * surv_obj_.Sc_time_(j)
            };
            double m { m1 + m2 + m3 };
            double w1 { m1 / m };
            double w2 { m2 / m };
            set_pmin_bound(w1, cure_obj_.control_.pmin_);
            set_pmin_bound(w2, cure_obj_.control_.pmin_);
            estep_m(j) = w1 + w2;
            estep_event(j) = w1;
        }
        estep_cured_ = 1 - estep_m.elem(rev_ord);
        estep_event_ = event(rev_ord);
        estep_censor_ = 1 - estep_cured_ - estep_event_;
        // compute weighted c-index over the certain records
        c_index_ = Concordance(
            surv_obj_.time_.elem(cer_ind_),
            surv_obj_.event_.elem(cer_ind_),
            surv_xbeta_.elem(cer_ind_),
            p_vec_event.elem(cer_ind_)
            ).index_;
    }

    // function to compute the observe data log-likelihood function
    // for given fitted model and estimates
    inline double CoxphCureMar::obs_log_likelihood() const
    {
        double obs_ell { 0 };
        arma::vec sus_prob { cure_obj_.predict() };
        // for case 1
        for (size_t j: case1_ind_) {
            obs_ell += std::log(mar_prob_(j)) +
                std::log(sus_prob(j)) +
                std::log(surv_obj_.h_time_(j)) -
                surv_obj_.H_time_(j) - surv_obj_.Hc_time_(j);
        }
        // for case 2
        for (size_t j: case2_ind_) {
            obs_ell += std::log(mar_prob_(j)) +
                std::log(sus_prob(j) * surv_obj_.S_time_(j) + 1 - sus_prob(j)) -
                surv_obj_.Hc_time_(j) + std::log(surv_obj_.hc_time_(j));
        }
        // for case 3
        for (size_t j: case3_ind_) {
            double m12_common {
                sus_prob(j) * surv_obj_.S0_time_(j) * surv_obj_.Sc_time_(j)
            };
            double m1 { surv_obj_.h_time_(j) * m12_common };
            double m2 { surv_obj_.hc_time_(j) * m12_common };
            double m3 {
                (1 - sus_prob(j)) * surv_obj_.hc_time_(j) *
                surv_obj_.Sc_time_(j)
            };
            obs_ell += std::log(1 - mar_prob_(j)) +
                std::log(m1 + m2 + m3);
        }
        return obs_ell;
    }

    // for given fitted model and a new set of data
    inline double CoxphCureMar::obs_log_likelihood(
        // all the inputs will be sorted inside of the function
        // so their copies are asked here
        arma::vec new_time,
        arma::vec new_event,
        arma::mat new_surv_x,
        arma::mat new_cure_x,
        arma::vec new_surv_offset = 0,
        arma::vec new_cure_offset = 0
        ) const
    {
        // check if the number of covariates matchs the fitted model
        if (new_surv_x.n_cols != surv_p_) {
            throw std::range_error(
                "The number of columns ('new_surv_x') must match the model."
                );
        }
        if (new_cure_x.n_cols != cure_p0_) {
            throw std::range_error(
                "The number of columns ('new_cure_x') must match the model."
                );
        }
        // number of observations
        unsigned int new_n_obs { new_surv_x.n_rows };
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
        new_surv_x = new_surv_x.rows(ord);
        new_cure_x = new_cure_x.rows(ord);
        // process offset terms
        if (new_surv_offset.n_elem != new_surv_x.n_rows) {
            new_surv_offset = arma::zeros(new_surv_x.n_rows);
        } else {
            new_surv_offset = new_surv_offset(ord);
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
        arma::vec new_surv_xbeta {
            mat2vec(new_surv_x * surv_coef_) + new_surv_offset
        };
        arma::vec exp_surv_xbeta { arma::exp(new_surv_xbeta) };
        h_vec %= exp_surv_xbeta;
        H_vec %= exp_surv_xbeta;
        S_vec = arma::exp(- H_vec);
        arma::vec new_cure_xgamma {
            mat2vec(new_cure_x * cure_coef_) + new_cure_offset
        };
        arma::vec p_vec { 1 / (1 + arma::exp(- new_cure_xgamma)) };
        set_pmin_bound(p_vec, cure_obj_.control_.pmin_);
        // for case 1
        for (size_t j: new_case1_ind) {
            obs_ell += std::log(mar_prob_(j)) +
                std::log(p_vec(j)) +
                std::log(h_vec(j)) -
                surv_obj_.H_time_(j) - surv_obj_.Hc_time_(j);
        }
        // for case 2
        for (size_t j: new_case2_ind) {
            obs_ell += std::log(mar_prob_(j)) +
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
            obs_ell += std::log(mar_prob_(j)) +
                std::log(m1 + m2 + m3);
        }
        return obs_ell;
    }

    inline double CoxphCureMar::obs_log_likelihood(
        const CoxphCureMar& new_object) const
    {
        return obs_log_likelihood(new_object.surv_obj_.time_,
                                  new_object.surv_obj_.event_,
                                  new_object.surv_obj_.get_x(true, true),
                                  new_object.cure_obj_.get_x(true, false),
                                  new_object.surv_obj_.control_.offset_,
                                  new_object.cure_obj_.control_.offset_);
    }

}


#endif
