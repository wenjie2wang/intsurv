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

#ifndef INTSURV_COXPH_CURE_H
#define INTSURV_COXPH_CURE_H

#include <RcppArmadillo.h>
#include <string>
#include "assessment.h"
#include "coxph_reg.h"
#include "logistic_reg.h"
#include "utils.h"


namespace Intsurv {

    class CoxphCure
    {
    protected:
        // caches
        unsigned int cure_p0_;  // coef df of cure part wo intercept
        arma::uvec case1_ind_;
        arma::uvec case2_ind_;
        unsigned int max_event_time_ind_; // index of the maximum event time
        double dn_obs_;         // double version of n_obs_

    public:
        CoxphReg cox_obj_;
        LogisticReg cure_obj_;
        bool cox_standardize_;
        bool cure_standardize_;

        // outputs
        unsigned int cox_p_;    // coef df of cox part
        unsigned int cure_p_;   // coef df of cure part wi intercept
        arma::vec cox_coef_;
        arma::vec cure_coef_;
        arma::vec cox_xbeta_;
        arma::vec cure_xbeta_;

        // for each subject and in the original order of X
        arma::vec susceptible_prob_; // probability of being susceptible
        // values in the last E-step
        arma::vec estep_cured_;
        arma::vec estep_susceptible_;
        unsigned int n_obs_;    // number of observations
        unsigned int n_event_;  // number of events
        unsigned int cox_coef_df_;
        unsigned int cure_coef_df_;
        unsigned int coef_df_;
        double bic1_;
        double bic2_;
        double aic_;
        double neg_ll_;
        double c_index_;

        double pmin_ { 1e-5 };
        unsigned int n_iter_;   // number of iterations

        // tail completion
        unsigned int tail_completion_ = 1;
        double tail_tau_ = - 1.0;

        // hazard and survival function estimates at unique time
        arma::vec unique_time_;
        arma::vec h0_est_;
        arma::vec H0_est_;
        arma::vec S0_est_;

        // regularized by particular lambdas
        double cox_l1_lambda_;
        double cox_l2_lambda_;
        double cure_l1_lambda_;
        double cure_l2_lambda_;

        // default constructor
        CoxphCure() {}

        // constructors
        CoxphCure(const arma::vec& time,
                  const arma::vec& event,
                  const arma::mat& cox_x,
                  const arma::mat& cure_x,
                  const bool cure_intercept = true,
                  const bool cox_standardize = true,
                  const bool cure_standardize = true,
                  const arma::vec& cox_offset = arma::vec(),
                  const arma::vec& cure_offset = arma::vec())
        {
            // create the CoxphReg object
            cox_obj_ = CoxphReg(time, event, cox_x, cox_standardize);
            cox_obj_.set_offset(cox_offset, false);
            // pre-process x and y
            cox_p_ = cox_x.n_cols;
            cure_p0_ = cure_x.n_cols;
            cure_p_ = cure_p0_ +
                static_cast<unsigned int>(cure_intercept);
            n_obs_ = cox_x.n_rows;
            dn_obs_ = static_cast<double>(n_obs_);
            const arma::uvec& cox_sort_ind { cox_obj_.ord_ };
            const arma::vec& s_event { cox_obj_.event_ };
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
            case1_ind_ = arma::find(s_event > 0);
            case2_ind_ = arma::find(s_event < 1);
            n_event_ = case1_ind_.n_elem;
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

        inline void set_pmin(const double pmin = 1e-5)
        {
            cure_obj_.set_pmin(pmin);
            pmin_ = pmin;
        }

        // function members
        // fit the Cox cure mode by EM algorithm
        inline void fit(
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const unsigned int max_iter,
            const double epsilon,
            const unsigned int cox_max_iter,
            const double cox_epsilon,
            const unsigned int cure_max_iter,
            const double cure_epsilon,
            const unsigned int tail_completion,
            const double tail_tau,
            const unsigned int verbose
            );

        // fit regularized Cox cure model with adaptive elastic net penalty
        // for perticular lambda's
        inline void net_fit(
            const double cox_l1_lambda,
            const double cox_l2_lambda,
            const double cure_l1_lambda,
            const double cure_l2_lambda,
            const arma::vec& cox_penalty_factor,
            const arma::vec& cure_penalty_factor,
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const bool cox_varying_active,
            const bool cure_varying_active,
            const unsigned int max_iter,
            const double epsilon,
            const unsigned int cox_max_iter,
            const double cox_epsilon,
            const unsigned int cure_max_iter,
            const double cure_epsilon,
            const unsigned int tail_completion,
            const double tail_tau,
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
            bic1_ = std::log(dn_obs_) * coef_df_ + 2 * neg_ll_;
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
    inline void CoxphCure::fit(
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int max_iter = 200,
        const double epsilon = 1e-4,
        const unsigned int cox_max_iter = 100,
        const double cox_epsilon = 1e-4,
        const unsigned int cure_max_iter = 100,
        const double cure_epsilon = 1e-4,
        const unsigned int tail_completion = 1,
        const double tail_tau = -1,
        const unsigned int verbose = 0
        )
    {
        // initialize cox_beta
        const arma::vec& time { cox_obj_.time_ };
        const arma::vec& event { cox_obj_.event_ };
        arma::vec cox_beta { arma::zeros(cox_p_) };
        if (cox_start.n_elem == cox_p_) {
            cox_beta = cox_obj_.gen_start(cox_start);
        } else {
            CoxphReg tmp_object {
                time.elem(case1_ind_),
                event.elem(case1_ind_),
                cox_obj_.x_.rows(case1_ind_),
                false
            };
            tmp_object.fit(cox_beta, cox_max_iter,
                           cox_epsilon, 0);
            cox_beta = tmp_object.coef_;
        }
        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(cure_p_) };
        if (cure_start.n_elem == cure_p_) {
            cure_beta = cure_start;
        } else {
            cure_obj_.fit(cure_beta, cure_max_iter,
                          cure_epsilon, 0);
            cure_beta = cure_obj_.coef_;
        }
        cox_obj_.coef_ = cox_beta;
        cure_obj_.coef_ = cure_beta;
        // initialization
        arma::vec p_vec { arma::zeros(n_obs_) };
        arma::vec estep_v { arma::ones(n_obs_) };
        double obs_ell {0.0}, obs_ell_old { - arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;
        // prepare for tail completion
        double max_event_time { time(max_event_time_ind_) };
        if (tail_tau < 0)
            tail_tau_ = arma::datum::inf;
        tail_completion_ = tail_completion;
        // for exp tail completion only
        double s0_tau, etail_lambda;
        // allow users to stop before the main loop
        Rcpp::checkUserInterrupt();
        // main loop of EM algorithm
        for (size_t i {0}; i < max_iter; ++i) {
            n_iter_ = i + 1;
            // update to the latest estimates
            p_vec = cure_obj_.predict();
            cox_obj_.compute_haz_surv_time();
            // record estimates from last step
            cox_beta = cox_obj_.coef_;
            cure_beta = cure_obj_.coef_;
            s0_wi_tail = cox_obj_.S0_time_;
            s_wi_tail = cox_obj_.S_time_;
            // tail completion for the conditional survival function
            switch(tail_completion_) {
                case 0:
                    for (size_t j: case2_ind_) {
                        // tail completion after the given tail_tau
                        // by default, it means no tail completion
                        if (time(j) > tail_tau_) {
                            cox_obj_.S_time_(j) = 0.0;
                            cox_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 1:
                    for (size_t j: case2_ind_) {
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            cox_obj_.S_time_(j) = 0.0;
                            cox_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 2:
                    s0_tau = cox_obj_.S0_time_(max_event_time_ind_);
                    etail_lambda = - std::log(s0_tau / max_event_time);
                    for (size_t j: case2_ind_) {
                        // exponential tail by Peng (2003)
                        arma::vec cox_xbeta { cox_obj_.get_xbeta() };
                        if (time(j) > max_event_time) {
                            cox_obj_.S0_time_(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            cox_obj_.S_time_(j) = std::pow(
                                cox_obj_.S0_time_(j),
                                std::exp(cox_xbeta(j))
                                );
                        }
                    }
                    break;
                default:    // do nothing, otherwise
                    break;
            }
            // E-step: compute v vector
            for (size_t j: case2_ind_) {
                double numer_j { p_vec(j) * cox_obj_.S_time_(j) };
                estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
            }
            Rcpp::checkUserInterrupt();
            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\nRunning M-step for the survival layer:";
            }
            cox_obj_.set_offset_haz(arma::log(estep_v));
            cox_obj_.fit(cox_beta, cox_max_iter, cox_epsilon,
                         less_verbose(verbose, 3));
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "The M-step for the survival layer was done.\n";
            }
            Rcpp::checkUserInterrupt();
            // M-step for the Cure layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\nRunning M-step for the cure layer:";
            }
            cure_obj_.update_y(estep_v);
            cure_obj_.fit(cure_beta, cure_max_iter, cure_epsilon,
                          less_verbose(verbose, 3));
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done.\n";
            }
            // check if has any `nan`
            if (cox_obj_.coef_.has_nan() || cure_obj_.coef_.has_nan()) {
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function might go to infinite.\n"
                            << "Try to return etimates from last step.\n";
                // take the estimates from the last step
                cox_obj_.coef_ = cox_beta;
                cure_obj_.coef_ = cure_beta;
                cox_obj_.compute_haz_surv_time();
                cox_obj_.S0_time_ = s0_wi_tail;
                cox_obj_.S_time_ = s_wi_tail;
                // cox_obj_.compute_censor_haz_surv_time();
                cox_obj_.est_haz_surv();
                break;
            }
            // update tolerance
            tol1 = rel_l1_norm(cox_obj_.coef_, cox_beta);
            tol2 = rel_l1_norm(cure_obj_.coef_, cure_beta);
            if (verbose > 0) {
                // compute observed log-likelihood
                obs_ell = 0;
                // for case 1
                for (size_t j: case1_ind_) {
                    obs_ell += std::log(p_vec(j)) +
                        std::log(cox_obj_.h_time_(j)) -
                        cox_obj_.H_time_(j);
                }
                // for case 2
                for (size_t j: case2_ind_) {
                    obs_ell += std::log(
                        p_vec(j) * cox_obj_.S_time_(j) + (1 - p_vec(j))
                        );
                }
                Rcpp::Rcout << "\n"
                            << std::string(50, '=')
                            << "\niteration: "
                            << n_iter_;
                if (verbose > 2) {
                    Rcpp::Rcout << "\n  Cox coef: "
                                << arma2rvec(cox_obj_.coef_)
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
            if (tol1 < epsilon && tol2 < epsilon) {
                if (verbose > 0) {
                    Rcpp::Rcout << "\n"
                                << std::string(50, '=')
                                << "\n"
                                << "reached convergence after "
                                << n_iter_
                                << " iterations\n\n";
                }
                // compute hazard and survival function estimates
                cox_obj_.est_haz_surv();
                break;
            }
        } // end of the EM algorithm
        if (verbose > 0 && n_iter_ == max_iter) {
            msg("Rearched maximum number of iterations");
        }
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
        unique_time_ = cox_obj_.unique_time_;
        h0_est_ = cox_obj_.h0_est_;
        H0_est_ = cox_obj_.H0_est_;
        S0_est_ = cox_obj_.S0_est_;
        neg_ll_ = - obs_log_likelihood();
        cox_coef_df_ = compute_coef_df(cox_obj_.coef_);
        cure_coef_df_ = compute_coef_df(cure_obj_.coef_);
        coef_df_ = cox_coef_df_ + cure_coef_df_;
        compute_bic1();
        compute_bic2();
        compute_aic();

        // prepare scores and prob in their original order
        const arma::uvec& rev_ord { cox_obj_.rev_ord_ };
        cox_xbeta_ = cox_obj_.get_xbeta().elem(rev_ord);
        cure_xbeta_ = cure_obj_.get_xbeta().elem(rev_ord);
        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj_.predict() };
        susceptible_prob_ = p_vec_event.elem(rev_ord);
        p_vec_event.elem(case1_ind_).ones();
        // compute posterior probabilities from E-step
        for (size_t j: case2_ind_) {
            double numer_j { p_vec(j) * cox_obj_.S_time_(j) };
            estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        estep_susceptible_ = estep_v.elem(rev_ord);
        estep_cured_ = 1 - estep_susceptible_;
        // compute weighted c-index
        c_index_ = Intsurv::Concordance(
            time, event, cox_xbeta_, p_vec_event
            ).index_;
    }

    // fit regularized Cox cure model with adaptive elastic net penalty
    // for a perticular lambda
    // lambda_1 * lasso * factors + lambda_2 * ridge
    inline void CoxphCure::net_fit(
        const double cox_l1_lambda,
        const double cox_l2_lambda,
        const double cure_l1_lambda,
        const double cure_l2_lambda,
        const arma::vec& cox_penalty_factor = arma::vec(),
        const arma::vec& cure_penalty_factor = arma::vec(),
        const arma::vec& cox_start = arma::vec(),
        const arma::vec& cure_start = arma::vec(),
        const bool cox_varying_active = true,
        const bool cure_varying_active = true,
        const unsigned int max_iter = 200,
        const double epsilon = 1e-4,
        const unsigned int cox_max_iter = 100,
        const double cox_epsilon = 1e-4,
        const unsigned int cure_max_iter = 100,
        const double cure_epsilon = 1e-4,
        const unsigned int tail_completion = 1,
        const double tail_tau = -1,
        const unsigned int verbose = 0
        )
    {
        // penalty factor
        cox_obj_.set_penalty_factor(cox_penalty_factor);
        cure_obj_.set_penalty_factor(cure_penalty_factor);
        // max l1 lambda
        cox_obj_.set_l1_lambda_max();
        cure_obj_.set_l1_lambda_max();
        cox_l1_lambda_ = cox_l1_lambda;
        cure_l1_lambda_ = cure_l1_lambda;
        cox_l2_lambda_ = cox_l2_lambda;
        cure_l2_lambda_ = cure_l2_lambda;
        // early stop: return lambda_max if max_iter = 0
        if (max_iter == 0) {
            return;
        }
        // set the start estimates
        arma::vec cox_beta { cox_obj_.gen_start(cox_start) };
        arma::vec cure_beta { cure_obj_.gen_start(cure_start) };
        cox_obj_.coef_ = cox_beta;
        cure_obj_.coef_ = cure_beta;
        // initialization
        arma::vec p_vec { arma::zeros(n_obs_) };
        const arma::vec& time { cox_obj_.time_ };
        const arma::vec& event { cox_obj_.event_ };
        arma::vec estep_v { event };
        double obs_ell { 0.0 };
        double reg_obj {0}, reg_obj_old { arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;
        // prepare for tail completion
        double max_event_time { time(max_event_time_ind_) };
        if (tail_tau < 0)
            tail_tau_ = arma::datum::inf;
        tail_completion_ = tail_completion;
        // prepare for exponential tail completion method
        double s0_tau {0}, etail_lambda {0};
        // allow users to stop here
        Rcpp::checkUserInterrupt();
        n_iter_ = 0;
        // main loop of EM algorithm
        for (size_t i {0}; i < max_iter; ++i) {
            n_iter_ = i + 1;
            // update to the latest estimates
            p_vec = cure_obj_.predict(cure_obj_.coef_);
            cox_obj_.compute_haz_surv_time();
            // record estimates from last step
            cox_beta = cox_obj_.coef_;
            cure_beta = cure_obj_.coef_;
            s0_wi_tail = cox_obj_.S0_time_;
            s_wi_tail = cox_obj_.S_time_;
            // tail completion for the conditional survival function
            switch(tail_completion_) {
                case 0:
                    for (size_t j: case2_ind_) {
                        // tail completion after the given tail_tau
                        // by default, it means no tail completion
                        if (time(j) > tail_tau_) {
                            cox_obj_.S_time_(j) = 0.0;
                            cox_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 1:
                    for (size_t j: case2_ind_) {
                        // zero-tail constraint
                        if (time(j) > max_event_time) {
                            cox_obj_.S_time_(j) = 0.0;
                            cox_obj_.S0_time_(j) = 0.0;
                        }
                    }
                    break;
                case 2:
                    s0_tau = cox_obj_.S0_time_(max_event_time_ind_);
                    etail_lambda = - std::log(s0_tau / max_event_time);
                    for (size_t j: case2_ind_) {
                        // exponential tail by Peng (2003)
                        arma::vec cox_xbeta { cox_obj_.get_xbeta() };
                        if (time(j) > max_event_time) {
                            cox_obj_.S0_time_(j) = std::exp(
                                - etail_lambda * time(j)
                                );
                            cox_obj_.S_time_(j) = std::pow(
                                cox_obj_.S0_time_(j),
                                std::exp(cox_xbeta(j))
                                );
                        }
                    }
                    break;
                default:    // do nothing, otherwise
                    break;
            }
            // E-step: compute v vector
            for (size_t j: case2_ind_) {
                double numer_j { p_vec(j) * cox_obj_.S_time_(j) };
                estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
            }
            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\nRunning the M-step for the survival layer:";
            }
            cox_obj_.set_offset_haz(arma::log(estep_v));
            cox_obj_.net_fit(cox_l1_lambda_,
                             cox_l2_lambda_,
                             cox_obj_.penalty_factor_,
                             cox_beta,
                             cox_varying_active,
                             cox_max_iter,
                             cox_epsilon,
                             less_verbose(verbose, 3));
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\n"
                            << "The M-step for the survival layer was done.\n";
            }
            Rcpp::checkUserInterrupt();
            // M-step for the Cure layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\nRunning the M-step for the cure layer:";
            }
            cure_obj_.update_y(estep_v);
            cure_obj_.net_fit(cure_l1_lambda_,
                              cure_l2_lambda_,
                              cure_obj_.penalty_factor_,
                              cure_beta,
                              cure_varying_active,
                              cure_max_iter,
                              cure_epsilon,
                              less_verbose(verbose, 3));
            if (verbose > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done.\n";
            }
            // check if has any `nan`
            if (cox_obj_.coef_.has_nan() || cure_obj_.coef_.has_nan()) {
                Rcpp::Rcout << "Warning: Found NA's in coef. "
                            << "The objective function might go to infinite.\n"
                            << "Try to return etimates from last step.\n";
                // take the estimates from the last step
                cox_obj_.coef_ = cox_beta;
                cure_obj_.coef_ = cure_beta;
                cox_obj_.compute_haz_surv_time();
                cox_obj_.S0_time_ = s0_wi_tail;
                cox_obj_.S_time_ = s_wi_tail;
                // cox_obj_.compute_censor_haz_surv_time();
                cox_obj_.est_haz_surv();
                break;
            }
            // update tolerance
            tol1 = rel_l1_norm(cox_obj_.coef_, cox_beta);
            tol2 = rel_l1_norm(cure_obj_.coef_, cure_beta);
            // verbose tracing for objective function
            if (verbose > 0) {
                // compuete the regularized objective function
                double reg_cox { cox_obj_.net_penalty() };
                double reg_cure { cure_obj_.net_penalty() };
                reg_obj = - obs_ell / dn_obs_ + reg_cox + reg_cure;
                Rcpp::Rcout << "\n"
                            << std::string(50, '=')
                            << "\niteration: "
                            << n_iter_;
                if (verbose > 2) {
                    Rcpp::Rcout << "\n  Cox coef: "
                                << arma2rvec(cox_obj_.coef_)
                                << "\n    relative diff: " << tol1
                                << "\n  cure coef: "
                                << arma2rvec(cure_obj_.coef_)
                                << "\n    relative diff: " << tol2;
                }
                Rcpp::Rcout << "\n  observed negative log-likelihood: "
                            << - obs_ell
                            << "\n  regularized objective function: "
                            << reg_obj
                            << "\n    penalty on Cox model: "
                            << reg_cox
                            << "\n    penalty on cure layer: "
                            << reg_cure
                            << "\n";
                if (reg_obj > reg_obj_old) {
                    Rcpp::Rcout << "Warning: the objective function"
                                << " somehow increased.\n";
                }
                reg_obj_old = reg_obj;
            }
            if (tol1 < epsilon && tol2 < epsilon) {
                if (verbose > 0) {
                    Rcpp::Rcout << "\n"
                                << std::string(50, '=')
                                << "\n"
                                << "reached convergence after "
                                << n_iter_
                                << " iterations\n\n";
                }
                cox_obj_.est_haz_surv();
                break;
            }
        } // end of the EM algorithm
        if (verbose > 0 && n_iter_ == max_iter) {
            msg("Rearched maximum number of iterations");
        }
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
        // cox_en_coef_ = cox_obj_.en_coef_;
        // cure_en_coef_ = cure_obj_.en_coef_;

        unique_time_ = cox_obj_.unique_time_;
        h0_est_ = cox_obj_.h0_est_;
        H0_est_ = cox_obj_.H0_est_;
        S0_est_ = cox_obj_.S0_est_;

        neg_ll_ = - obs_ell;
        cox_coef_df_ = compute_coef_df(cox_obj_.coef_);
        cure_coef_df_ = compute_coef_df(cure_obj_.coef_);
        coef_df_ = cox_coef_df_ + cure_coef_df_;
        cox_l1_lambda_ = cox_l1_lambda;
        cox_l2_lambda_ = cox_l2_lambda;
        cure_l1_lambda_ = cure_l1_lambda;
        cure_l2_lambda_ = cure_l2_lambda;

        // compute BIC
        compute_bic1();
        compute_bic2();
        compute_aic();

        // record tail completion
        // tail_completion = tail_completion;
        // tail_tau = tail_tau;

        // prepare scores and prob in their original order
        const arma::uvec& rev_ord { cox_obj_.rev_ord_ };
        cox_xbeta_ = cox_obj_.get_xbeta().elem(rev_ord);
        cure_xbeta_ = cure_obj_.get_xbeta().elem(rev_ord);

        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj_.predict() };
        susceptible_prob_ = p_vec_event.elem(rev_ord);
        p_vec_event.elem(case1_ind_).ones();

        // compute posterior probabilities from E-step
        for (size_t j: case2_ind_) {
            double numer_j { p_vec(j) * cox_obj_.S_time_(j) };
            estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        estep_susceptible_ = estep_v.elem(rev_ord);
        estep_cured_ = 1 - estep_susceptible_;
        // compute weight C-index
        c_index_ = Intsurv::Concordance(
            time, event, cox_xbeta_, p_vec_event
            ).index_;
    }

    // function to compute the observe data log-likelihood function
    // for given fitted model and estimates
    inline double CoxphCure::obs_log_likelihood() const
    {
        double obs_ell { 0 };
        arma::vec sus_prob { cure_obj_.predict() };
        // for case 1
        for (size_t j: case1_ind_) {
            obs_ell += std::log(sus_prob(j)) +
                std::log(cox_obj_.h_time_(j)) -
                cox_obj_.H_time_(j);
        }
        // for case 2
        for (size_t j: case2_ind_) {
            obs_ell += std::log(
                sus_prob(j) * cox_obj_.S_time_(j) + (1 - sus_prob(j))
                );
        }
        return obs_ell;
    }

    // for given fitted model and a new set of data
    inline double CoxphCure::obs_log_likelihood(
        // all the inputs will be sorted inside of the function
        // so their copies are asked here
        arma::vec new_time,
        arma::vec new_event,
        arma::mat new_cox_x,
        arma::mat new_cure_x,
        arma::vec new_cox_offset = arma::vec(),
        arma::vec new_cure_offset = arma::vec()
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
        if (new_cure_x.n_rows != new_n_obs) {
            throw std::range_error(
                "The number of rows of the new data must be the same."
                );
        }

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
        arma::uvec new_case1_ind { arma::find(new_event > 0) };
        arma::uvec new_case2_ind { arma::find(new_event < 1) };
        // construct the baseline survival curve
        // tail completion has already been applied to S0_est_
        arma::vec S0_vec {
            arma::join_cols(arma::ones<arma::vec>(1), S0_est_)
        };
        // baseline estimates
        arma::vec S_vec {
            step_fun(new_time, unique_time_, S0_vec)
        };
        arma::vec H_vec { - arma::log(S_vec) };
        // only consider positive values
        arma::uvec which_h { arma::find(h0_est_ > 0) };
        arma::vec h_vec {
            step_fun2(new_time,
                      unique_time_.elem(which_h),
                      h0_est_.elem(which_h))
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
        set_pmin_bound(p_vec, pmin_);
        // for case 1
        for (size_t j: new_case1_ind) {
            obs_ell += std::log(p_vec(j)) +
                std::log(h_vec(j)) - H_vec(j);
        }
        // for case 2
        for (size_t j: new_case2_ind) {
            obs_ell += std::log(
                p_vec(j) * S_vec(j) + (1 - p_vec(j))
                );
        }
        return obs_ell;
    }

}

#endif
