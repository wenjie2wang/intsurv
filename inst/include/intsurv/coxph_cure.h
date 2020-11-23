//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2020  Wenjie Wang <wang@wwenjie.org>
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

#ifndef COXPH_CURE_H
#define COXPH_CURE_H

#include <RcppArmadillo.h>
#include <string>
#include "assessment.h"
#include "coxph_reg.h"
#include "logistic_reg.h"
#include "utils.h"


namespace Intsurv {

    class CoxphCure {
    private:
        CoxphReg cox_obj;
        LogisticReg cure_obj;
        unsigned int cox_p;     // coef df of cox part
        unsigned int cure_p;    // coef df of cure part wi intercept
        unsigned int cure_p0;   // coef df of cure part wo intercept
        arma::uvec case1_ind;
        arma::uvec case2_ind;
        unsigned int max_event_time_ind; // index of the maximum event time

    public:
        arma::vec cox_coef;
        arma::vec cure_coef;
        unsigned int coef_df;     // degree of freedom of coef estimates
        double negLogL;           // negative log-likelihood
        unsigned int nObs;        // number of observations
        unsigned int nEvent;      // number of events
        unsigned int num_iter;    // number of iterations
        double bic1;              // BIC: log(num_obs) * coef_df + 2 * negLogL
        double bic2;              // BIC: log(num_event) * coef_df + 2 * negLogL
        double aic;               // AIC: 2 * coef_df + 2 * negLogL
        double c_index;           // weighted C-index

        // for each subject and in the original order of X
        arma::vec cox_xBeta;        // score from the survival layer
        arma::vec cure_xBeta;       // score from the cure layer
        arma::vec susceptible_prob; // probability of being susceptible
        // values in the last E-step
        arma::vec estep_cured;
        arma::vec estep_susceptible;

        // tail completion
        // unsigned int tail_completion = 1;
        // double tail_tau = - 1.0;

        // hazard and survival function estimates at unique time
        arma::vec unique_time;
        arma::vec h0_est;
        arma::vec H0_est;
        arma::vec S0_est;

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
        CoxphCure() {}

        // constructors
        CoxphCure(const arma::vec& time,
                  const arma::vec& event,
                  const arma::mat& cox_x,
                  const arma::mat& cure_x,
                  const bool& cure_intercept = true,
                  const bool& cox_standardize = true,
                  const bool& cure_standardize = true)
        {
            // create the CoxphReg object
            this->cox_obj = CoxphReg(time, event, cox_x, cox_standardize);
            // pre-process x and y
            this->cox_p = cox_x.n_cols;
            this->cure_p0 = cure_x.n_cols;
            this->cure_p = this->cure_p0 +
                static_cast<unsigned int>(cure_intercept);
            this->nObs = cox_x.n_rows;
            arma::uvec cox_sort_ind { cox_obj.get_sort_index() };
            arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
            arma::vec s_event { event.elem(cox_sort_ind) };
            this->case1_ind = arma::find(s_event > 0);
            this->case2_ind = arma::find(s_event < 1);
            this->nEvent = this->case1_ind.n_elem;
            this->max_event_time_ind = arma::max(this->case1_ind);
            // create the LogisticReg object
            this->cure_obj = LogisticReg(cure_xx, s_event, cure_intercept,
                                         cure_standardize);
        }

        // function members
        // helper functions
        inline unsigned int get_cox_p() const {
            return this->cox_p;
        }
        inline unsigned int get_cure_p() const {
            return this->cure_p;
        }

        // fit the Cox cure mode by EM algorithm
        inline void fit(
            const arma::vec& cox_start,
            const arma::vec& cure_start,
            const unsigned int& em_max_iter,
            const double& em_rel_tol,
            const unsigned int& cox_mstep_max_iter,
            const double& cox_mstep_rel_tol,
            const unsigned int& cure_mstep_max_iter,
            const double& cure_mstep_rel_tol,
            const bool& firth,
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
            const unsigned int& tail_completion,
            double tail_tau,
            const double& pmin,
            const unsigned int& early_stop,
            const unsigned int& verbose
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
            const double pmin
            ) const;

        // compute BIC
        inline void compute_bic1() {
            this->bic1 = std::log(nObs) * coef_df + 2 * negLogL;
        }
        inline void compute_bic2() {
            this->bic2 = std::log(case1_ind.n_elem) *
                coef_df + 2 * negLogL;
        }
        inline void compute_aic() {
            this->aic = 2 * (coef_df + negLogL);
        }

    };                          // end of class definition


    // fit the Cox cure mode by EM algorithm
    inline void CoxphCure::fit(
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int& em_max_iter = 300,
        const double& em_rel_tol = 1e-5,
        const unsigned int& cox_mstep_max_iter = 200,
        const double& cox_mstep_rel_tol = 1e-5,
        const unsigned int& cure_mstep_max_iter = 200,
        const double& cure_mstep_rel_tol = 1e-5,
        const bool& firth = false,
        const unsigned int& tail_completion = 1,
        double tail_tau = -1,
        const double& pmin = 1e-5,
        const unsigned int& early_stop = 0,
        const unsigned int& verbose = 0
        )
    {
        // initialize cox_beta
        const arma::vec time { cox_obj.get_time() };
        const arma::vec event { cox_obj.get_event() };
        arma::vec cox_beta { arma::zeros(this->cox_p) };
        if (cox_start.n_elem == this->cox_p) {
            cox_beta = cox_start;
        } else {
            CoxphReg tmp_object {
                time.elem(case1_ind),
                event.elem(case1_ind),
                cox_obj.get_x().rows(case1_ind)
            };
            tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
            cox_beta = tmp_object.coef;
        }
        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(this->cure_p) };
        if (cure_start.n_elem == this->cure_p) {
            cure_beta = cure_start;
        } else {
            if (firth) {
                cure_obj.firth_fit(cure_beta, cure_mstep_max_iter,
                                   cure_mstep_rel_tol, pmin);
            } else {
                cure_obj.fit(cure_beta, cure_mstep_max_iter,
                             cure_mstep_rel_tol, pmin);
            }
            cure_beta = cure_obj.coef;
        }
        cox_obj.coef = cox_beta;
        cure_obj.coef = cure_beta;

        // initialization
        arma::vec p_vec { arma::zeros(nObs) };
        arma::vec estep_v { cox_obj.get_event() };
        size_t i {0};
        double obs_ell {0}, obs_ell_old { - arma::datum::inf };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;

        // prepare for tail completion
        double max_event_time { time(this->max_event_time_ind) };
        if (tail_tau < 0)
            tail_tau = arma::datum::inf;

        // allow users to stop the main loop
        Rcpp::checkUserInterrupt();

        // main loop of EM algorithm
        while (true) {

            // update to the latest estimates
            p_vec = cure_obj.predict(cure_obj.coef, pmin);
            cox_obj.compute_haz_surv_time();

            // prepare for exponential tail completion method
            double s0_tau {0}, etail_lambda {0};
            if (tail_completion == 2) {
                s0_tau = cox_obj.S0_time(max_event_time_ind);
                etail_lambda = - std::log(s0_tau / max_event_time);
            }
            // tail completion for case 2
            for (size_t j: case2_ind) {
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

            // compute observed log-likelihood
            obs_ell = 0;
            // for case 1
            for (size_t j: case1_ind) {
                obs_ell += std::log(p_vec(j)) +
                    std::log(cox_obj.h_time(j)) +
                    std::log(cox_obj.S_time(j));
            }
            // for case 2
            for (size_t j: case2_ind) {
                obs_ell += std::log(
                    p_vec(j) * cox_obj.S_time(j) + (1 - p_vec(j))
                    );
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
                break;
            }
            // early exit if the observed data log-likelihood decreased, which
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
                // cox_obj.compute_censor_haz_surv_time();
                cox_obj.est_haz_surv();
                // use old obs likelihood
                obs_ell = obs_ell_old;
                // break here
                break;
            }

            // check convergence
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

            // E-step: compute v vector
            for (size_t j: case2_ind) {
                double numer_j { p_vec(j) * cox_obj.S_time(j) };
                estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning M-step for the survival layer:";
            }
            cox_obj.set_offset(arma::log(estep_v));
            cox_obj.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol,
                        early_stop == 1, verbose > 2);
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
            cure_obj.update_y(estep_v);
            if (firth) {
                cure_obj.firth_fit(cure_beta, cure_mstep_max_iter,
                                   cure_mstep_rel_tol, pmin);
            } else {
                cure_obj.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol,
                             pmin, early_stop == 1, verbose > 2);
            }
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done."
                            << std::endl;
            }

            // update tolerance
            tol1 = rel_l1_norm(cox_obj.coef, cox_beta);
            tol2 = rel_l1_norm(cure_obj.coef, cure_beta);

        } // end of the EM algorithm

        // reset cox_obj and cure_obj in case of further usage
        cox_obj.reset_offset();
        cure_obj.update_y(cox_obj.get_event());

        // prepare outputs
        this->cox_coef = cox_obj.coef;
        this->cure_coef = cure_obj.coef;
        this->unique_time = cox_obj.unique_time;
        this->h0_est = cox_obj.h0_est;
        this->H0_est = cox_obj.H0_est;
        this->S0_est = cox_obj.S0_est;
        this->negLogL = - obs_ell;
        this->coef_df = cox_obj.coef_df + cure_obj.coef_df;
        this->num_iter = i;
        this->compute_bic1();
        this->compute_bic2();
        this->compute_aic();

        // // record tail completion
        // this->tail_completion = tail_completion;
        // this->tail_tau = tail_tau;

        // prepare scores and prob in their original order
        arma::uvec rev_ord { cox_obj.get_rev_sort_index() };
        this->cox_xBeta = cox_obj.xBeta.elem(rev_ord);
        this->cure_xBeta = cure_obj.xBeta.elem(rev_ord);
        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj.prob_vec };
        p_vec_event.elem(case1_ind).ones();
        this->susceptible_prob = cure_obj.prob_vec.elem(rev_ord);
        // compute posterior probabilities from E-step
        for (size_t j: case2_ind) {
            double numer_j { p_vec(j) * cox_obj.S_time(j) };
            estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        this->estep_susceptible = estep_v.elem(rev_ord);
        this->estep_cured = 1 - this->estep_susceptible;
        // compute weighted c-index
        this->c_index = Intsurv::Concordance(
            time, event, cox_obj.xBeta, p_vec_event
            ).index;
    }


    // fit regularized Cox cure model with adaptive elastic net penalty
    // for a perticular lambda
    // lambda_1 * lasso * factors + lambda_2 * ridge
    inline void CoxphCure::regularized_fit(
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
        const unsigned int& tail_completion = 1,
        double tail_tau = -1,
        const double& pmin = 1e-5,
        const unsigned int& early_stop = 0,
        const unsigned int& verbose = 0
        )
    {
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
        arma::vec cox_beta { arma::zeros(cox_p) };
        arma::vec cure_beta { arma::zeros(cure_p) };

        // compute the large enough lambdas that result in all-zero estimates
        arma::vec cox_grad_zero { arma::abs(cox_obj.gradient(cox_beta)) };
        arma::vec cure_grad_zero {
            arma::abs(cure_obj.gradient(cure_beta, pmin))
        };
        cure_grad_zero = cure_grad_zero.tail(cure_l1_penalty.n_elem);
        // excluding variable with zero penalty factor
        arma::uvec cox_active_l1_penalty { arma::find(cox_l1_penalty > 0) };
        arma::uvec cure_active_l1_penalty { arma::find(cure_l1_penalty > 0) };
        this->cox_l1_lambda_max = arma::max(
            cox_grad_zero.elem(cox_active_l1_penalty) /
            cox_l1_penalty.elem(cox_active_l1_penalty)
            ) / this->nObs;
        this->cure_l1_lambda_max = arma::max(
            cure_grad_zero.elem(cure_active_l1_penalty) /
            cure_l1_penalty.elem(cure_active_l1_penalty)
            ) / this->nObs;

        // early stop: return lambda_max if em_max_iter = 0
        if (em_max_iter == 0) {
            return;
        }

        // set the start estimates
        if (cox_start.n_elem == cox_p) {
            cox_beta = cox_start;
        }
        if (cure_start.n_elem == cure_p) {
            cure_beta = cure_start;
        }
        cure_obj.coef = cure_beta;
        cox_obj.coef = cox_beta;
        cure_obj.coef_df = get_coef_df(cure_beta);
        cox_obj.coef_df = get_coef_df(cox_beta);

        // initialization
        arma::vec p_vec { arma::zeros(nObs) };
        const arma::vec time { cox_obj.get_time() };
        const arma::vec event { cox_obj.get_event() };
        arma::vec estep_v { event };
        size_t i {0};
        double obs_ell {0};
        double reg_obj {0};
        double reg_obj_old { arma::datum::inf }, obs_ell_old { reg_obj_old };
        double bic1_old { arma::datum::inf }, bic2_old { bic1_old };
        double tol1 { arma::datum::inf }, tol2 { tol1 };
        arma::vec s0_wi_tail, s_wi_tail;
        bool verbose_mstep { verbose > 2 };

        // prepare for tail completion
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

            // prepare for exponential tail completion method
            double s0_tau {0}, etail_lambda {0};
            if (tail_completion == 2) {
                s0_tau = cox_obj.S0_time(max_event_time_ind);
                etail_lambda = - std::log(s0_tau / max_event_time);
            }
            // tail completion
            for (size_t j: case2_ind) {
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

            // compute observed log-likelihood
            obs_ell = 0;
            // for case 1
            for (size_t j: case1_ind) {
                obs_ell += std::log(p_vec(j)) +
                    std::log(cox_obj.h_time(j)) +
                    std::log(cox_obj.S_time(j));
            }
            // for case 2
            for (size_t j: case2_ind) {
                obs_ell += std::log(
                    p_vec(j) * cox_obj.S_time(j) + (1 - p_vec(j))
                    );
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

            // compute bic
            this->negLogL = - obs_ell;
            this->coef_df = cox_obj.coef_df + cure_obj.coef_df;
            this->compute_bic1();
            this->compute_bic2();
            this->compute_aic();

            // verbose tracing for objective function
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
                    Rcpp::Rcout << "The objective function increased."
                                << std::endl;
                }
                early_exit = early_exit || early_stop == 1;
            }
            // early exit if bic increased
            if (this->bic1 > bic1_old && this->bic2 > bic2_old) {
                if (verbose) {
                    Rcpp::Rcout << "The BIC increased."
                                << std::endl;
                }
                early_exit = early_exit || early_stop == 2;
            }
            // return the estimates from last step
            if (early_exit) {
                if (verbose) {
                    Rcpp::Rcout << "Ended the EM algorithm after iteration "
                                << i
                                << " with estimates from last step."
                                << std::endl;
                }
                // compute hazard and survival function estimates
                cox_obj.coef = cox_beta;
                cure_obj.coef = cure_beta;
                // update coef_df and en_coef
                cox_obj.update_from_coef(cox_l2_lambda);
                cure_obj.update_from_coef(cure_l2_lambda);
                // update hazard and survival function estimates
                cox_obj.compute_censor_haz_surv_time();
                cox_obj.S0_time = s0_wi_tail;
                cox_obj.S_time = s_wi_tail;
                cox_obj.est_haz_surv();
                // convert back obs_ell
                obs_ell = obs_ell_old;
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
            reg_obj_old = reg_obj;
            bic1_old = this->bic1;
            bic2_old = this->bic2;
            s0_wi_tail = cox_obj.S0_time;
            s_wi_tail = cox_obj.S_time;

            // update iter for the next iteration
            ++i;

            // E-step: compute v vector
            for (size_t j: case2_ind) {
                double numer_j { p_vec(j) *  cox_obj.S_time(j)};
                estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
                // special care prevents coef diverging
                // if (estep_v(j) < pmin) {
                //     estep_v(j) = 0;
                // } else if (estep_v(j) > 1 - pmin) {
                //     estep_v(j) = 1;
                // }
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nRunning the M-step for the survival layer:";
            }
            cox_obj.set_offset(arma::log(estep_v));
            cox_obj.regularized_fit(
                cox_l1_lambda, cox_l2_lambda, cox_l1_penalty_factor,
                cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol,
                early_stop == 1, verbose_mstep
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
                            << "\nRunning the M-step for the cure layer:";
            }
            cure_obj.update_y(estep_v);
            cure_obj.regularized_fit(
                cure_l1_lambda, cure_l2_lambda, cure_l1_penalty_factor,
                cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol,
                pmin, early_stop == 1, verbose_mstep
                );
            if (verbose > 1) {
                Rcpp::Rcout << "\n" << std::string(40, '-')
                            << "\nThe M-step for the cure layer was done."
                            << std::endl;
            }

            // update tolerance
            tol1 = rel_l1_norm(cox_obj.coef, cox_beta);
            tol2 = rel_l1_norm(cure_obj.coef, cure_beta);

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

        this->negLogL = - obs_ell;
        this->coef_df = cox_obj.coef_df + cure_obj.coef_df;
        this->cox_l1_lambda = cox_l1_lambda;
        this->cox_l2_lambda = cox_l2_lambda;
        this->cure_l1_lambda = cure_l1_lambda;
        this->cure_l2_lambda = cure_l2_lambda;
        this->num_iter = i;

        // compute BIC
        this->compute_bic1();
        this->compute_bic2();
        this->compute_aic();

        // record tail completion
        // this->tail_completion = tail_completion;
        // this->tail_tau = tail_tau;

        // prepare scores and prob in their original order
        arma::uvec rev_ord { cox_obj.get_rev_sort_index() };
        this->cox_xBeta = cox_obj.xBeta.elem(rev_ord);
        this->cure_xBeta = cure_obj.xBeta.elem(rev_ord);

        // set prob to be 1 for events for computing C-index
        arma::vec p_vec_event { cure_obj.prob_vec };
        p_vec_event.elem(case1_ind).ones();
        this->susceptible_prob = cure_obj.prob_vec.elem(rev_ord);

        // compute posterior probabilities from E-step
        for (size_t j: case2_ind) {
            double numer_j { p_vec(j) * cox_obj.S_time(j) };
            estep_v(j) = numer_j / (1 - p_vec(j) + numer_j);
        }
        this->estep_susceptible = estep_v.elem(rev_ord);
        this->estep_cured = 1 - this->estep_susceptible;
        // compute weight C-index
        this->c_index = Intsurv::Concordance(
            time, event, cox_obj.xBeta, p_vec_event
            ).index;
    }

    // function to compute the observe data log-likelihood function
    // for given fitted model and estimates
    inline double CoxphCure::obs_log_likelihood() const
    {
        double obs_ell { 0 };
        arma::vec sus_prob { cure_obj.prob_vec };
        // for case 1
        for (size_t j: case1_ind) {
            obs_ell += std::log(sus_prob(j)) +
                std::log(cox_obj.h_time(j)) -
                cox_obj.H_time(j);
        }
        // for case 2
        for (size_t j: case2_ind) {
            obs_ell += std::log(
                sus_prob(j) * cox_obj.S_time(j) +
                (1 - sus_prob(j))
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
        const double pmin = 1e-5
        ) const
    {
        // check if the number of covariates matchs the fitted model
        if (new_cox_x.n_cols != this->cox_p) {
            throw std::range_error(
                "The number of columns ('new_cox_x') must match the model."
                );
        }
        if (new_cure_x.n_cols != this->cure_p0) {
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

        // add intercept if needed
        if (this->cure_p > this->cure_p0) {
            new_cure_x = arma::join_horiz(
                arma::ones(new_cure_x.n_rows), new_cure_x
                );
        }
        arma::uvec new_case1_ind { arma::find(new_event > 0) };
        arma::uvec new_case2_ind { arma::find(new_event < 1) };
        // construct the baseline survival curve
        // tail completion has already been applied to this->S0_est
        arma::vec S0_vec {
            arma::join_cols(arma::ones<arma::vec>(1), this->S0_est)
        };
        // baseline estimates
        arma::vec S_vec {
            step_fun(new_time, this->unique_time, S0_vec)
        };
        arma::vec H_vec { - arma::log(S_vec) };
        // only consider positive values
        arma::uvec which_h { arma::find(this->h0_est > 0) };
        arma::vec h_vec {
            step_fun2(new_time,
                      this->unique_time.elem(which_h),
                      this->h0_est.elem(which_h))
        };
        // apply x * beta
        // compute parts for the new data
        arma::vec new_cox_xbeta { mat2vec(new_cox_x * this->cox_coef) };
        arma::vec exp_cox_xbeta { arma::exp(new_cox_xbeta) };
        h_vec %= exp_cox_xbeta;
        H_vec %= exp_cox_xbeta;
        S_vec = arma::exp(- H_vec);
        arma::vec new_cure_xgamma { mat2vec(new_cure_x * this->cure_coef) };
        arma::vec p_vec { 1 / (1 + arma::exp(- new_cure_xgamma)) };
        for (size_t i {0}; i < p_vec.n_elem; ++i) {
            if (p_vec(i) < pmin) {
                p_vec(i) = pmin;
            } else if (p_vec(i) > 1 - pmin) {
                p_vec(i) = 1 - pmin;
            }
        }
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
