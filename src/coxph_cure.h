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

#ifndef COXPH_CURE_H
#define COXPH_CURE_H

#include <RcppArmadillo.h>
#include <string>
#include "coxph_reg.h"
#include "logistic_reg.h"
#include "utils.h"


namespace Intsurv {

    class CoxphCure {
    private:
        CoxphReg cox_obj;
        LogisticReg cure_obj;
        unsigned int cox_p;
        unsigned int cure_p;
        arma::uvec case1_ind;
        arma::uvec case2_ind;

    public:
        arma::vec cox_coef;
        arma::vec cure_coef;
        double negLogL;
        unsigned int nObs;        // number of observations
        unsigned int coef_df;     // degree of freedom of coef estimates
        unsigned int num_iter;    // number of iterations

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
            this->cure_p = cure_x.n_cols +
                static_cast<unsigned int>(cure_intercept);
            this->nObs = cox_x.n_rows;
            arma::uvec cox_sort_ind { cox_obj.get_sort_index() };
            arma::mat cure_xx { cure_x.rows(cox_sort_ind) };
            arma::vec s_event { event.elem(cox_sort_ind) };
            this->case1_ind = arma::find(s_event > 0);
            this->case2_ind = arma::find(s_event < 1);
            // create the LogisticReg object
            this->cure_obj = LogisticReg(cure_xx, s_event, cure_intercept,
                                         cure_standardize);
        }

        // function members
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
            const bool verbose
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
            const bool verbose
            );

    };                          // end of class definition


    // fit the Cox cure mode by EM algorithm
    inline void CoxphCure::fit(
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const unsigned int& em_max_iter = 500,
        const double& em_rel_tol = 1e-5,
        const unsigned int& cox_mstep_max_iter = 100,
        const double& cox_mstep_rel_tol = 1e-3,
        const unsigned int& cure_mstep_max_iter = 100,
        const double& cure_mstep_rel_tol = 1e-3,
        const bool verbose = false
        )
    {
        // initialize cox_beta
        arma::vec cox_beta { arma::zeros(this->cox_p) };
        if (cox_start.n_elem == this->cox_p) {
            cox_beta = cox_start;
        } else {
            arma::mat tmp_cox_x { cox_obj.get_x() };
            arma::vec tmp_time { cox_obj.get_time() };
            arma::vec tmp_event { cox_obj.get_event() };
            arma::uvec tmp_idx { arma::find(tmp_event > 0) };
            CoxphReg tmp_object {
                CoxphReg(tmp_time.elem(tmp_idx),
                         tmp_event.elem(tmp_idx),
                         tmp_cox_x.rows(tmp_idx))
            };
            tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
            cox_beta = tmp_object.coef;
        }
        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(this->cure_p) };
        if (cure_start.n_elem == this->cure_p) {
            cure_beta = cure_start;
        } else {
            cure_obj.fit(cure_beta, cure_mstep_max_iter,
                         cure_mstep_rel_tol);
            cure_beta = cure_obj.coef;
        }
        cure_obj.coef = cure_beta;
        arma::vec p_vec { cure_obj.predict(cure_beta) };
        cox_obj.coef = cox_beta;
        cox_obj.compute_haz_surv_time(cox_beta);
        size_t i {0};
        double obs_ell {0}, obs_ell_old { - arma::datum::inf };

        // allow users to stop the main loop
        Rcpp::checkUserInterrupt();

        // if verbose, print starting value
        if (verbose) {
            Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                        << "starting values:\n"
                        << "  Cox coef: "
                        << arma2rvec(cox_beta) << "\n"
                        << "  cure coef: "
                        << arma2rvec(cure_beta)
                        << std::endl;
        }

        // main loop of EM algorithm
        while (true) {
            // update iter
            ++i;

            arma::vec estep_v { cox_obj.get_event() };
            double numer_j {0}, tol1 {0}, tol2 {0};
            // reset values
            obs_ell = 0;

            // E-step: compute v vector
            for (size_t j: case2_ind) {
                numer_j = p_vec(j) * cox_obj.S_time(j);
                // more numerically robust
                estep_v(j) = 1 / ((1 - p_vec(j)) / numer_j + 1);
            }
            arma::vec log_v { arma::log(estep_v) };

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            cox_obj.set_offset(log_v);
            cox_obj.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the Cure layer
            cure_obj.update_y(estep_v);
            cure_obj.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);

            // check convergence
            tol1 = rel_l1_norm(cox_obj.coef, cox_beta);
            tol2 = rel_l1_norm(cure_obj.coef, cure_beta);

            // prepare the E-step for next iteration
            cox_obj.compute_haz_surv_time();
            p_vec = cure_obj.predict(cure_beta);

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // compute observed log-likelihood
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

            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                            << "iteration: " << i << "\n"
                            << "  Cox coef: "
                            << arma2rvec(cox_beta) << "\n"
                            << "    relative diff: " << tol1 << "\n"
                            << "  cure coef: "
                            << arma2rvec(cure_beta) << "\n"
                            << "    relative diff: " << tol2
                            << std::endl;
                Rcpp::Rcout << "  observed negative log-likelihood: "
                            << - obs_ell << std::endl;
            }

            // early exit if has any `nan`
            if (cox_beta.has_nan() || cure_beta.has_nan()) {
                obs_ell = - arma::datum::inf;
                throw std::range_error(
                    "The negative log-likelihood function went to infinite."
                    );
                break;
            }

            // early exit if negative log-likelihood increased
            if (obs_ell < obs_ell_old) {
                Rcpp::Rcout << "\nWarning: "
                            << "The negative log-likelihood increased.\n"
                            << "Ended the EM algorithm after "
                            << i
                            << " iterations with estimates from the last step."
                            << std::endl;
                // take the estimates from the last step
                cox_obj.coef = cox_beta;
                cure_obj.coef = cure_beta;
                // compute hazard and survival function estimates
                cox_obj.compute_haz_surv_time();
                cox_obj.est_haz_surv();
                // break here
                break;
            }

            // check convergence
            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i >= em_max_iter) {
                if (i < em_max_iter) {
                    Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                                << "reached convergence after " << i
                                << " iterations\n" << std::endl;
                } else {
                    Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                                << "reached the max iteration number."
                                << std::endl;
                }
                // compute hazard and survival function estimates
                cox_obj.est_haz_surv();
                // get out of the loop here
                break;
            }

            // update to last estimates
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;

        } // end of the EM algorithm
        // prepare outputs
        this->cox_coef = cox_obj.coef;
        this->cure_coef = cure_obj.coef;
        this->negLogL = - obs_ell;
        this->coef_df = cox_obj.coef_df + cure_obj.coef_df;
        this->num_iter = i;
        this->unique_time = cox_obj.unique_time;
        this->h0_est = cox_obj.h0_est;
        this->H0_est = cox_obj.H0_est;
        this->S0_est = cox_obj.S0_est;
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
        const unsigned int& cox_mstep_max_iter = 50,
        const double& cox_mstep_rel_tol = 1e-3,
        const unsigned int& cure_mstep_max_iter = 50,
        const double& cure_mstep_rel_tol = 1e-3,
        const bool verbose = false
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
        // L1 penalty factor for Cure model
        arma::vec cure_l1_penalty { arma::ones(this->cure_p) };
        if (cure_l1_penalty_factor.n_elem == cure_p) {
            // re-scale so that sum(factor) = number of predictors
            cure_l1_penalty = cure_l1_penalty_factor * cure_p /
                arma::sum(cure_l1_penalty_factor);
        }
        this->cure_l1_penalty_factor = cure_l1_penalty;

        // initialized with all zeros coef
        arma::vec cox_beta { arma::zeros(cox_p) };
        arma::vec cure_beta { arma::zeros(cure_p) };

        // compute the large enough lambdas that result in all-zero estimates
        arma::vec cox_grad_zero { arma::abs(cox_obj.gradient(cox_beta)) };
        this->cox_l1_lambda_max =
            arma::max(cox_grad_zero / cox_l1_penalty) / this->nObs;
        arma::vec cure_grad_zero { arma::abs(cure_obj.gradient(cure_beta)) };
        this->cure_l1_lambda_max =
            arma::max(cure_grad_zero.tail(cure_l1_penalty.n_elem) /
                      cure_l1_penalty) / this->nObs;

        // set the start estimates
        if (cox_start.n_elem == cox_p) {
            cox_beta = cox_start;
        }
        if (cure_start.n_elem == cure_p) {
            cure_beta = cure_start;
        }
        cure_obj.coef = cure_beta;
        cox_obj.coef = cox_beta;

        // set iterators
        arma::vec p_vec { arma::zeros(this->nObs) };
        size_t i {0};
        double obs_ell {0};
        double reg_obj {0}, reg_obj_old { arma::datum::inf };

        // allow users to stop here
        Rcpp::checkUserInterrupt();

        // main loop of EM algorithm
        while (true) {
            // update to last estimates
            p_vec = cure_obj.predict(cure_beta);
            cox_obj.compute_haz_surv_time();
            obs_ell = 0;

            // compute observed log-likelihood
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
                cox_l1_lambda *
                    arma::sum(arma::abs(cox_obj.coef) % cox_l1_penalty) +
                cox_l2_lambda * arma::as_scalar(crossprod(cox_obj.coef))
            };
            double reg_cure {
                cure_l1_lambda *
                arma::sum(arma::abs(cure_obj.coef.tail_rows(cure_p)) %
                          cure_l1_penalty) +
                cure_l2_lambda *
                arma::as_scalar(crossprod(cure_obj.coef.tail_rows(cure_p)))
            };
            reg_obj = - obs_ell / this->nObs + reg_cox + reg_cure;

            // verbose tracing for objective function
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                            << "iteration: " << i << "\n"
                            << "  Cox coef: "
                            << arma2rvec(cox_obj.coef) << "\n"
                            << "  cure coef: "
                            << arma2rvec(cure_obj.coef) << "\n"
                            << std::endl;
                Rcpp::Rcout << "  observed negative log-likelihood: "
                            << - obs_ell << std::endl;
                Rcpp::Rcout << "  regularized objective function: "
                            << reg_obj << std::endl;
                Rcpp::Rcout << "    penalty on Cox model: "
                            << reg_cox << std::endl;
                Rcpp::Rcout << "    penalty on cure layer: "
                            << reg_cure << std::endl;
            }

            // prepare for the E-step
            arma::vec estep_v { cox_obj.get_event() };
            double numer_j {0}, denom_j {0};

            // E-step: compute v vector
            for (size_t j: case2_ind) {
                numer_j = p_vec(j) * cox_obj.S_time(j);
                denom_j = numer_j + 1 - p_vec(j);
                estep_v(j) = numer_j / denom_j;
            }
            arma::vec log_v { arma::log(estep_v) };

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the survival layer
            cox_obj.set_offset(log_v);
            cox_obj.regularized_fit(
                cox_l1_lambda, cox_l2_lambda, cox_l1_penalty_factor,
                cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol
                );

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // M-step for the Cure layer
            cure_obj.update_y(estep_v);
            cure_obj.regularized_fit(
                cure_l1_lambda, cure_l2_lambda, cure_l1_penalty_factor,
                cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol
                );

            bool early_exit { false };
            // early exit if has any `nan`
            if (cox_obj.coef.has_nan() || cure_obj.coef.has_nan()) {
                obs_ell = - arma::datum::inf;
                Rcpp::Rcout << "Warning: "
                            << "Found NA's. "
                            << "The objective function went to infinite."
                            << std::endl;
                Rcpp::Rcout << "Ended the EM algorithm after "
                            << i
                            << " iterations with estimates from the last step."
                            << std::endl;
                early_exit = true;
            }
            // early exit if the regularized objective function increased, which
            // is technically impossible but may suggest a convergence
            // numerically
            if (reg_obj > reg_obj_old) {
                if (verbose) {
                    Rcpp::Rcout
                        << "\nThe regularized objective function increased.\n"
                        << "Ended the EM algorithm after "
                        << i
                        << " iterations with estimates from the last step."
                        << std::endl;
                }
                early_exit = true;
            }
            // return the estimates from last step
            if (early_exit) {
                // compute hazard and survival function estimates
                cox_obj.coef = cox_beta;
                cure_obj.coef = cure_beta;
                cox_obj.compute_censor_haz_surv_time();
                cox_obj.est_haz_surv();
                // break here
                break;
            }

            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // check convergence
            double tol1 { rel_l1_norm(cox_obj.coef, cox_beta) };
            double tol2 { rel_l1_norm(cure_obj.coef, cure_beta) };

            if (verbose) {
                Rcpp::Rcout << "\n  relative tolerance of\n"
                            << "    Cox coef: " << tol1 << "\n"
                            << "    cure coef: " << tol2 << std::endl;
            }

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i > em_max_iter) {
                if (verbose) {
                    if (i < em_max_iter) {
                        Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                                    << "reached convergence after " << i
                                    << " iterations\n" << std::endl;
                    } else {
                        Rcpp::Rcout << "\n" << std::string(40, '=') << "\n"
                                    << "reached the max iteration number."
                                    << std::endl;
                    }
                }
                // compute hazard and survival function estimates
                cox_obj.est_haz_surv();
                // get out of the loop here
                break;
            }
            // update fitted objects
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;
            reg_obj_old = reg_obj;

            // update iter
            ++i;

        } // end of the EM algorithm
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
    }

}

#endif
