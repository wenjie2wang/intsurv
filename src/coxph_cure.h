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
        unsigned int nObs;        // number of observations
        // the "big enough" L1 lambda => zero coef
        double cox_l1_lambda_max;
        double cure_l1_lambda_max;

        arma::vec cox_coef;
        arma::vec cure_coef;
        double negLogL;

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
        CoxphCure() {}

        // constructors
        CoxphCure(const arma::vec& time,
                  const arma::vec& event,
                  const arma::mat& cox_x,
                  const arma::mat& cure_x,
                  const bool cure_intercept = true,
                  const bool cox_standardize = true,
                  const bool cure_standardize = true)
        {
            // create the CoxphReg object
            this->cox_obj = CoxphReg(time, event, cox_x, cox_standardize);
            // pre-process x and y
            this->cox_p = cox_x.n_cols;
            this->cure_p = cure_x.n_cols +
                static_cast<unsigned int>(cure_intercept);
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
            const double& cure_mstep_rel_tol
            );

        // fit regularized Cox cure model with adaptive elastic net penalty
        // for a perticular lambda
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
            const double& cure_mstep_rel_tol
            );

    };
    // end of class definition


    // fit the Cox cure mode by EM algorithm
    inline void CoxphCure::fit(
        const arma::vec& cox_start,
        const arma::vec& cure_start,
        const unsigned int& em_max_iter,
        const double& em_rel_tol,
        const unsigned int& cox_mstep_max_iter,
        const double& cox_mstep_rel_tol,
        const unsigned int& cure_mstep_max_iter,
        const double& cure_mstep_rel_tol
        )
    {
        // initialize cox_beta
        arma::vec cox_beta { arma::zeros(cox_p) };
        if (cox_start.n_elem == cox_p) {
            cox_beta = cox_start;
        } else {
            arma::mat tmp_cox_x { cox_obj.get_x() };
            arma::vec tmp_time { cox_obj.get_time() };
            arma::vec tmp_event { cox_obj.get_event() };
            arma::uvec tmp_idx { arma::find(tmp_event > 0) };
            Intsurv::CoxphReg tmp_object {
                Intsurv::CoxphReg(tmp_time.elem(tmp_idx),
                                  tmp_event.elem(tmp_idx),
                                  tmp_cox_x.rows(tmp_idx))
            };
            tmp_object.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);
            cox_beta = tmp_object.coef;
        }
        // initialize cure_beta
        arma::vec cure_beta { arma::zeros(cure_p) };
        if (cure_start.n_elem == cure_p) {
            cure_beta = cure_start;
        } else {
            cure_obj.fit(cure_beta, cure_mstep_max_iter,
                         cure_mstep_rel_tol);
            cure_beta = cure_obj.coef;
        }
        arma::vec p_vec { cure_obj.predict(cure_beta) };
        cox_obj.compute_haz_surv_time(cox_beta);
        size_t i {0};
        double obs_ell {0};

        // main loop of EM algorithm
        while (true) {
            arma::vec estep_v { cox_obj.get_event() };
            double numer_j {0}, denom_j {0}, tol1 {0}, tol2 {0};

            // E-step: compute v vector
            for (size_t j: case2_ind) {
                numer_j = p_vec(j) * cox_obj.S_time(j);
                denom_j = numer_j + 1 - p_vec(j);
                estep_v(j) = numer_j / denom_j;
            }
            arma::vec log_v { arma::log(estep_v) };

            // M-step for the survival layer
            cox_obj.set_offset(log_v);
            cox_obj.fit(cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol);

            // M-step for the Cure layer
            cure_obj.update_y(estep_v);
            cure_obj.fit(cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol);

            // check convergence
            tol1 = Intsurv::rel_l2_norm(cox_obj.coef, cox_beta);
            tol2 = Intsurv::rel_l2_norm(cure_obj.coef, cure_beta);
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;

            // update to last estimates
            cox_obj.compute_haz_surv_time();
            p_vec = cure_obj.predict(cure_beta);

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i > em_max_iter) {
                // compute observed data likelihood
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
                break;
            }
            // update iter
            ++i;
        } // end of the EM algorithm
        // prepare outputs
        this->cox_coef = cox_beta;
        this->cure_coef = cure_beta;
        this->negLogL = - obs_ell;
    }

    // fit regularized Cox cure model with adaptive elastic net penalty
    // for a perticular lambda
    // lambda_1 * lasso * factors + lambda_2 * ridge
    inline void CoxphCure::regularized_fit(
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
        const double& cure_mstep_rel_tol
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
            arma::max(cox_grad_zero / cox_l1_penalty) / cox_obj.nObs;
        arma::vec cure_grad_zero { arma::abs(cure_obj.gradient(cure_beta)) };
        this->cure_l1_lambda_max =
            arma::max(cure_grad_zero.tail(cure_l1_penalty.n_elem) /
                      cure_l1_penalty) / cure_obj.nObs;

        // set the start estimates
        if (cox_start.n_elem == cox_p) {
            cox_beta = cox_start;
        }
        if (cure_start.n_elem == cure_p) {
            cure_beta = cure_start;
        }

        // set iterator and prepare the E-step
        size_t i {0};
        double obs_ell {0};
        arma::vec p_vec { cure_obj.predict(cure_beta) };
        cox_obj.compute_haz_surv_time(cox_beta);

        // main loop of EM algorithm
        while (true) {
            arma::vec estep_v { cox_obj.get_event() };
            double numer_j {0}, denom_j {0}, tol1 {0}, tol2 {0};

            // E-step: compute v vector
            for (size_t j: case2_ind) {
                numer_j = p_vec(j) * cox_obj.S_time(j);
                denom_j = numer_j + 1 - p_vec(j);
                estep_v(j) = numer_j / denom_j;
            }
            arma::vec log_v { arma::log(estep_v) };

            // M-step for the survival layer
            cox_obj.set_offset(log_v);
            cox_obj.regularized_fit(
                cox_l1_lambda, cox_l2_lambda, cox_l1_penalty_factor,
                cox_beta, cox_mstep_max_iter, cox_mstep_rel_tol
                );

            // M-step for the Cure layer
            cure_obj.update_y(estep_v);
            cure_obj.regularized_fit(
                cure_l1_lambda, cure_l2_lambda, cure_l1_penalty_factor,
                cure_beta, cure_mstep_max_iter, cure_mstep_rel_tol
                );

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

            // check convergence
            tol1 = Intsurv::rel_l2_norm(cox_obj.coef, cox_beta);
            tol2 = Intsurv::rel_l2_norm(cure_obj.coef, cure_beta);
            cox_beta = cox_obj.coef;
            cure_beta = cure_obj.coef;

            // update to last estimates
            cox_obj.compute_haz_surv_time();
            p_vec = cure_obj.predict(cure_beta);

            if ((tol1 < em_rel_tol && tol2 < em_rel_tol) || i > em_max_iter) {
                // compute observed data likelihood
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
                break;
            }
            // update iter
            ++i;
        } // end of the EM algorithm
        // prepare outputs
        this->cox_coef = cox_beta;
        this->cure_coef = cure_beta;
        this->en_cox_coef = cox_obj.en_coef;
        this->en_cure_coef = cure_obj.en_coef;
        this->negLogL = - obs_ell;
        this->cox_l1_lambda = cox_l1_lambda;
        this->cox_l2_lambda = cox_l2_lambda;
        this->cure_l1_lambda = cure_l1_lambda;
        this->cure_l2_lambda = cure_l2_lambda;
    }

}

#endif
