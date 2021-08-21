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

#ifndef CV_COXPH_CURE_MCAR_H
#define CV_COXPH_CURE_MCAR_H

#include <RcppArmadillo.h>
#include "cross-validation.h"
#include "coxph_cure_mcar.h"

namespace Intsurv {
    // estimation from cross-validation
    inline arma::vec cv_coxph_cure_mcar(
        const arma::vec& time,
        const arma::vec& event,
        const arma::mat& cox_x,
        const arma::mat& cure_x,
        const bool cure_intercept = true,
        const unsigned long n_folds = 10,
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const arma::vec& cox_offset = 0,
        const arma::vec& cure_offset = 0,
        const bool& cox_standardize = true,
        const bool& cure_standardize = true,
        const unsigned int& em_max_iter = 300,
        const double& em_rel_tol = 1e-5,
        const unsigned int& cox_mstep_max_iter = 100,
        const double& cox_mstep_rel_tol = 1e-5,
        const unsigned int& cure_mstep_max_iter = 100,
        const double& cure_mstep_rel_tol = 1e-5,
        const unsigned int& tail_completion = 1,
        double tail_tau = -1,
        const double& pmin = 1e-5,
        const unsigned int& early_stop = 0,
        const unsigned int& verbose = 0
        )
    {
        // stratify
        arma::vec event0na { event };
        const double const4na { 0.5 };
        event0na.replace(arma::datum::nan, const4na);

        const arma::uvec case1_ind { arma::find(event0na > const4na) };
        const arma::uvec case2_ind { arma::find(event0na < const4na) };
        const arma::uvec case3_ind { arma::find(event0na == const4na) };

        const arma::vec time_case1 { time.elem(case1_ind) };
        const arma::vec time_case2 { time.elem(case2_ind) };
        const arma::vec time_case3 { time.elem(case3_ind) };

        const arma::vec event_case1 { event0na.elem(case1_ind) };
        const arma::vec event_case2 { event0na.elem(case2_ind) };
        const arma::vec event_case3 { event0na.elem(case3_ind) };

        const arma::mat cox_x_case1 { cox_x.rows(case1_ind) };
        const arma::mat cox_x_case2 { cox_x.rows(case2_ind) };
        const arma::mat cox_x_case3 { cox_x.rows(case3_ind) };

        const arma::mat cure_x_case1 { cure_x.rows(case1_ind) };
        const arma::mat cure_x_case2 { cure_x.rows(case2_ind) };
        const arma::mat cure_x_case3 { cure_x.rows(case3_ind) };

        // process offset terms
        arma::vec cox_offset_case1 { arma::zeros(1) },
            cox_offset_case2 { cox_offset_case1 },
            cox_offset_case3 { cox_offset_case1 },
            cure_offset_case1 { cox_offset_case1 },
            cure_offset_case2 { cox_offset_case1 },
            cure_offset_case3 { cox_offset_case1 };
        if (cox_offset.n_elem == cox_x.n_rows) {
            cox_offset_case1 = cox_offset.elem(case1_ind);
            cox_offset_case2 = cox_offset.elem(case2_ind);
            cox_offset_case3 = cox_offset.elem(case3_ind);
        }
        if (cure_offset.n_elem == cure_x.n_rows) {
            cure_offset_case1 = cure_offset.elem(case1_ind);
            cure_offset_case2 = cure_offset.elem(case2_ind);
            cure_offset_case3 = cure_offset.elem(case3_ind);
        }

        // get the index of the largest event time
        const arma::uvec which_time_max { time_case1.index_max() };

        // cross-validation
        const unsigned long n_case1 { case1_ind.n_elem };
        const unsigned long n_case2 { case2_ind.n_elem };
        const unsigned long n_case3 { case3_ind.n_elem };

        CrossValidation cv_obj_case1 { n_case1, n_folds, which_time_max };
        CrossValidation cv_obj_case2 { n_case2, n_folds };
        CrossValidation cv_obj_case3 { n_case3, n_folds };

        arma::vec cv_vec { arma::zeros(n_folds) };

        for (size_t i {0}; i < n_folds; ++i) {
            // training set
            arma::vec train_time {
                arma::join_vert(
                    time_case1.elem(cv_obj_case1.train_index.at(i)),
                    time_case2.elem(cv_obj_case2.train_index.at(i)),
                    time_case3.elem(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::vec train_event {
                arma::join_vert(
                    event_case1.elem(cv_obj_case1.train_index.at(i)),
                    event_case2.elem(cv_obj_case2.train_index.at(i)),
                    event_case3.elem(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::mat train_cox_x {
                arma::join_vert(
                    cox_x_case1.rows(cv_obj_case1.train_index.at(i)),
                    cox_x_case2.rows(cv_obj_case2.train_index.at(i)),
                    cox_x_case3.rows(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::mat train_cure_x {
                arma::join_vert(
                    cure_x_case1.rows(cv_obj_case1.train_index.at(i)),
                    cure_x_case2.rows(cv_obj_case2.train_index.at(i)),
                    cure_x_case3.rows(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::vec train_cox_offset {
                arma::join_vert(
                    cox_offset_case1.elem(cv_obj_case1.train_index.at(i)),
                    cox_offset_case2.elem(cv_obj_case2.train_index.at(i)),
                    cox_offset_case3.elem(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::vec train_cure_offset {
                arma::join_vert(
                    cure_offset_case1.elem(cv_obj_case1.train_index.at(i)),
                    cure_offset_case2.elem(cv_obj_case2.train_index.at(i)),
                    cure_offset_case3.elem(cv_obj_case2.train_index.at(i))
                    )
            };
            // testing set
            arma::vec test_time {
                arma::join_vert(
                    time_case1.elem(cv_obj_case1.test_index.at(i)),
                    time_case2.elem(cv_obj_case2.test_index.at(i)),
                    time_case3.elem(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::vec test_event {
                arma::join_vert(
                    event_case1.elem(cv_obj_case1.test_index.at(i)),
                    event_case2.elem(cv_obj_case2.test_index.at(i)),
                    event_case3.elem(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::mat test_cox_x {
                arma::join_vert(
                    cox_x_case1.rows(cv_obj_case1.test_index.at(i)),
                    cox_x_case2.rows(cv_obj_case2.test_index.at(i)),
                    cox_x_case3.rows(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::mat test_cure_x {
                arma::join_vert(
                    cure_x_case1.rows(cv_obj_case1.test_index.at(i)),
                    cure_x_case2.rows(cv_obj_case2.test_index.at(i)),
                    cure_x_case3.rows(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::vec test_cox_offset {
                arma::join_vert(
                    cox_offset_case1.elem(cv_obj_case1.test_index.at(i)),
                    cox_offset_case2.elem(cv_obj_case2.test_index.at(i)),
                    cox_offset_case3.elem(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::vec test_cure_offset {
                arma::join_vert(
                    cure_offset_case1.elem(cv_obj_case1.test_index.at(i)),
                    cure_offset_case2.elem(cv_obj_case2.test_index.at(i)),
                    cure_offset_case3.elem(cv_obj_case2.test_index.at(i))
                    )
            };
            // define object
            CoxphCureMcar cc_obj {
                train_time, train_event, train_cox_x, train_cure_x,
                cure_intercept, cox_standardize, cure_standardize,
                train_cox_offset, train_cure_offset
            };
            // model-fitting
            cc_obj.fit(
                cox_start, cure_start,
                em_max_iter, em_rel_tol,
                cox_mstep_max_iter, cox_mstep_rel_tol,
                cure_mstep_max_iter, cure_mstep_rel_tol,
                tail_completion, tail_tau,
                pmin, early_stop, verbose
                );
            // compute observed log-likelihood function for the test data
            cv_vec(i) = cc_obj.obs_log_likelihood(
                test_time, test_event, test_cox_x, test_cure_x,
                test_cox_offset, test_cure_offset, pmin
                );
        }
        return cv_vec;
    }

    // reguaralized fit
    inline arma::vec cv_coxph_cure_mcar_reg(
        const arma::vec& time,
        const arma::vec& event,
        const arma::mat& cox_x,
        const arma::mat& cure_x,
        const bool cure_intercept = true,
        const unsigned long n_folds = 10,
        const double& cox_l1_lambda = 0,
        const double& cox_l2_lambda = 0,
        const arma::vec& cox_l1_penalty_factor = 0,
        const double& cure_l1_lambda = 0,
        const double& cure_l2_lambda = 0,
        const arma::vec& cure_l1_penalty_factor = 0,
        const arma::vec& cox_start = 0,
        const arma::vec& cure_start = 0,
        const arma::vec& cox_offset = 0,
        const arma::vec& cure_offset = 0,
        const bool cox_standardize = true,
        const bool cure_standardize = true,
        const unsigned int& em_max_iter = 500,
        const double& em_rel_tol = 1e-4,
        const unsigned int& cox_mstep_max_iter = 200,
        const double& cox_mstep_rel_tol = 1e-4,
        const unsigned int& cure_mstep_max_iter = 200,
        const double& cure_mstep_rel_tol = 1e-4,
        const unsigned int& tail_completion = 1,
        double tail_tau = -1,
        const double& pmin = 1e-5,
        const unsigned int& early_stop = 0,
        const unsigned int& verbose = 0
        )
    {
        // stratify
        arma::vec event0na { event };
        const double const4na { 0.5 };
        event0na.replace(arma::datum::nan, const4na);

        const arma::uvec case1_ind { arma::find(event0na > const4na) };
        const arma::uvec case2_ind { arma::find(event0na < const4na) };
        const arma::uvec case3_ind { arma::find(event0na == const4na) };

        const arma::vec time_case1 { time.elem(case1_ind) };
        const arma::vec time_case2 { time.elem(case2_ind) };
        const arma::vec time_case3 { time.elem(case3_ind) };

        const arma::vec event_case1 { event0na.elem(case1_ind) };
        const arma::vec event_case2 { event0na.elem(case2_ind) };
        const arma::vec event_case3 { event0na.elem(case3_ind) };

        const arma::mat cox_x_case1 { cox_x.rows(case1_ind) };
        const arma::mat cox_x_case2 { cox_x.rows(case2_ind) };
        const arma::mat cox_x_case3 { cox_x.rows(case3_ind) };

        const arma::mat cure_x_case1 { cure_x.rows(case1_ind) };
        const arma::mat cure_x_case2 { cure_x.rows(case2_ind) };
        const arma::mat cure_x_case3 { cure_x.rows(case3_ind) };

        // process offset terms
        arma::vec cox_offset_case1 { arma::zeros(1) },
            cox_offset_case2 { cox_offset_case1 },
            cox_offset_case3 { cox_offset_case1 },
            cure_offset_case1 { cox_offset_case1 },
            cure_offset_case2 { cox_offset_case1 },
            cure_offset_case3 { cox_offset_case1 };
        if (cox_offset.n_elem == cox_x.n_rows) {
            cox_offset_case1 = cox_offset.elem(case1_ind);
            cox_offset_case2 = cox_offset.elem(case2_ind);
            cox_offset_case3 = cox_offset.elem(case3_ind);
        }
        if (cure_offset.n_elem == cure_x.n_rows) {
            cure_offset_case1 = cure_offset.elem(case1_ind);
            cure_offset_case2 = cure_offset.elem(case2_ind);
            cure_offset_case3 = cure_offset.elem(case3_ind);
        }

        // get the index of the largest event time
        const arma::uvec which_time_max { time_case1.index_max() };

        // cross-validation
        const unsigned long n_case1 { case1_ind.n_elem };
        const unsigned long n_case2 { case2_ind.n_elem };
        const unsigned long n_case3 { case3_ind.n_elem };

        CrossValidation cv_obj_case1 { n_case1, n_folds, which_time_max };
        CrossValidation cv_obj_case2 { n_case2, n_folds };
        CrossValidation cv_obj_case3 { n_case3, n_folds };

        arma::vec cv_vec { arma::zeros(n_folds) };

        for (size_t i {0}; i < n_folds; ++i) {
            // training set
            arma::vec train_time {
                arma::join_vert(
                    time_case1.elem(cv_obj_case1.train_index.at(i)),
                    time_case2.elem(cv_obj_case2.train_index.at(i)),
                    time_case3.elem(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::vec train_event {
                arma::join_vert(
                    event_case1.elem(cv_obj_case1.train_index.at(i)),
                    event_case2.elem(cv_obj_case2.train_index.at(i)),
                    event_case3.elem(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::mat train_cox_x {
                arma::join_vert(
                    cox_x_case1.rows(cv_obj_case1.train_index.at(i)),
                    cox_x_case2.rows(cv_obj_case2.train_index.at(i)),
                    cox_x_case3.rows(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::mat train_cure_x {
                arma::join_vert(
                    cure_x_case1.rows(cv_obj_case1.train_index.at(i)),
                    cure_x_case2.rows(cv_obj_case2.train_index.at(i)),
                    cure_x_case3.rows(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::vec train_cox_offset {
                arma::join_vert(
                    cox_offset_case1.elem(cv_obj_case1.train_index.at(i)),
                    cox_offset_case2.elem(cv_obj_case2.train_index.at(i)),
                    cox_offset_case3.elem(cv_obj_case3.train_index.at(i))
                    )
            };
            arma::vec train_cure_offset {
                arma::join_vert(
                    cure_offset_case1.elem(cv_obj_case1.train_index.at(i)),
                    cure_offset_case2.elem(cv_obj_case2.train_index.at(i)),
                    cure_offset_case3.elem(cv_obj_case2.train_index.at(i))
                    )
            };
            // testing set
            arma::vec test_time {
                arma::join_vert(
                    time_case1.elem(cv_obj_case1.test_index.at(i)),
                    time_case2.elem(cv_obj_case2.test_index.at(i)),
                    time_case3.elem(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::vec test_event {
                arma::join_vert(
                    event_case1.elem(cv_obj_case1.test_index.at(i)),
                    event_case2.elem(cv_obj_case2.test_index.at(i)),
                    event_case3.elem(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::mat test_cox_x {
                arma::join_vert(
                    cox_x_case1.rows(cv_obj_case1.test_index.at(i)),
                    cox_x_case2.rows(cv_obj_case2.test_index.at(i)),
                    cox_x_case3.rows(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::mat test_cure_x {
                arma::join_vert(
                    cure_x_case1.rows(cv_obj_case1.test_index.at(i)),
                    cure_x_case2.rows(cv_obj_case2.test_index.at(i)),
                    cure_x_case3.rows(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::vec test_cox_offset {
                arma::join_vert(
                    cox_offset_case1.elem(cv_obj_case1.test_index.at(i)),
                    cox_offset_case2.elem(cv_obj_case2.test_index.at(i)),
                    cox_offset_case3.elem(cv_obj_case3.test_index.at(i))
                    )
            };
            arma::vec test_cure_offset {
                arma::join_vert(
                    cure_offset_case1.elem(cv_obj_case1.test_index.at(i)),
                    cure_offset_case2.elem(cv_obj_case2.test_index.at(i)),
                    cure_offset_case3.elem(cv_obj_case2.test_index.at(i))
                    )
            };
            // define object
            CoxphCureMcar cc_obj {
                train_time, train_event, train_cox_x, train_cure_x,
                cure_intercept, cox_standardize, cure_standardize,
                train_cox_offset, train_cure_offset
            };
            // model-fitting
            cc_obj.regularized_fit(
                cox_l1_lambda, cox_l2_lambda,
                cure_l1_lambda, cure_l2_lambda,
                cox_l1_penalty_factor, cure_l1_penalty_factor,
                cox_start, cure_start,
                em_max_iter, em_rel_tol,
                cox_mstep_max_iter, cox_mstep_rel_tol,
                cure_mstep_max_iter, cure_mstep_rel_tol,
                tail_completion, tail_tau,
                pmin, early_stop, verbose
                );
            // compute observed log-likelihood function for the test data
            cv_vec(i) = cc_obj.obs_log_likelihood(
                test_time, test_event, test_cox_x, test_cure_x,
                test_cox_offset, test_cure_offset, pmin
                );
        }
        return cv_vec;
    }


}  // Intsurv


#endif /* CV_COXPH_CURE_MCAR_H */
