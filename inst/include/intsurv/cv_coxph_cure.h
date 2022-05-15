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

#ifndef INTSURV_CV_COXPH_CURE_H
#define INTSURV_CV_COXPH_CURE_H

#include <RcppArmadillo.h>
#include "cross-validation.h"
#include "coxph_cure.h"
#include "control.h"
#include "subset.h"

namespace intsurv {

    // estimation from cross-validation
    inline arma::vec cv_coxph_cure(
        const CoxphCure& object,
        const unsigned int nfolds = 5
        )
    {
        // stratify
        const arma::uvec& case1_ind { object.case1_ind_ };
        const arma::uvec& case2_ind { object.case2_ind_ };
        // FIXME: if there are ties among the max event time,
        // a better way is to form a vector;
        // currently, the first max event time stays in the training set
        const arma::uvec which_time_max { object.max_event_time_ind_ };
        // cross-validation
        const unsigned int n_case1 { case1_ind.n_elem };
        const unsigned int n_case2 { case2_ind.n_elem };
        CrossValidation cv_obj_case1 { n_case1, nfolds, which_time_max };
        CrossValidation cv_obj_case2 { n_case2, nfolds };
        arma::vec cv_vec { arma::zeros(nfolds) };
        for (size_t i {0}; i < nfolds; ++i) {
            arma::uvec train_idx {
                arma::join_vert(cv_obj_case1.train_index_.at(i),
                                cv_obj_case2.train_index_.at(i))
            };
            arma::uvec valid_idx {
                arma::join_vert(cv_obj_case1.test_index_.at(i),
                                cv_obj_case2.test_index_.at(i))
            };
            CoxphCure train_obj { subset(object, train_idx) };
            CoxphCure valid_obj { subset(object, valid_idx) };
            // set verbose to zero
            train_obj.control_.set_verbose(0);
            train_obj.surv_obj_.control_.set_verbose(0);
            train_obj.cure_obj_.control_.set_verbose(0);
            // model-fitting
            train_obj.fit();
            // compute observed log-likelihood function for the test data
            cv_vec(i) = train_obj.obs_log_likelihood(valid_obj);
        }
        return cv_vec;
    }

    // regularized Cox cure model
    inline arma::vec cv_coxph_cure_reg(
        const CoxphCure& object,
        const unsigned int nfolds = 5
        )
    {
        // stratify
        const arma::uvec& case1_ind { object.case1_ind_ };
        const arma::uvec& case2_ind { object.case2_ind_ };
        const arma::uvec which_time_max { object.max_event_time_ind_ };
        // cross-validation
        const unsigned int n_case1 { case1_ind.n_elem };
        const unsigned int n_case2 { case2_ind.n_elem };
        CrossValidation cv_obj_case1 { n_case1, nfolds, which_time_max };
        CrossValidation cv_obj_case2 { n_case2, nfolds };
        arma::vec cv_vec { arma::zeros(nfolds) };
        for (size_t i {0}; i < nfolds; ++i) {
            // training set
            arma::uvec train_idx {
                arma::join_vert(cv_obj_case1.train_index_.at(i),
                                cv_obj_case2.train_index_.at(i))
            };
            arma::uvec valid_idx {
                arma::join_vert(cv_obj_case1.test_index_.at(i),
                                cv_obj_case2.test_index_.at(i))
            };
            CoxphCure train_obj { subset(object, train_idx) };
            CoxphCure valid_obj { subset(object, valid_idx) };
            // set verbose to zero
            train_obj.control_.set_verbose(0);
            train_obj.surv_obj_.control_.set_verbose(0);
            train_obj.cure_obj_.control_.set_verbose(0);
            // model-fitting
            train_obj.net_fit();
            // compute observed log-likelihood function for the test data
            cv_vec(i) = train_obj.obs_log_likelihood(valid_obj);
        }
        return cv_vec;
    }

}  // intsurv


#endif
