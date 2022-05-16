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
        const arma::uvec which_time_max { object.max_event_time_ind_ };
        // stratified cross-validation
        arma::uvec cv_strata { arma::zeros<arma::uvec>(object.n_obs_) };
        cv_strata.elem(object.case2_ind_).ones();
        CrossValidation cv_obj { object.n_obs_, nfolds,
            which_time_max, cv_strata };
        arma::vec cv_vec { arma::zeros(nfolds) };
        for (size_t i {0}; i < nfolds; ++i) {
            CoxphCure train_obj {
                subset(object, cv_obj.train_index_.at(i)) };
            CoxphCure valid_obj {
                subset(object, cv_obj.test_index_.at(i)) };
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
        const arma::uvec which_time_max { object.max_event_time_ind_ };
        // stratified cross-validation
        arma::uvec cv_strata { arma::zeros<arma::uvec>(object.n_obs_) };
        cv_strata.elem(object.case2_ind_).ones();
        CrossValidation cv_obj { object.n_obs_, nfolds,
            which_time_max, cv_strata };
        arma::vec cv_vec { arma::zeros(nfolds) };
        for (size_t i {0}; i < nfolds; ++i) {
            CoxphCure train_obj {
                subset(object, cv_obj.train_index_.at(i)) };
            CoxphCure valid_obj {
                subset(object, cv_obj.test_index_.at(i)) };
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
