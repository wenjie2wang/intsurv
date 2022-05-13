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

#ifndef INTSURV_CROSS_VALIDATION_H
#define INTSURV_CROSS_VALIDATION_H

#include <vector>
#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {
    class CrossValidation {
    protected:
        unsigned int n_obs_;
        unsigned int n_folds_ = 5;

    public:
        std::vector<arma::uvec> train_index_;
        std::vector<arma::uvec> test_index_;

        // default constructor
        CrossValidation();

        // major constructor
        CrossValidation(const unsigned int n_obs,
                        const unsigned int n_folds) :
            n_obs_ { n_obs },
            n_folds_ { n_folds }
        {
            test_index_ = get_cv_test_index(n_obs_, n_folds_);
            arma::uvec all_index {
                arma::regspace<arma::uvec>(0, n_obs - 1)
            };
            for (size_t i {0}; i < n_folds_; ++i) {
                train_index_.push_back(
                    vec_diff(all_index, test_index_.at(i))
                    );
            }
        }

        // explicit constructor
        explicit CrossValidation(const unsigned int n_obs) :
            n_obs_ { n_obs }
        {
            CrossValidation(n_obs_, 5);
        }

        // a special constructor
        CrossValidation(const unsigned int n_obs,
                        const unsigned int n_folds,
                        // keep static index always in the training set
                        const arma::uvec static_train_index) :
            n_obs_ { n_obs },
            n_folds_ { n_folds }
        {
            test_index_ = get_cv_test_index(n_obs_, n_folds_);
            // remove static train index from test index
            for (size_t i {0}; i < n_folds_; ++i) {
                test_index_.at(i) = vec_diff(
                    test_index_.at(i), static_train_index
                    );
            }
            arma::uvec all_index {
                arma::regspace<arma::uvec>(0, n_obs - 1)
            };
            for (size_t i {0}; i < n_folds_; ++i) {
                train_index_.push_back(
                    vec_diff(all_index, test_index_.at(i))
                    );
            }
        }

        // helper function
        unsigned int get_n_folds() const
        {
            return n_folds_;
        }
        unsigned int get_n_obs() const
        {
            return n_obs_;
        }

    };
}

#endif
