//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2020  Wenjie Wang <wjwang.stat@gmail.com>
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

#ifndef CROSS_VALIDATION_H
#define CROSS_VALIDATION_H

#include <vector>
#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {
    class CrossValidation {
    private:
        unsigned long n_obs_;
        unsigned long n_folds_ = 10;
    public:
        std::vector<arma::uvec> train_index;
        std::vector<arma::uvec> test_index;

        // default constructor
        CrossValidation();

        // major constructor
        CrossValidation(const unsigned long n_obs,
                        const unsigned long n_folds) :
            n_obs_ { n_obs },
            n_folds_ { n_folds }
        {
            this->test_index = get_cv_test_index(n_obs_, n_folds_);
            arma::uvec all_index {
                arma::regspace<arma::uvec>(0, n_obs - 1)
            };
            for (size_t i {0}; i < n_folds_; ++i) {
                this->train_index.push_back(
                    vec_diff(all_index, this->test_index.at(i))
                    );
            }
        }

        // explicit constructor
        explicit CrossValidation(const unsigned long n_obs) :
            n_obs_ { n_obs }
        {
            CrossValidation(n_obs_, 10);
        }

        // a special constructor
        CrossValidation(const unsigned long n_obs,
                        const unsigned long n_folds,
                        // keep static index always in the training set
                        const arma::uvec static_train_index) :
            n_obs_ { n_obs },
            n_folds_ { n_folds }
        {
            this->test_index = get_cv_test_index(n_obs_, n_folds_);
            // remove static train index from test index
            for (size_t i {0}; i < n_folds_; ++i) {
                this->test_index.at(i) = vec_diff(
                    this->test_index.at(i), static_train_index
                    );
            }
            arma::uvec all_index {
                arma::regspace<arma::uvec>(0, n_obs - 1)
            };
            for (size_t i {0}; i < n_folds_; ++i) {
                this->train_index.push_back(
                    vec_diff(all_index, this->test_index.at(i))
                    );
            }
        }

        // helper function
        unsigned long get_n_folds() const
        {
            return n_folds_;
        }
        unsigned long get_n_obs() const
        {
            return n_obs_;
        }

    };
}

#endif /* CROSS_VALIDATION_H */
