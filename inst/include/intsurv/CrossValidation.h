//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2025  Wenjie Wang <wang@wwenjie.org>
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

namespace intsurv {

    typedef std::vector<arma::uvec> cv_index;
    typedef std::vector<cv_index> cv_indices;

    class CrossValidation
    {
    protected:

        unsigned int n_obs_;
        unsigned int n_folds_ = 5;

        // generate cross-validation indices
        // for given number of folds and number of observations
        inline cv_index get_cv_test_index(const unsigned int n_obs,
                                          const unsigned int n_folds) const
        {
            // number of observations must be at least two
            if (n_obs < 2) {
                throw std::range_error(
                    "Cross-validation needs at least two observations."
                    );
            }
            // number of folds is at most number of observations
            if (n_folds > n_obs) {
                throw std::range_error(
                    "Number of folds should be <= number of observations."
                    );
            }
            // define output
            cv_index out;
            // observation indices random permuted
            arma::uvec obs_idx { arma::randperm(n_obs) };
            // remaining number of observations
            size_t re_n_obs { n_obs };
            // determine the size of folds and indices one by one
            for (size_t i {0}; i < n_folds; ++i) {
                size_t fold_i { re_n_obs / (n_folds - i) };
                size_t j { n_obs - re_n_obs };
                arma::uvec idx_i { obs_idx.subvec(j, j + fold_i - 1) };
                out.push_back(idx_i);
                re_n_obs -= fold_i;
            }
            return out;
        }

        // random split
        inline cv_indices split(const int n_obs,
                                const int n_folds) const
        {
            cv_index test_idx, train_idx;
            test_idx = get_cv_test_index(n_obs, n_folds);
            arma::uvec all_index {
                arma::regspace<arma::uvec>(0, n_obs - 1)
            };
            for (size_t i {0}; i < n_folds_; ++i) {
                train_idx.push_back(
                    vec_diff(all_index, test_idx.at(i))
                    );
            }
            return { train_idx, test_idx };
        }

        // strata takes values from {0, ..., k}
        inline cv_indices stratified_split(const arma::uvec& strata,
                                           const unsigned int n_folds) const
        {
            const unsigned int n_strata { arma::max(strata) + 1 };
            // for the first strata
            arma::uvec k_idx { arma::find(strata == 0) };
            unsigned int n_k { k_idx.n_elem };
            cv_indices tmp_cv { split(n_k, n_folds) };
            cv_index train_index, test_index;
            for (size_t ii { 0 }; ii < n_folds; ++ii) {
                train_index.push_back(k_idx.elem(tmp_cv.at(0).at(ii)));
                test_index.push_back(k_idx.elem(tmp_cv.at(1).at(ii)));
            }
            // for the remaining strata
            for (size_t j { 1 }; j < n_strata; ++j) {
                k_idx = arma::find(strata == j);
                n_k = k_idx.n_elem;
                tmp_cv = split(n_k, n_folds);
                for (size_t ii { 0 }; ii < n_folds; ++ii) {
                    train_index.at(ii) = arma::join_cols(
                        train_index.at(ii),
                        k_idx.elem(tmp_cv.at(0).at(ii))
                        );
                    test_index.at(ii) = arma::join_cols(
                        test_index.at(ii),
                        k_idx.elem(tmp_cv.at(1).at(ii))
                        );
                }
            }
            cv_indices out { train_index, test_index };
            return out;
        }

    public:
        cv_index train_index_;
        cv_index test_index_;

        // default constructor
        CrossValidation();

        // explicit constructor
        explicit CrossValidation(const unsigned int n_obs) :
            n_obs_ { n_obs }
        {
            CrossValidation(n_obs_, 5, arma::uvec(), arma::uvec());
        }

        // a special constructor
        CrossValidation(const unsigned int n_obs,
                        const unsigned int n_folds,
                        const arma::uvec& static_train_index = arma::uvec(),
                        const arma::uvec& strata = arma::uvec()) :
            n_obs_ { n_obs },
            n_folds_ { n_folds }
        {
            const bool stratified { ! strata.is_empty() };
            cv_indices res;
            if (stratified) {
                res = stratified_split(strata, n_folds);
            } else {
                res = split(n_obs, n_folds);
            }
            train_index_ = res.at(0);
            test_index_ = res.at(1);
            if (! static_train_index.is_empty()) {
                // remove static train index from test index
                for (size_t i {0}; i < n_folds_; ++i) {
                    test_index_.at(i) = vec_diff(test_index_.at(i),
                                                 static_train_index);
                }
                for (size_t i {0}; i < n_folds_; ++i) {
                    train_index_.at(i) = vec_union(train_index_.at(i),
                                                   static_train_index);
                }
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
