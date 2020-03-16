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

        // explicit constructor
        explicit CrossValidation(const unsigned long n_obs) :
            n_obs_ { n_obs }
        {
            CrossValidation(n_obs_, 10);
        }

        // constructor
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
