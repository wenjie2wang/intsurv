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

#ifndef INTSURV_CONTROL_H
#define INTSURV_CONTROL_H

#include <RcppArmadillo.h>
#include <stdexcept>
#include "utils.h"

namespace intsurv {

    class Control
    {
    public:
        // common
        unsigned int max_iter_ { 1000 };
        double epsilon_ { 1e-4 };
        bool standardize_ { true };
        arma::vec start_ { arma::vec() };
        arma::vec offset_ { arma::vec() };
        unsigned int verbose_ { 0 };

        // for logistic models
        bool intercept_ { true };
        double pmin_ { 1e-5 };

        // for cure rate models
        unsigned int tail_completion_ { 1 };
        double tail_tau_ { - 1.0 };

        // for elastic-net models
        // common
        arma::vec penalty_factor_ { arma::vec() };
        bool varying_active_ { true };
        // single lambda's
        double l1_lambda_ { 1000.0 };
        double l2_lambda_ { 1000.0 };
        // solution paths
        unsigned int nlambda_ { 2 };
        double lambda_min_ratio_ { 1e-2 };
        double alpha_ { 1.0 };
        arma::vec lambda_ { arma::vec() };

        // default constructor
        Control() {}

        Control(const unsigned int max_iter,
                const double epsilon,
                const bool standardize = true,
                const unsigned int verbose = 0)
        {
            if (is_lt(epsilon, 0.0)) {
                throw std::range_error("The 'epsilon' cannot be negative.");
            }
            max_iter_ = max_iter;
            epsilon_ = epsilon;
            standardize_ = standardize;
            verbose_ = verbose;
        }

        // some setters
        Control* set_start(const arma::vec& start = arma::vec())
        {
            start_ = start;
            return this;
        }
        Control* set_offset(const arma::vec& offset = arma::vec())
        {
            offset_ = offset;
            return this;
        }
        Control* set_standardize(const bool standardize = true)
        {
            standardize_ = standardize;
            return this;
        }
        Control* set_verbose(const unsigned int verbose = 0)
        {
            verbose_ = verbose;
            return this;
        }

        // methods to set values for different models
        Control* logistic(const bool intercept = true,
                          const double pmin = 1e-5)
        {
            if (is_gt(pmin, 0.01) || is_le(pmin, 0.0)) {
                throw std::range_error(
                    "The 'pmin' must be between 0 and 0.01.");
            }
            intercept_ = intercept;
            pmin_ = pmin;
            return this;
        }
        Control* cure(const unsigned int tail_completion = 1,
                      const double tail_tau = -1.0)
        {
            tail_completion_ = tail_completion;
            tail_tau_ = tail_tau;
            return this;
        }
        Control* net(const arma::vec& penalty_factor = arma::vec(),
                     const bool varying_active = true)
        {
            if (! penalty_factor.is_empty() &&
                arma::any(penalty_factor < 0.0)) {
                throw std::range_error(
                    "The 'penalty_factor' cannot be negative.");
            }
            penalty_factor_ = penalty_factor;
            varying_active_ = varying_active;
            return this;
        }
        Control* net_fit(const double l1_lambda,
                         const double l2_lambda)
        {
            if (is_lt(l1_lambda, 0.0) || is_lt(l2_lambda, 0.0)) {
                throw std::range_error("The 'lambda' cannot be negative.");
            }
            l1_lambda_ = l1_lambda;
            l2_lambda_ = l2_lambda;
            return this;
        }
        Control* net_path(const unsigned int nlambda,
                          const double lambda_min_ratio,
                          const double alpha = 1.0,
                          const arma::vec lambda = arma::vec())
        {
            if (is_gt(lambda_min_ratio, 1.0) || is_le(lambda_min_ratio, 0.0)) {
                throw std::range_error(
                    "The 'lambda_min_ratio' must be between 0 and 1.");
            }
            if (is_gt(alpha, 1.0) || is_lt(alpha, 0.0)) {
                throw std::range_error(
                    "The 'alpha' must be between 0 and 1.");
            }
            if (! lambda.is_empty() && arma::any(lambda < 0.0)) {
                throw std::range_error("The 'lambda' cannot be negative.");
            }
            lambda_min_ratio_ = lambda_min_ratio;
            nlambda_ = nlambda;
            alpha_ = alpha;
            lambda_ = lambda;
            return this;
        }
    };

}

#endif /* INTSURV_CONTROL_H */
