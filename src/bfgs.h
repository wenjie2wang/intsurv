//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
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

#ifndef BFGS_H
#define BFGS_H

#include <stdexcept>
#include <RcppArmadillo.h>
// #include <armadillo>

#include "utils.h"
#include "hessian_approximation.h"
#include "line_search.h"

namespace Intsurv {

    class control_bfgs
    {
    public:
        // data members
        double epsilon {1e-5};
        unsigned int max_iter {1000};
        control_line_search line_search {control_line_search()};

        // constructors
        control_bfgs() {}
        control_bfgs(double epsilon_,
                     unsigned int max_iter_,
                     control_line_search line_search_)
        {
            if (epsilon_ <= 0)
                throw std::range_error(
                    "Error: the 'epsilon' must be positive.");
            epsilon = epsilon_;
            if (max_iter == 0)
                throw std::range_error(
                    "Error: the 'max_iter' must be positive.");
            max_iter = max_iter_;
            // control for line search
            line_search = line_search_;
        }
    };


    template <typename T>
    int bfgs(arma::vec& x, const T& object,
             const control_bfgs& control)
    {
        unsigned int n_parm {x.n_rows};
        int status_code {0};
        arma::mat H_k {arma::eye(n_parm, n_parm)};
        double alpha {1};
        arma::vec g_k, p_k, s_k, y_k;
        // line search control object
        control_line_search control_ls {control.line_search};

        for (size_t i {0}; i < control.max_iter; ++i) {
            g_k = object.gradient(x);
            if (norm(g_k) < control.epsilon * norm(x)) {
                return status_code;
            }
            p_k = - H_k * g_k;
            status_code = LineSearchBacktracking(
                alpha, x, object, p_k, control_ls);
            s_k = alpha * p_k;
            y_k = object.gradient(x) - g_k;
            status_code = update_inverse_hessian(H_k, y_k, s_k);
        }
        if (status_code != 0) {
            std::runtime_error("Error: BFGS failed.");
        }
        return status_code;
    }

}


#endif
