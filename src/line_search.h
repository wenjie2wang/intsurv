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

#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

#include <RcppArmadillo.h>
// #include <armadillo>
#include <stdexcept>
#include "utils.h"

namespace Intsurv {

    // TODO: add the line search method based on

    // Mor\'e, Jorge J, & Thuente, D. J. (1994). Line search algorithms with
    // guaranteed sufficient decrease. ACM Transactions on Mathematical
    // Software (TOMS), 20(3), 286–307.

    //! @class Class for Control Parameters
    class control_line_search
    {
    public:
        double c1 {1e-4};
        double c2 {0.9};
        unsigned int max_iter {1000};
        double dec {0.5};
        double inc {2.1};
        unsigned int condition {2};
        double min_step {1e-20};
        double max_step {1e+20};
        // constructors
        control_line_search() {}
        control_line_search(double c1_, double c2_, unsigned int max_iter_,
                            double dec_, double inc_, unsigned int condition_,
                            double min_step_, double max_step_)
        {
            // check: 0 < c1 < c2 < 1
            if (c1_ > 0 && c1_ < c2_ && c2_ < 1) {
                c1 = c1_;
                c2 = c2_;
            } else {
                throw std::range_error(
                    "Error: '0 < c1 < c2 < 1' is not satisfied.");
            }
            // check: max_iter > 0
            if (max_iter_ == 0)
                throw std::range_error(
                    "Error: the 'max_iter' must be positive.");
            max_iter = max_iter_;
            // check: 0 < dec_ < 1 and inc_ > 1
            if (dec_ <= 0 || dec_ > 1)
                throw std::range_error(
                    "Error: '0 < dec < 1' is not satisfied.");
            if (inc <= 1)
                throw std::range_error(
                    "Error: 'inc > 1' is not satisfied.");
            dec = dec_;
            inc = inc_;
            // check: condition = 1, 2, or 3
            if (condition_ > 0 && condition_ < 4) {
                condition = condition_;
            } else {
                throw std::range_error(
                    "Error: the 'condition' has to be 1, 2, or 3.");
            }
            // check: min_step_ > 0 and max_step > 0
            if (min_step_ < 0)
                throw std::range_error(
                    "Error: the 'min_step' cannot be negative.");
            if (max_step_ <= 0)
                throw std::range_error(
                    "Error: the 'max_step' must be positive.");
            min_step = min_step_;
            max_step = max_step_;
        }
    };

    //! Function Template for Backtracking Line Search
    //!
    //! @param step The step size/length variable. It will be overwritten by the
    //! selected step size/length.
    //!
    //! @param x The solution from the current iteration.  It will be
    //! overwritten by the next solution based on search direction and the
    //! selected step size.
    //!
    //! @param object Objects with
    //! "objective" function: std::function<double(arma::vec)>
    //! having argument for overwriting gradient.
    //!
    //! @param drt The search direction at the current iteration.
    //!
    //! @param control A object from the `control_line_search` class
    //! specifying all the control paramaters for line search.
    //!
    //! @reference
    //! backtracking line search from the lbfgs C-library
    template <typename T>
    inline int LineSearchBacktracking(double& step,
                                      arma::vec& x,
                                      const T& object,
                                      const arma::vec& drt,
                                      const control_line_search& control)
    {
        // initialize step size to be 1
        step = 1.0;
        // the function value at the current x
        const arma::vec xp {x};
        arma::vec f_grad {x};
        const double fx_init {object.objective(xp, f_grad)};
        // the projection of gradient on the search direction
        const double dg_init {
            vec2num(crossprod(f_grad, drt))
        };
        if (dg_init > 0)
            std::logic_error("Error: the moving direction is not decreasing.");
        const double dg_test {control.c1 * dg_init};
        double size_factor, f_x, dg;

        for (size_t i {0}; i < control.max_iter; ++i) {
            // update x, objective function value, and it gradient
            x = xp + step * drt;
            f_x = object.objective(x, f_grad);

            // check the first Wolfe condition
            if (f_x > fx_init + step * dg_test) {
                // if not satisfied, decrease the step size
                size_factor = control.dec;
            } else {
                // otherwise, the first Wolfe (Armijo) condition is satisfied
                if (control.condition == 1)
                    return 0;
                // update projection
                dg = vec2num(crossprod(f_grad, drt));
                // check the second Wolfe (curvature) condition
                if (dg < control.c2 * dg_init) {
                    // if not satisfied, increase the step size
                    size_factor = control.inc;
                } else {
                    // otherwise, the regular Wolfe conditions are satisfied
                    if (control.condition == 2)
                        return 0;
                    // check strong Wolfe conditions
                    if (dg > - control.c2 * dg_init) {
                        // if not satisfied, decrease the step size
                        // exclude points far from stationary points
                        size_factor = control.dec;
                    } else {
                        // strong Wolfe conditions are satisfied
                        return 0;
                    }
                }
            }
            if (step < control.min_step) {
                throw std::runtime_error(
                    "Error: step size is less than the min value.");
            }
            if (step > control.max_step) {
                throw std::runtime_error(
                    "Error: step size is greater than the max value.");
            }
            step *= size_factor;
        }
        throw std::runtime_error(
            "Error: failed to find step size with iterations.");
    }

}

#endif
