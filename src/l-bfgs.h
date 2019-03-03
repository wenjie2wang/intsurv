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

#ifndef L_BFGS_H
#define L_BFGS_H

#include <stdexcept>
// #include <armadillo>
#include <RcppArmadillo.h>

#include "utils.hpp"
#include "hessian_approximation.hpp"
#include "line_search.hpp"

namespace Intsurv {

    class control_lbfgs
    {
    public:
        // data members
        unsigned int m {8};
        double epsilon {1e-5};
        unsigned int max_iter {1000};
        control_line_search line_search {control_line_search()};

        // constructors
        control_lbfgs() {}
        control_lbfgs(unsigned int m_,
                      double epsilon_,
                      unsigned int max_iter_,
                      control_line_search line_search_)
        {
            // m usually takes from 3 to 20
            if (m_ == 0)
                throw std::range_error("Error: the 'm' must be positive.");
            m = m_;
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

    // reference:
    // Nocedal, J., & Wright, S. (2006). Numerical Optimization. :
    // Springer-Verlag New York.

    // algorithm 7.4 on Page 178, Chapter 7.2
    inline arma::vec twoLoop(const int& k,
                             const int& m,
                             arma::vec q,
                             const arma::mat& s_mat,
                             const arma::mat& y_mat,
                             const arma::mat& H_k)
    {
        unsigned long int n_param {s_mat.n_rows};
        int j {0};
        arma::vec y_i {arma::vec(n_param)};
        arma::vec s_i {arma::vec(n_param)};
        double denom_i {0};

        // define cache variables
        arma::rowvec rho_rvec {arma::rowvec(m)};
        arma::rowvec alpha_rvec {arma::rowvec(m)};
        double beta {0};

        // define the output
        arma::vec r {arma::vec(n_param)};

        // the first loop, i = k - 1, k - 2, ..., k - m.
        for (int i {k - 1}; i >= std::max(0, k - m); --i) {
            j = i % m;
            y_i = y_mat.col(j);
            s_i = s_mat.col(j);
            // compute rho_i
            denom_i = vec2num(crossprod(s_i, y_i));
            if (! isAlmostEqual(denom_i, 0.0)) {
                rho_rvec(j) = 1 / denom_i;
                alpha_rvec(j) = rho_rvec(j) * vec2num(crossprod(s_i, q));
                q -= alpha_rvec(j) * y_i;
            }
        }
        // compute r from the first loop
        r = H_k * q;
        // the second loop, i = k - m, k - m + 1, ..., k - 1.
        for (int i {std::max(0, k - m)}; i < k; ++i) {
            j = i % m;
            y_i = y_mat.col(j);
            s_i = s_mat.col(j);
            beta = rho_rvec(j) * vec2num(crossprod(y_i, r));
            r += s_i * (alpha_rvec(j) - beta);
        }
        return r;
    }

    //! Function Template for Backtracking Line Search
    //!
    //! @param x The initial guesses or starting values of solution.  It will
    //! be overwritten by the solution.
    //!
    //! @param object Objects with
    //! "objective" function: std::function<double(arma::vec)>
    //! having argument for overwriting gradient.
    //!
    //! @param control A object from the `control_lbfgs` class
    //! specifying all the control paramaters.
    //!
    //! @reference
    //! Nocedal, J., & Wright, S. (2006). Numerical Optimization. :
    //! Springer-Verlag New York.
    template <typename T>
    inline int lbfgs(arma::vec& x,
                     const T& object,
                     const control_lbfgs& control)
    {
        unsigned long int n_param {x.n_rows};
        int status_code {0};
        // initialize inverse hessian approximation H_k
        arma::mat H_k {arma::eye(n_param, n_param)};
        double alpha_k {1};
        // initialize for k = 0
        arma::vec g_k, p_k;
        // initialize variables for recursion formula
        arma::mat s_mat {arma::zeros(n_param, control.m)};
        arma::mat y_mat {s_mat};
        int j0, j1;
        // line search control object
        control_line_search control_ls {control.line_search};

        // algorithm 7.5 on Page 179, Chapter 7.2
        for (size_t k {0}; k < control.max_iter; ++k) {
            // update gradient function
            g_k = object.gradient(x);
            // check whether convergence is reached by relative difference
            if (norm(g_k) < control.epsilon * std::max(norm(x), 1.0)) {
                return status_code;
            }
            // compute column index
            j0 = (k - 1) % control.m;
            j1 = k % control.m;
            // initialize H_k^0
            H_k = initialize_inverse_hessian(y_mat.col(j0), s_mat.col(j0));
            // call algorithm 7.4 to compute p_k
            p_k = - twoLoop(k, control.m, g_k, s_mat, y_mat, H_k);
            // compute alpha_k
            status_code = LineSearchBacktracking(
                alpha_k, x, object, p_k, control_ls);
            if (status_code != 0)
                throw std::runtime_error("Error: line search failed.");
            // update s_mat and y_mat
            s_mat.col(j1) = alpha_k * p_k;
            y_mat.col(j1) = object.gradient(x) - g_k;
        }
        if (status_code != 0) {
            throw std::runtime_error("Error: L-BFGS failed.");
        }
        return status_code;
    }

}

#endif
