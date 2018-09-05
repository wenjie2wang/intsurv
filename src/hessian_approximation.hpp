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

#ifndef HESSIAN_APPROXIMATION_H
#define HESSIAN_APPROXIMATION_H

#include <stdexcept>
#include <RcppArmadillo.h>
#include "utils.hpp"

namespace Intsurv {

    // reference:
    // Nocedal, J., & Wright, S. (2006). Numerical Optimization. :
    // Springer-Verlag New York.

    // function that updates the inverse of hessian approximation directly
    inline int update_inverse_hessian(arma::mat& H_k,
                                      const arma::vec& y_k,
                                      const arma::vec& s_k)
    {
        double rho_k_denom {arma::as_scalar(crossprod(y_k, s_k))};
        if (isAlmostEqual(rho_k_denom, 0.0)) {
            return - 1;
        }
        double rho_k {1 / rho_k_denom};
        arma::mat s_y_t {tcrossprod(s_k, y_k)};
        arma::mat I_mat {arma::eye<arma::mat>(size(s_y_t))};
        H_k = (I_mat - rho_k * s_y_t) * H_k * (I_mat - rho_k * s_y_t.t()) +
            rho_k * tcrossprod(s_k);
        return 0;
    }

    // Chapter 7.2, Page 179
    // function that initializes the inverse of hessian approximation
    inline arma::mat initialize_inverse_hessian(const arma::vec& y_k,
                                                const arma::vec& s_k)
    {
        double numer {vec2num(crossprod(s_k, y_k))};
        double denom {vec2num(crossprod(y_k))};
        if (isAlmostEqual(denom, 0.0)) {
            // initialize to a identity matrix
            return arma::eye(y_k.n_rows, y_k.n_rows);
        } else {
            double gamma_k {std::exp(std::log(numer) - std::log(denom))};
            return arma::eye(y_k.n_rows, y_k.n_rows) * gamma_k;
        }
    }

}

#endif
