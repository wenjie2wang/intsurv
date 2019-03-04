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

#ifndef LOGISTIC_REG_H
#define LOGISTIC_REG_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {

    // define class for inputs and outputs
    class LogisticReg
    {
    private:
        arma::mat x;
        arma::vec y;

    public:
        // constructors
        LogisticReg(const arma::mat& x_, const arma::vec& y_) :
            x(x_), y(y_) {}

        // function members
        inline double objective(const arma::vec& beta) const;
        inline arma::vec gradient(const arma::vec& beta) const;
        inline double objective(const arma::vec& beta, arma::vec& grad) const;
        inline arma::vec linkinv(const arma::vec& eta) const;

    };

    // define objective function
    inline double LogisticReg::objective(const arma::vec& beta) const
    {
        arma::vec x_beta {x * beta};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        arma::vec y_x_beta {y % x_beta};
        double res {
            arma::as_scalar(arma::sum(arma::log(1 + exp_x_beta) - y_x_beta))
        };
        return res;
    }

    // define gradient function
    inline arma::vec LogisticReg::gradient(const arma::vec& beta) const
    {
        arma::vec x_beta {x * beta};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        arma::vec res
            { x.t() * (exp_x_beta / (1 + exp_x_beta) - y)};
        return res;
    }

    // define objective function and overwrites graidient
    inline double LogisticReg::objective(const arma::vec& beta,
                                         arma::vec& grad) const
    {
        arma::vec x_beta {x * beta};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        grad = x.t() * (exp_x_beta / (1 + exp_x_beta) - y);
        arma::vec y_x_beta {y % x_beta};
        double res {
            arma::as_scalar(arma::sum(arma::log(1 + exp_x_beta) - y_x_beta))
        };
        return res;
    }

    // define inverse link function
    inline arma::vec LogisticReg::linkinv(const arma::vec& eta) const
    {
        return 1 / (1 + arma::exp(- eta));
    }

}


#endif
