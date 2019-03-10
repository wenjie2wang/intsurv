//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
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
        bool intercept;
        arma::rowvec cmd_lowerbound;

    public:
        arma::vec coef;
        double negative_logL {0};

        // constructors
        LogisticReg(const arma::mat& x_, const arma::vec& y_,
                    const bool intercept_ = true)
        {
            intercept = intercept_;
            if (intercept_) {
                x = arma::join_horiz(arma::ones(x_.n_rows), x_);
            } else {
                x = x_;
            }
            y = y_;
        }

        // function members
        inline arma::vec linkinv(const arma::vec& eta) const;

        inline double objective(const arma::vec& beta) const;

        inline arma::vec gradient(const arma::vec& beta) const;

        inline double gradient(const arma::vec& beta,
                               const unsigned int k) const;

        inline double objective(const arma::vec& beta, arma::vec& grad) const;

        // fit regular logistic regression model
        void fit(const arma::vec& start,
                 const unsigned int& max_iter,
                 const double& rel_tol);

        // compute cov lowerbound used in regularied model
        void compute_cmd_lowerbound(bool force_update);

        // update step for regularized logistic regression model
        void regularized_fit_update(arma::vec& beta,
                                    arma::uvec& is_active,
                                    const arma::rowvec& b_vec,
                                    const arma::vec& penalty,
                                    const bool& update_active);

        // fit regularized logistic regression model
        void regularized_fit(const double& lambda,
                             const arma::vec& penalty_factor,
                             const arma::vec& start,
                             const unsigned int& max_iter,
                             const double& rel_tol);

        // some simple functions
        unsigned int sample_size() const
        {
            return y.n_elem;
        }

        // helper function members to access some private members
        arma::mat get_x() const { return x; }
        arma::vec get_y() const { return y; }

    };

    // define inverse link function
    inline arma::vec LogisticReg::linkinv(const arma::vec& eta) const
    {
        return 1 / (1 + arma::exp(- eta));
    }

    // define objective function (negative log-likehood function)
    inline double LogisticReg::objective(const arma::vec& beta) const
    {
        arma::vec x_beta { x * beta };
        arma::vec exp_x_beta { arma::exp(x_beta) };
        arma::vec y_x_beta { y % x_beta };
        return arma::as_scalar(arma::sum(arma::log(1 + exp_x_beta) - y_x_beta));
    }

    // define gradient function
    inline arma::vec LogisticReg::gradient(const arma::vec& beta) const
    {
        arma::vec y_hat { linkinv(x * beta) };
        return x.t() * (y_hat - y);
    }
    // define gradient function at k-th dimension
    inline double LogisticReg::gradient(const arma::vec& beta,
                                        const unsigned int k) const
    {
        arma::vec y_hat { linkinv(x * beta) };
        return arma::sum((y_hat - y) % x.col(k));
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

    // fitting regular logistic model by monotonic quadratic approximation
    // algorithm non-integer y vector is allowed
    // reference: Bohning and Lindsay (1988) SIAM
    void LogisticReg::fit(const arma::vec& start = 0,
                          const unsigned int& max_iter = 1000,
                          const double& rel_tol = 1e-6)
    {
        arma::vec beta0 { arma::zeros(x.n_cols) };
        if (start.n_elem == x.n_cols) {
            beta0 = start;
        }
        arma::vec beta { beta0 };
        arma::vec eta { arma::zeros(y.n_elem) };
        arma::vec y_hat { eta };
        arma::mat iter_mat { 4 * arma::inv_sympd(x.t() * x) * x.t() };
        size_t i {0};
        while (i < max_iter) {
            eta = x * beta0;
            y_hat = linkinv(eta);
            beta = beta0 + iter_mat * (y - y_hat);
            if (Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
                break;
            }
            // update beta
            beta0 = beta;
            i++;
        }
        this->coef = beta;
    }

    // compute CMD lowerbound vector
    void LogisticReg::compute_cmd_lowerbound(bool force_update = false)
    {
        if (force_update || cmd_lowerbound.is_empty()) {
            this->cmd_lowerbound = arma::sum(arma::square(x), 0) /
                (4 * y.n_elem);
        }
    }

    // run one cycle of coordinate descent over a given active set
    void LogisticReg::regularized_fit_update(arma::vec& beta,
                                             arma::uvec& is_active,
                                             const arma::rowvec& b_vec,
                                             const arma::vec& penalty,
                                             const bool& update_active = false)
    {
        double dlj { 0 };
        unsigned int n { y.n_elem };
        double n_y { static_cast<double>(n) };
        for (size_t j {0}; j < beta.n_elem; ++j) {
            if (is_active[j]) {
                dlj = this->gradient(beta, j) / n_y;
                // update beta
                beta[j] = Intsurv::soft_threshold(
                    b_vec[j] * beta[j] - dlj, penalty[j]) / b_vec[j];
                if (update_active) {
                    // check if it has been shrinkaged to zero
                    if (Intsurv::isAlmostEqual(beta[j], 0)) {
                        is_active[j] = 0;
                    } else {
                        is_active[j] = 1;
                    }
                }
            }
        }
    }

    // regularized logistic model by coordinate-majorization-descent algorithm
    void LogisticReg::regularized_fit(const double& lambda = 0,
                                      const arma::vec& penalty_factor = 0,
                                      const arma::vec& start = 0,
                                      const unsigned int& max_iter = 1000,
                                      const double& rel_tol = 1e-6)
    {
        arma::vec beta0 { arma::zeros(x.n_cols) };
        if (start.n_elem == x.n_cols) {
            beta0 = start;
        }
        arma::vec beta { beta0 };
        arma::uvec is_active { arma::ones<arma::uvec>(x.n_cols) };
        arma::uvec is_active_stored { is_active };
        // set penalty terms
        unsigned int int_intercept { static_cast<unsigned int>(intercept) };
        arma::vec penalty { arma::ones(x.n_cols - int_intercept) };
        if (penalty_factor.n_elem == penalty.n_elem) {
            penalty = penalty_factor / arma::sum(penalty_factor);
        }
        penalty *= lambda;
        if (this->intercept) {
            penalty = arma::join_vert(arma::zeros(1), penalty);
        }
        // compute lowerbound vector used in CMD algorithm
        this->compute_cmd_lowerbound();
        arma::rowvec b_vec { this->cmd_lowerbound };
        size_t i {0};

        // use active-set if p > n ("helps when p >> n")
        if (x.n_cols > x.n_rows) {
            size_t ii {0};
            while (i < max_iter) {
                // cycles over the active set
                while (ii < max_iter) {
                    regularized_fit_update(beta, is_active,
                                           b_vec, penalty, true);
                    if (Intsurv::isAlmostEqual(Intsurv::l2_norm(beta), 0) ||
                        Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                is_active_stored = is_active;
                // run a full cycle over the converged beta
                is_active = arma::ones<arma::uvec>(x.n_cols);
                regularized_fit_update(beta, is_active,
                                       b_vec, penalty, true);
                // check two active sets coincide
                if (arma::sum(arma::abs(is_active - is_active_stored)) > 0) {
                    // if different, repeat this process
                    ii = 0;
                    i++;
                } else {
                    break;
                }
            }
        } else {
            // regular coordinate descent
            while (i < max_iter) {
                regularized_fit_update(beta, is_active,
                                       b_vec, penalty, false);
                if (Intsurv::isAlmostEqual(Intsurv::l2_norm(beta), 0) ||
                    Intsurv::rel_l2_norm(beta, beta0) < rel_tol) {
                    break;
                }
                beta0 = beta;
                i++;
            }
        }
        this->coef = beta;
    }

}


#endif
