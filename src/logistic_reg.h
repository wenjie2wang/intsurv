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
        arma::mat x;            // (standardized) x
        arma::vec y;
        bool intercept;
        bool standardize;       // is x standardized
        arma::rowvec x_center;  // the column center of x
        arma::rowvec x_scale;   // the scale of x


        // iteration matrix in Bohning and Lindsay (1988)
        arma::mat bl_iter_mat;
        // for regularized regression
        arma::rowvec cmd_lowerbound;

    public:
        arma::vec coef0;        // coef before rescaling
        arma::vec coef;         // coef rescaled
        double neg_LogL;        // negative likelihood

        // regularized
        arma::mat reg_coef;     // coef rescaled
        arma::vec reg_lambda;   // lambda sequence
        arma::vec reg_neg_logL; // negative likelihood

        // constructors
        LogisticReg(const arma::mat& x_,
                    const arma::vec& y_,
                    const bool intercept_ = true,
                    const bool standardize_ = true)
        {
            intercept = intercept_;
            standardize = standardize_;
            x = x_;
            if (standardize) {
                x_center = arma::mean(x);
                x_scale = arma::stddev(x, 1);
                for (size_t j {0}; j < x.n_cols; ++j) {
                    if (x_scale(j) > 0) {
                        x.col(j) = (x.col(j) - x_center(j)) / x_scale(j);
                    } else {
                        throw std::range_error(
                            "The design 'x' contains constant column."
                            );
                    }
                }
            }
            if (intercept) {
                x = arma::join_horiz(arma::ones(x.n_rows), x);
            }
            y = y_;
        }

        // function members
        // transfer coef for standardized data to coef for non-standardized data
        inline void rescale_coef()
        {
            this->coef = coef0;
            if (this->standardize) {
                if (this->intercept) {
                    arma::uvec non_int_ind {
                        arma::regspace<arma::uvec>(1, coef0.n_elem - 1)
                    };
                    this->coef[0] = coef0(0) -
                        arma::as_scalar((x_center / x_scale) *
                                        coef0.elem(non_int_ind));
                }
                for (size_t j {1}; j < coef0.n_elem; ++j) {
                    this->coef[j] = coef0[j] / x_scale[j - 1];
                }
            }
        }
        // transfer coef for non-standardized data to coef for standardized data
        inline arma::vec rev_rescale_coef(const arma::vec& beta) const
        {
            arma::vec beta0 { beta };
            double tmp {0};
            for (size_t j {1}; j < beta.n_elem; ++j) {
                beta0(j) *= x_scale(j - 1);
                tmp += beta(j) * x_center(j - 1);
            }
            beta0(0) += tmp;
            return beta0;
        }

        inline arma::vec linkinv(const arma::vec& eta) const;

        // here beta is coef vector for non-standardized data
        inline arma::vec predict(const arma::vec& beta) const;

        inline double objective(const arma::vec& beta) const;

        inline arma::vec gradient(const arma::vec& beta) const;

        inline double gradient(const arma::vec& beta,
                               const unsigned int k) const;

        inline double objective(const arma::vec& beta, arma::vec& grad) const;

        // compute iteration matrix in Bohning and Lindsay (1988)
        inline void compute_bl_iter_mat(bool force_update);

        // fit regular logistic regression model
        inline void fit(const arma::vec& start,
                        const unsigned int& max_iter,
                        const double& rel_tol);

        // compute cov lowerbound used in regularied model
        inline void compute_cmd_lowerbound(bool force_update);

        // update step for regularized logistic regression model
        inline void regularized_fit_update(arma::vec& beta,
                                           arma::uvec& is_active,
                                           const arma::vec& penalty,
                                           const bool& update_active);

        inline void reg_active_fit(arma::vec& beta,
                                   const arma::uvec& is_active,
                                   const arma::vec& penalty,
                                   const bool& varying_active_set,
                                   const unsigned int& max_iter,
                                   const double& rel_tol);

        // fit regularized logistic regression model
        // for a perticular lambda
        inline void regularized_fit(const double& lambda,
                                    const arma::vec& penalty_factor,
                                    const arma::vec& start,
                                    const unsigned int& max_iter,
                                    const double& rel_tol);

        // overload for a sequence of lambda's
        inline void regularized_fit(arma::vec lambda,
                                    const unsigned int& nlambda,
                                    double lambda_min_ratio,
                                    const arma::vec& penalty_factor,
                                    const unsigned int& max_iter,
                                    const double& rel_tol);

        // function that helps update y
        inline void update_y(const arma::vec& y_) { this->y = y_; }

        // some simple functions
        inline unsigned int sample_size() const
        {
            return y.n_elem;
        }

        // helper function members to access some private members
        inline arma::mat get_x(bool include_intercept = true) const
        {
            arma::mat out {this->x};
            if (include_intercept && this-> intercept) {
                out.shed_col(0);
            }
            if (this->standardize) {
                for (size_t j {0}; j < out.n_cols; ++j) {
                    out.col(j) = x_scale(j) * out.col(j) + x_center(j);
                }
            }
            return out;
        }
        inline arma::vec get_y() const { return y; }

    };

    // define inverse link function
    inline arma::vec LogisticReg::linkinv(const arma::vec& beta) const
    {
        return 1 / (1 + arma::exp(- mat2vec(x * beta)));
    }
    inline arma::vec LogisticReg::predict(const arma::vec& beta) const
    {
        arma::vec beta0 { beta };
        if (this->standardize) {
            beta0 = rev_rescale_coef(beta0);
        }
        return linkinv(beta0);
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
        return x.t() * (linkinv(beta) - y);
    }
    // define gradient function at k-th dimension
    inline double LogisticReg::gradient(const arma::vec& beta,
                                        const unsigned int k) const
    {
        return arma::sum((linkinv(beta) - y) % x.col(k));
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

    // compute iteration matrix in Bohning and Lindsay (1988)
    inline void LogisticReg::compute_bl_iter_mat(bool force_update = false)
    {
        if (force_update || this->bl_iter_mat.is_empty()) {
            this->bl_iter_mat = 4 * arma::inv_sympd(x.t() * x) * x.t();
        }
    }

    // fitting regular logistic model by monotonic quadratic approximation
    // algorithm non-integer y vector is allowed
    // reference: Bohning and Lindsay (1988) SIAM
    inline void LogisticReg::fit(const arma::vec& start = 0,
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
        this->compute_bl_iter_mat();
        size_t i {0};
        while (i < max_iter) {
            y_hat = linkinv(beta0);
            beta = beta0 + this->bl_iter_mat * (y - y_hat);
            if (rel_l2_norm(beta, beta0) < rel_tol) {
                break;
            }
            // update beta
            beta0 = beta;
            i++;
        }
        this->coef0 = beta;
        // rescale coef back
        this->rescale_coef();
    }

    // compute CMD lowerbound vector
    inline void LogisticReg::compute_cmd_lowerbound(bool force_update = false)
    {
        if (force_update || cmd_lowerbound.is_empty()) {
            this->cmd_lowerbound = arma::sum(arma::square(x), 0) /
                (4 * x.n_rows);
        }
    }

    // run one cycle of coordinate descent over a given active set
    inline void LogisticReg::regularized_fit_update(
        arma::vec& beta,
        arma::uvec& is_active,
        const arma::vec& penalty,
        const bool& update_active = false
        )
    {
        // compute lowerbound vector used in CMD algorithm
        this->compute_cmd_lowerbound();
        double dlj { 0 };
        for (size_t j {0}; j < beta.n_elem; ++j) {
            if (is_active[j]) {
                dlj = this->gradient(beta, j) / y.n_elem;
                // update beta
                beta[j] = soft_threshold(
                    cmd_lowerbound[j] * beta[j] - dlj, penalty[j]) /
                    cmd_lowerbound[j];
                if (update_active) {
                    // check if it has been shrinkaged to zero
                    if (isAlmostEqual(beta[j], 0)) {
                        is_active[j] = 0;
                    } else {
                        is_active[j] = 1;
                    }
                }
            }
        }
    }

    // run a complete cycle of CMD for a given active set and lambda
    inline void LogisticReg::reg_active_fit(
        arma::vec& beta,
        const arma::uvec& is_active,
        const arma::vec& penalty,
        const bool& varying_active_set = false,
        const unsigned int& max_iter = 1000,
        const double& rel_tol = 1e-6
        )
    {
        size_t i {0};
        arma::vec beta0 { beta };
        arma::uvec is_active_stored { is_active };

        // use active-set if p > n ("helps when p >> n")
        if (varying_active_set) {
            arma::uvec is_active_new { is_active };
            size_t ii {0};
            while (i < max_iter) {
                // cycles over the active set
                while (ii < max_iter) {
                    regularized_fit_update(beta, is_active_stored,
                                           penalty, true);
                    if (rel_l2_norm(beta, beta0) < rel_tol) {
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                // run a full cycle over the converged beta
                regularized_fit_update(beta, is_active_new, penalty, true);
                // check two active sets coincide
                if (arma::sum(arma::abs(is_active_new - is_active_stored))) {
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
                regularized_fit_update(beta, is_active_stored, penalty, false);
                if (rel_l2_norm(beta, beta0) < rel_tol) {
                    break;
                }
                beta0 = beta;
                i++;
            }
        }
    }

    // regularized logistic model by coordinate-majorization-descent algorithm
    // for a perticular lambda
    inline void LogisticReg::regularized_fit(
        const double& lambda = 0,
        const arma::vec& penalty_factor = 0,
        const arma::vec& start = 0,
        const unsigned int& max_iter = 1000,
        const double& rel_tol = 1e-6
        )
    {
        // set penalty terms
        unsigned int int_intercept { static_cast<unsigned int>(intercept) };
        unsigned int n_predictor { x.n_cols - int_intercept };
        if (n_predictor < 1) {
            throw std::range_error("Predictors not found.");
        }
        arma::vec penalty { arma::ones(n_predictor) };
        if (penalty_factor.n_elem == penalty.n_elem) {
            penalty = penalty_factor * x.n_cols / arma::sum(penalty_factor);
        }

        arma::vec beta { arma::zeros(x.n_cols) };
        arma::vec grad_zero { arma::abs(this->gradient(beta)) };
        arma::vec grad_beta { grad_zero }, strong_rhs { grad_beta };

        // large enough lambda for all-zero coef (except intercept)
        double lambda_max {
            arma::max(grad_zero.tail(n_predictor) /
                      penalty) / this->x.n_rows
        };

        if (this->intercept) {
            penalty = arma::join_vert(arma::zeros(1), penalty);
        }

        // get solution of lambda_max for a warm start
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x.n_cols) };
        if (this->intercept) {
            // only needs to estimate intercept
            is_active_strong(0) = 1;
            regularized_fit_update(beta, is_active_strong,
                                   penalty * lambda_max, false);
        }
        this->coef0 = beta;
        // rescale coef back
        this->rescale_coef();

        // early exit for lambda greater than lambda_max
        if (lambda >= lambda_max) {
            return;
        }

        // use the input start if correctly specified
        if (start.n_elem == x.n_cols) {
            beta = start;
        }

        // update active set by strong rule
        grad_beta = arma::abs(this->gradient(this->coef0)) / this->x.n_rows;
        strong_rhs = (2 * lambda - lambda_max) * penalty;

        for (size_t j {1}; j < n_predictor + 1; ++j) {
            if (grad_beta(j) > strong_rhs(j)) {
                is_active_strong(j) = 1;
            } else {
                beta(j) = 0;
            }
        }
        arma::uvec is_active_strong_new { is_active_strong };

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x.n_cols > x.n_rows) {
            varying_active_set = true;
        }

        strong_rhs = lambda * penalty;
        bool kkt_failed { true };
        // eventually, strong rule will guess correctly
        while (kkt_failed) {
            // update beta
            reg_active_fit(beta, is_active_strong, strong_rhs,
                           varying_active_set, max_iter, rel_tol);
            // check kkt condition
            for (size_t j {1}; j < n_predictor + 1; ++j) {
                if (is_active_strong(j)) {
                    continue;
                }
                if (std::abs(this->gradient(beta, j)) / this->x.n_rows >
                    strong_rhs(j)) {
                    // update active set
                    is_active_strong_new(j) = 1;
                }
            }
            if (arma::sum(arma::abs(is_active_strong - is_active_strong_new))) {
                is_active_strong = is_active_strong_new;
            } else {
                kkt_failed = false;
            }
        }
        this->coef0 = beta;
        // rescale coef back
        this->rescale_coef();
    }


    // for a sequence of lambda's
    inline void LogisticReg::regularized_fit(
        arma::vec lambda = 0,
        const unsigned int& nlambda = 100,
        double lambda_min_ratio = 1e-4,
        const arma::vec& penalty_factor = 0,
        const unsigned int& max_iter = 1000,
        const double& rel_tol = 1e-6
        )
    {
        // set penalty terms
        unsigned int int_intercept { static_cast<unsigned int>(intercept) };
        unsigned int n_predictor { x.n_cols - int_intercept };
        if (n_predictor < 1) {
            throw std::range_error("Predictors not found.");
        }
        arma::vec penalty { arma::ones(n_predictor) };
        if (penalty_factor.n_elem == penalty.n_elem) {
            penalty = penalty_factor * x.n_cols / arma::sum(penalty_factor);
        }

        // construct lambda sequence
        arma::vec beta { arma::zeros(x.n_cols) };
        arma::vec grad_zero { arma::abs(this->gradient(beta)) };
        double lambda_max {
            arma::max(grad_zero.tail(penalty.n_elem) /
                      penalty) / this->x.n_rows
        };
        double log_lambda_max { std::log(lambda_max) };

        // take unique lambda and sort descendingly
        lambda = arma::reverse(arma::unique(lambda));
        arma::vec lambda_seq {
            arma::zeros(std::max(nlambda, lambda.n_elem))
                };
        if (lambda.n_elem == 1 && nlambda > 1) {
            if (this->x.n_cols > this->x.n_rows) {
                lambda_min_ratio = std::sqrt(lambda_min_ratio);
            }
            lambda_seq = arma::exp(
                arma::linspace(log_lambda_max,
                               log_lambda_max + std::log(lambda_min_ratio),
                               nlambda)
                );
        } else {
            lambda_seq = lambda;
        }
        this->reg_lambda = lambda_seq;

        // update penalty for intercept
        if (this->intercept) {
            penalty = arma::join_vert(arma::zeros(1), penalty);
        }

        // set up the estimate matrix
        arma::mat beta_mat { arma::zeros(x.n_cols, lambda_seq.n_elem) };
        arma::vec grad_beta { grad_zero }, strong_rhs { grad_zero };

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x.n_cols) };

        // get solution of lambda_max for a warm start
        if (this->intercept) {
            // only needs to estimate intercept
            is_active_strong(0) = 1;
            regularized_fit_update(beta, is_active_strong,
                                   penalty * lambda_max, false);
        }
        this->coef0 = beta;
        // rescale coef back
        this->rescale_coef();

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x.n_cols > x.n_rows) {
            varying_active_set = true;
        }

        // outer loop for the lambda sequence
        for (size_t k {0}; k < lambda_seq.n_elem; ++k) {
            // early exit for large lambda greater than lambda_max
            if (lambda_seq(k) >= lambda_max) {
                beta_mat.col(k) = this->coef;
                continue;
            }
            // update acitve set by strong rule (for lambda < lamda_max)
            grad_beta = arma::abs(this->gradient(beta)) / this->x.n_rows;
            if (k == 0) {
                // use lambda_max
                strong_rhs = (2 * lambda_seq(k) - lambda_max) * penalty;
            } else {
                // use the last lambda
                strong_rhs = (2 * lambda_seq(k) - lambda_seq(k - 1)) * penalty;
            }
            for (size_t j {1}; j < n_predictor + 1; ++j) {
                if (grad_beta(j) > strong_rhs(j)) {
                    is_active_strong(j) = 1;
                } else {
                    beta(j) = 0;
                }
            }
            arma::uvec is_active_strong_new { is_active_strong };
            strong_rhs = lambda_seq(k) * penalty;
            bool kkt_failed { true };
            // eventually, strong rule will guess correctly
            while (kkt_failed) {
                // update beta
                reg_active_fit(beta, is_active_strong, strong_rhs,
                               varying_active_set, max_iter, rel_tol);
                // check kkt condition
                grad_beta = arma::abs(this->gradient(beta)) / this->x.n_rows;
                for (size_t j {1}; j < n_predictor + 1; ++j) {
                    if (is_active_strong(j)) {
                        continue;
                    }
                    if (grad_beta(j) > strong_rhs(j)) {
                        // update active set
                        is_active_strong_new(j) = 1;
                    }
                }
                if (arma::sum(arma::abs(is_active_strong -
                                        is_active_strong_new))) {
                    is_active_strong = is_active_strong_new;
                } else {
                    kkt_failed = false;
                }
            }
            this->coef0 = beta;
            // rescale coef back
            this->rescale_coef();
            beta_mat.col(k) = this->coef;
        }
        this->reg_coef = beta_mat;
    }

}


#endif
