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
        arma::vec coef0;        // coef before rescaling

        // iteration matrix in Bohning and Lindsay (1988)
        arma::mat bl_iter_mat;
        // for regularized coordinate majorization descent
        arma::rowvec cmd_lowerbound;

    public:
        unsigned int nObs;           // number of observations
        arma::vec l1_penalty_factor; // adaptive weights for lasso penalty
        double l1_lambda_max;   // the "big enough" lambda => zero coef

        // for a sinle l1_lambda and l2_lambda
        double l1_lambda;       // tuning parameter for lasso penalty
        double l2_lambda;       // tuning parameter for ridge penalty
        arma::vec coef;         // coef (rescaled for origin x)
        arma::vec en_coef;      // (rescaled) elastic net estimates
        double negLogL;         // negative log-likelihood
        unsigned int coef_df;   // number of non-zero coef estimates

        // for a lambda sequence
        double alpha;           // tuning parameter
        arma::vec lambda_vec;   // lambda sequence
        arma::mat coef_mat;     // coef matrix (rescaled for origin x)
        arma::mat en_coef_mat;  // elastic net estimates
        arma::vec negLogL_vec;  // negative log-likelihood vector
        arma::uvec coef_df_vec; // coef df vector

        // default constructor
        LogisticReg() {}

        // constructors
        LogisticReg(const arma::mat& x_,
                    const arma::vec& y_,
                    const bool intercept_ = true,
                    const bool standardize_ = true)
        {
            intercept = intercept_;
            standardize = standardize_;
            x = x_;
            this->nObs = x.n_rows;
            if (standardize) {
                if (intercept) {
                    x_center = arma::mean(x);
                } else {
                    x_center = arma::zeros<arma::rowvec>(x.n_cols);
                }
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
                    for (size_t j {1}; j < coef0.n_elem; ++j) {
                        this->coef[j] = coef0[j] / x_scale[j - 1];
                    }
                } else {
                    for (size_t j {0}; j < coef0.n_elem; ++j) {
                        this->coef[j] = coef0[j] / x_scale[j];
                    }
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

        inline arma::vec linkinv(const arma::vec& eta,
                                 const double pmin) const;

        // here beta is coef vector for non-standardized data
        inline arma::vec predict(const arma::vec& beta) const;

        inline double objective() const;

        inline arma::vec gradient(const arma::vec& beta) const;

        inline double gradient(const arma::vec& beta,
                               const unsigned int k) const;

        // Firth-type score function
        inline arma::vec firth_score(const arma::vec& beta) const;

        inline arma::vec firth_score(const arma::vec& beta,
                                     const unsigned int k) const;

        inline double objective(const arma::vec& beta, arma::vec& grad) const;

        // compute iteration matrix in Bohning and Lindsay (1988)
        inline void compute_bl_iter_mat(bool force_update);

        // fit regular logistic regression model
        inline void fit(const arma::vec& start,
                        const unsigned int& max_iter,
                        const double& rel_tol);

        // fit firth logistic regression model
        inline void firth_fit(const arma::vec& start,
                              const unsigned int& max_iter,
                              const double& rel_tol);

        // compute cov lowerbound used in regularied model
        inline void compute_cmd_lowerbound(bool force_update);

        // update step for regularized logistic regression model
        inline void regularized_fit_update(arma::vec& beta,
                                           arma::uvec& is_active,
                                           const double& l1_lambda,
                                           const double& l2_lambda,
                                           const arma::vec& l1_penalty_factor,
                                           const bool& update_active);

        inline void reg_active_fit(arma::vec& beta,
                                   const arma::uvec& is_active,
                                   const double& l1_lambda,
                                   const double& l2_lambda,
                                   const arma::vec& l1_penalty_factor,
                                   const bool& varying_active_set,
                                   const unsigned int& max_iter,
                                   const double& rel_tol);

        // fit regularized logistic regression model
        // for a perticular lambda
        inline void regularized_fit(const double& l1_lambda,
                                    const double& l2_lambda,
                                    const arma::vec& l1_penalty_factor,
                                    const arma::vec& start,
                                    const unsigned int& max_iter,
                                    const double& rel_tol);

        // overload for a sequence of lambda's
        inline void regularized_fit(arma::vec lambda,
                                    const double& alpha,
                                    const unsigned int& nlambda,
                                    double lambda_min_ratio,
                                    const arma::vec& l1_penalty_factor,
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
    inline arma::vec LogisticReg::linkinv(const arma::vec& beta,
                                          const double pmin = 1e-5) const
    {
        arma::vec p_vec { 1 / (1 + arma::exp(- mat2vec(x * beta))) };
        // special care prevents coef diverging
        // reference: Friedman, J., Hastie, T., & Tibshirani, R. (2010)
        for (size_t i {0}; i < p_vec.n_elem; ++i) {
            if (p_vec(i) < pmin) {
                p_vec(i) = pmin;
            } else if (p_vec(i) > 1 - pmin) {
                p_vec(i) = 1 - pmin;
            }
        }
        return p_vec;
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
    inline double LogisticReg::objective() const
    {
        double res { 0 };
        arma::vec tmp { arma::zeros(2) };
        for (size_t i { 0 }; i < x.n_rows; ++i) {
            double x_beta { arma::as_scalar(x.row(i) * this->coef0) };
            tmp[1] = x_beta;
            res += log_sum_exp(tmp) - y(i) * x_beta;
        }
        return res;
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

    // Firth-type score function
    inline arma::vec LogisticReg::firth_score(const arma::vec& pi_vec) const
    {
        arma::mat b_mat {x};
        for (size_t i {0}; i < x.n_rows; ++i) {
            b_mat.row(i) *= std::sqrt(pi_vec(i) * (1 - pi_vec(i)));
        }
        // QR decomposition
        arma::mat q_mat, r_mat;
        arma::qr_econ(q_mat, r_mat, b_mat);
        // compute diagonal of hat matrix
        arma::vec hat_vec { arma::zeros(q_mat.n_rows) };
        for (size_t i {0}; i < q_mat.n_rows; ++i) {
            hat_vec(i) = arma::sum(arma::square(q_mat.row(i)));
        }
        // compute score function
        arma::vec res { arma::zeros(x.n_cols) };
        for (size_t i {0}; i < x.n_rows; ++i) {
            arma::rowvec tmp {
                (y(i) - pi_vec(i) + hat_vec(i) * (0.5 - pi_vec(i))) * x.row(i)
            };
            res += tmp.t();
        }
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

    // compute iteration matrix in Bohning and Lindsay (1988)
    inline void LogisticReg::compute_bl_iter_mat(bool force_update = false)
    {
        if (force_update || this->bl_iter_mat.is_empty()) {
            this->bl_iter_mat = 4 * arma::inv_sympd(x.t() * x);
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
        arma::mat iter_mat { this->bl_iter_mat * x.t() };
        size_t i {0};
        while (i < max_iter) {
            y_hat = linkinv(beta0);
            beta = beta0 + iter_mat * (y - y_hat);
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
        // compute negative log-likelihood
        this->negLogL = this->objective();
        this->coef_df = beta.n_elem;
    }

    // fitting Firth logistic model with monotonic quadratic approximation
    // algorithm;  non-integer y vector is allowed.
    inline void LogisticReg::firth_fit(const arma::vec& start = 0,
                                       const unsigned int& max_iter = 1000,
                                       const double& rel_tol = 1e-6)
    {
        arma::vec beta0 { arma::zeros(x.n_cols) };
        if (start.n_elem == x.n_cols) {
            beta0 = start;
        }
        arma::vec beta { beta0 }, score_vec { beta0 };
        arma::vec eta { arma::zeros(y.n_elem) };
        arma::vec y_hat { eta };
        this->compute_bl_iter_mat();
        size_t i {0};
        while (i < max_iter) {
            y_hat = linkinv(beta0);
            score_vec = firth_score(y_hat);
            beta = beta0 + this->bl_iter_mat * score_vec;
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
        // compute negative log-likelihood
        this->negLogL = this->objective();
        this->coef_df = beta.n_elem;
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
        const double& l1_lambda,
        const double& l2_lambda,
        const arma::vec& l1_penalty_factor,
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
                beta[j] = soft_threshold(cmd_lowerbound[j] * beta[j] - dlj,
                                         l1_penalty_factor[j] * l1_lambda) /
                    (cmd_lowerbound[j] +
                     2 * l2_lambda * static_cast<double>(j > 0));
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
        const double& l1_lambda,
        const double& l2_lambda,
        const arma::vec& l1_penalty_factor,
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
                    regularized_fit_update(beta, is_active_stored, l1_lambda,
                                           l2_lambda, l1_penalty_factor, true);
                    // double d_tol {
                    //     arma::max(cmd_lowerbound.t() %
                    //               arma::pow(beta - beta0, 2))
                    // };
                    // if (d_tol < rel_tol * rel_tol) {
                    //     break;
                    // }
                    if (rel_l2_norm(beta, beta0) < rel_tol) {
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                // run a full cycle over the converged beta
                regularized_fit_update(beta, is_active_new, l1_lambda,
                                       l2_lambda, l1_penalty_factor, true);
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
                regularized_fit_update(beta, is_active_stored, l1_lambda,
                                       l2_lambda, l1_penalty_factor, false);
                // double d_tol {
                //     arma::max(cmd_lowerbound.t() %
                //               arma::pow(beta - beta0, 2))
                // };
                // if (d_tol < rel_tol * rel_tol) {
                //     break;
                // }
                if (rel_l2_norm(beta, beta0) < rel_tol) {
                    break;
                }
                beta0 = beta;
                i++;
            }
        }
    }

    // regularized logistic model by coordinate-majorization-descent algorithm
    // for perticular lambda's for lasso penalty and ridge penalty
    // lambda_1 * factor * lasso + lambda_2 * ridge
    inline void LogisticReg::regularized_fit(
        const double& l1_lambda = 0,
        const double& l2_lambda = 0,
        const arma::vec& l1_penalty_factor = 0,
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
        arma::vec l1_penalty { arma::ones(n_predictor) };
        if (l1_penalty_factor.n_elem == l1_penalty.n_elem) {
            l1_penalty = l1_penalty_factor * x.n_cols /
                arma::sum(l1_penalty_factor);
        }
        this->l1_penalty_factor = l1_penalty;

        arma::vec beta { arma::zeros(x.n_cols) };
        arma::vec grad_zero { arma::abs(this->gradient(beta)) };
        arma::vec grad_beta { grad_zero }, strong_rhs { grad_beta };

        // large enough lambda for all-zero coef (except intercept)
        this->l1_lambda_max =
            arma::max(grad_zero.tail(n_predictor) /
                      l1_penalty) / this->x.n_rows;

        if (this->intercept) {
            l1_penalty = arma::join_vert(arma::zeros(1), l1_penalty);
        }

        // get solution of l1_lambda_max for a warm start
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x.n_cols) };
        if (this->intercept) {
            // only needs to estimate intercept
            is_active_strong(0) = 1;
            regularized_fit_update(beta, is_active_strong, this->l1_lambda_max,
                                   l2_lambda, l1_penalty, false);
        }
        this->coef0 = beta;
        // rescale coef back
        this->rescale_coef();

        // early exit for lambda greater than lambda_max
        if (l1_lambda >= this->l1_lambda_max) {
            this->en_coef = this->coef;
            this->coef_df = 0;
            return;
        }

        // use the input start if correctly specified
        if (start.n_elem == x.n_cols) {
            beta = start;
        }

        // update active set by strong rule
        grad_beta = arma::abs(this->gradient(this->coef0)) / this->x.n_rows;
        strong_rhs = (2 * l1_lambda - this->l1_lambda_max) * l1_penalty;

        for (size_t j { int_intercept }; j < n_predictor + int_intercept; ++j) {
            if (grad_beta(j) >= strong_rhs(j)) {
                is_active_strong(j) = 1;
            } else {
                beta(j) = 0;
            }
        }
        arma::uvec is_active_strong_new { is_active_strong };

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x.n_cols > x.n_rows || x.n_cols > 50) {
            varying_active_set = true;
        }

        bool kkt_failed { true };
        strong_rhs = l1_lambda * l1_penalty;
        // eventually, strong rule will guess correctly
        while (kkt_failed) {
            // update beta
            reg_active_fit(beta, is_active_strong, l1_lambda, l2_lambda,
                           l1_penalty, varying_active_set, max_iter, rel_tol);
            // check kkt condition
            for (size_t j { int_intercept };
                 j < n_predictor + int_intercept; ++j) {
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
        // compute elastic net estimates, then rescale them back
        this->coef0 = (1 + l2_lambda) * beta;
        this->rescale_coef();
        this->en_coef = this->coef;
        // overwrite the naive elastic net estimate
        this->coef0 = beta;
        this->rescale_coef();
        // compute negative log-likelihood
        this->negLogL = this->objective();
        // compute degree of freedom
        this->coef_df = get_coef_df(beta);
        // record other inputs
        this->l1_lambda = l1_lambda;
        this->l2_lambda = l2_lambda;
    }


    // for a sequence of lambda's
    // lambda * (penalty_factor * alpha * lasso + (1 - alpha) / 2 * ridge)
    inline void LogisticReg::regularized_fit(
        arma::vec lambda = 0,
        const double& alpha = 1,
        const unsigned int& nlambda = 1,
        double lambda_min_ratio = 1e-4,
        const arma::vec& l1_penalty_factor = 0,
        const unsigned int& max_iter = 1000,
        const double& rel_tol = 1e-6
        )
    {
        // check alpha
        if ((alpha < 0) || (alpha > 1)) {
            throw std::range_error("Alpha must be between 0 and 1.");
        }
        this->alpha = alpha;

        // set penalty terms
        unsigned int int_intercept { static_cast<unsigned int>(intercept) };
        unsigned int n_predictor { x.n_cols - int_intercept };
        if (n_predictor < 1) {
            throw std::range_error("Predictors not found.");
        }
        arma::vec l1_penalty { arma::ones(n_predictor) };
        if (l1_penalty_factor.n_elem == l1_penalty.n_elem) {
            l1_penalty = l1_penalty_factor * x.n_cols /
                arma::sum(l1_penalty_factor);
        }
        this->l1_penalty_factor = l1_penalty;

        // construct lambda sequence
        arma::vec beta { arma::zeros(x.n_cols) };
        arma::vec grad_beta { arma::abs(this->gradient(beta)) };
        arma::vec strong_rhs { beta };
        l1_lambda_max =
            arma::max(grad_beta.tail(l1_penalty.n_elem) / l1_penalty) /
            this->x.n_rows / std::max(alpha, 1e-10);

        // take unique lambda and sort descendingly
        lambda = arma::reverse(arma::unique(lambda));
        // construct lambda sequence
        arma::vec lambda_seq;
        if (nlambda > 1) {
            double log_lambda_max { std::log(this->l1_lambda_max) };
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
        this->lambda_vec = lambda_seq;

        // update penalty for intercept
        if (this->intercept) {
            l1_penalty = arma::join_vert(arma::zeros(1), l1_penalty);
        }

        // initialize the estimate matrix
        this->coef_mat = arma::zeros(x.n_cols, lambda_seq.n_elem);
        this->en_coef_mat = this->coef_mat;
        this->negLogL_vec = arma::zeros(lambda_seq.n_elem);
        this->coef_df_vec = arma::zeros<arma::uvec>(lambda_seq.n_elem);

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x.n_cols) };

        // get solution of lambda_max for a warm start
        if (this->intercept) {
            // only needs to estimate intercept
            is_active_strong(0) = 1;
            regularized_fit_update(
                beta, is_active_strong, this->l1_lambda_max * alpha,
                this->l1_lambda_max * (1 - alpha) / 2, l1_penalty, false
                );
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
            if (alpha * lambda_seq(k) >= this->l1_lambda_max) {
                this->coef_mat.col(k) = this->coef;
                continue;
            }
            // update acitve set by strong rule (for lambda < lamda_max)
            grad_beta = arma::abs(this->gradient(beta)) / this->x.n_rows;
            if (k == 0) {
                // use lambda_max
                strong_rhs = alpha *
                    (2 * lambda_seq(k) - this->l1_lambda_max) * l1_penalty;
            } else {
                // use the last lambda
                strong_rhs = alpha *
                    (2 * lambda_seq(k) - lambda_seq(k - 1)) * l1_penalty;
            }
            for (size_t j { int_intercept };
                 j < n_predictor + int_intercept; ++j) {
                if (grad_beta(j) > strong_rhs(j)) {
                    is_active_strong(j) = 1;
                } else {
                    beta(j) = 0;
                }
            }
            arma::uvec is_active_strong_new { is_active_strong };
            strong_rhs = alpha * lambda_seq(k) * l1_penalty;
            bool kkt_failed { true };
            // eventually, strong rule will guess correctly
            while (kkt_failed) {
                // update beta
                reg_active_fit(beta, is_active_strong, lambda_seq(k) * alpha,
                               lambda_seq(k) * (1 - alpha) / 2, l1_penalty,
                               varying_active_set, max_iter, rel_tol);
                // check kkt condition
                grad_beta = arma::abs(this->gradient(beta)) / this->x.n_rows;
                for (size_t j { int_intercept };
                     j < n_predictor + int_intercept; ++j) {
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
            // compute elastic net estimates
            this->coef0 = (1 + (1 - alpha) * lambda_seq(k) / 2) * beta;
            this->rescale_coef();
            this->en_coef_mat.col(k) = this->coef;
            // compute naive elastic net estimates
            this->coef0 = beta;
            this->rescale_coef();
            this->coef_mat.col(k) = this->coef;
            // compute negative log-likelihood
            this->negLogL_vec(k) = this->objective();
            // compute degree of freedom
            this->coef_df_vec(k) = get_coef_df(beta);
        }
        // prepare outputs

    }

}


#endif
