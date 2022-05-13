//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2021  Wenjie Wang <wang@wwenjie.org>
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

#ifndef INTSURV_LOGISTIC_REG_H
#define INTSURV_LOGISTIC_REG_H

#include <utility>
#include <RcppArmadillo.h>
#include "control.h"
#include "utils.h"

namespace Intsurv {

    class LogisticReg
    {
    protected:
        // internals
        double dn_obs_;         // double version of n_obs
        arma::rowvec x_center_; // the column center of x
        arma::rowvec x_scale_;  // the scale of x
        unsigned int int_intercept_;
        unsigned int p0_;      // number of covariates without intercept
        // iteration matrix in Bohning and Lindsay (1988)
        arma::mat bl_iter_mat_;
        // for regularized coordinate majorization descent
        arma::rowvec cmd_lowerbound_;

        // process starting values
        inline arma::vec gen_start(const arma::vec& start = arma::vec()) const
        {
            if (start.n_elem == p_) {
                return rev_rescale_coef(start);
            }
            return arma::zeros(p_);
        }

        // process offset
        inline arma::vec gen_offset(const arma::vec& offset = arma::vec()) const
        {
            if (offset.n_elem == n_obs_) {
                return offset;
            }
            return arma::zeros(n_obs_);
            // if (offset.n_elem == 1 || offset.empty()) {
            //     return arma::zeros(n_obs_);
            // }
            // throw std::length_error(
            //     "The length of the specified offset must match sample size."
            //     );
        }

        // process penalty factor
        inline arma::vec gen_penalty_factor(
            const arma::vec& penalty_factor = arma::vec()
            ) const
        {
            if (penalty_factor.n_elem == p_) {
                if (arma::any(penalty_factor < 0.0)) {
                    throw std::range_error(
                        "The 'penalty_factor' cannot be negative.");
                }
                return penalty_factor;
            }
            arma::vec out { arma::ones(p_) };
            out[0] = 0.0;
            if (penalty_factor.n_elem == p0_) {
                for (size_t j {1}; j < p_; ++j) {
                    out[j] = penalty_factor[j - 1];
                }
                return out;
            }
            // else
            return out;
        }

        // make sure the following are well set
        inline LogisticReg* check_start()
        {
            if (control_.start_.n_elem != p_) {
                control_.start_ = gen_start(control_.start_);
            }
            return this;
        }

        inline LogisticReg* check_offset()
        {
            if (control_.offset_.n_elem != n_obs_) {
                control_.offset_ = gen_offset(control_.offset_);
            }
            return this;
        }

        inline LogisticReg* check_penalty_factor()
        {
            if (control_.penalty_factor_.n_elem != p_) {
                control_.penalty_factor_ = gen_penalty_factor(
                    control_.penalty_factor_);
            }
            return this;
        }

        //! @param beta coef estimates for the standardized x
        inline arma::vec predict0(const arma::vec& beta) const
        {
            arma::vec p_vec {
                1.0 / (1.0 + arma::exp(- mat2vec(x_ * beta) - control_.offset_))
            };
            // special care prevents coef diverging
            // reference: Friedman, J., Hastie, T., & Tibshirani, R. (2010)
            set_pmin_bound(p_vec, control_.pmin_);
            return p_vec;
        }

        inline double objective0(const arma::vec& beta) const
        {
            double res { 0.0 };
            for (size_t i { 0 }; i < x_.n_rows; ++i) {
                double x_beta {
                    arma::as_scalar(x_.row(i) * beta + control_.offset_(i))
                };
                res += std::log(std::exp(x_beta) + 1.0) - y_(i) * x_beta;
            }
            return res / dn_obs_;
        }

        inline arma::vec gradient0(const arma::vec& beta) const
        {
            return x_.t() * (predict0(beta) - y_) / dn_obs_;
        }

        // define gradient0 function at k-th dimension
        inline double gradient0(const arma::vec& beta,
                                const unsigned int k) const
        {
            return arma::accu((predict0(beta) - y_) % x_.col(k)) / dn_obs_;
        }

        // define objective function and overwrites graidient
        inline double objective0(const arma::vec& beta,
                                 arma::vec& grad) const
        {
            arma::vec x_beta {x_ * beta + control_.offset_};
            arma::vec exp_x_beta {arma::exp(x_beta)};
            grad = x_.t() * (exp_x_beta / (1.0 + exp_x_beta) - y_);
            arma::vec y_x_beta {y_ % x_beta};
            double res {
                arma::as_scalar(
                    arma::mean(arma::log(1.0 + exp_x_beta) - y_x_beta)
                    )
            };
            return res;
        }

        // compute iteration matrix in Bohning and Lindsay (1988)
        inline void set_bl_iter_mat(const bool force_update = false)
        {
            if (force_update || bl_iter_mat_.is_empty()) {
                bl_iter_mat_ = 4 * arma::inv_sympd(x_.t() * x_);
            }
        }

        // compute cov lowerbound used in regularied model
        inline void set_cmd_lowerbound(const bool force_update = false)
        {
            if (force_update || cmd_lowerbound_.is_empty()) {
                cmd_lowerbound_ = arma::sum(arma::square(x_), 0) /
                    (4.0 * dn_obs_);
            }
        }

        // update step for regularized logistic regression model
        inline void net_one_update(arma::vec& beta,
                                   arma::uvec& is_active,
                                   const double l1_lambda,
                                   const double l2_lambda,
                                   const arma::vec& penalty_factor,
                                   const bool update_active,
                                   const unsigned int verbose);

        inline void net_active_update(arma::vec& beta,
                                      arma::uvec& is_active,
                                      const double l1_lambda,
                                      const double l2_lambda,
                                      const arma::vec& penalty_factor,
                                      const bool varying_active,
                                      const unsigned int max_iter,
                                      const double epsilon,
                                      const unsigned int verbose);

    public:
        // model =============================================================
        arma::mat x_;           // (standardized) x
        arma::vec y_;
        unsigned int n_obs_;    // number of observations
        unsigned int p_;        // number of covariates with possible intercept
        // for a soltuon path
        double l1_lambda_max_;  // the "big enough" l1 lambda => zero coef
        double lambda_max_;     // l1_lambda_max / alpha
        // outputs ===========================================================
        arma::vec coef0_;       // coef before rescaling
        // for a single set of l1_lambda and l2_lambda
        arma::vec coef_;        // coef (rescaled for origin x)
        // for a lambda sequence
        arma::mat coef_mat_;    // coef matrix (rescaled for origin x)

        // controls
        Control control_;

        // default constructor
        LogisticReg() {}

        // constructors
        LogisticReg(const arma::mat& x,
                    const arma::vec& y,
                    const Control& control = Control()) :
            x_ (x),
            y_ (y),
            control_ (control)
        {
            int_intercept_ = static_cast<unsigned int>(control_.intercept_);
            n_obs_ = x_.n_rows;
            p0_ = x_.n_cols;
            p_ = p0_ + int_intercept_;
            dn_obs_ = static_cast<double>(n_obs_);
            if (control_.standardize_) {
                if (control_.intercept_) {
                    x_center_ = arma::mean(x_);
                } else {
                    x_center_ = arma::zeros<arma::rowvec>(p0_);
                }
                x_scale_ = arma::stddev(x_, 1);
                for (size_t j {0}; j < p0_; ++j) {
                    if (x_scale_(j) > 0) {
                        x_.col(j) = (x_.col(j) - x_center_(j)) / x_scale_(j);
                    } else {
                        // coef will be zero, set non-zero for rescaling
                        x_.col(j) = arma::zeros(n_obs_);
                        x_scale_(j) = - 1.0;
                    }
                }
            }
            if (control_.intercept_) {
                x_ = arma::join_horiz(arma::ones(n_obs_), x_);
            }
            set_offset();
        }

        inline LogisticReg* set_start(const arma::vec& start)
        {
            control_.start_ = gen_start(start);
            return this;
        }
        inline LogisticReg* set_start()
        {
            return set_start(control_.start_);
        }
        // set offset
        inline LogisticReg* set_offset(const arma::vec& offset)
        {
            control_.offset_ = gen_offset(offset);
            return this;
        }
        inline LogisticReg* set_offset()
        {
            return set_offset(control_.offset_);
        }
        // reset offset to zeros
        inline LogisticReg* reset_offset()
        {
            control_.offset_ = arma::zeros(n_obs_);
            return this;
        }
        inline LogisticReg* set_penalty_factor(const arma::vec& penalty_factor)
        {
            control_.penalty_factor_ = gen_penalty_factor(penalty_factor);
            return this;
        }
        inline LogisticReg* set_penalty_factor()
        {
            return set_penalty_factor(control_.penalty_factor_);
        }

        // transform coef for standardized data to the one for original data
        inline void rescale_coef()
        {
            coef_ = coef0_;
            if (control_.standardize_) {
                if (control_.intercept_) {
                    arma::uvec non_int_ind {
                        arma::regspace<arma::uvec>(1, p0_)
                    };
                    coef_[0] = coef0_(0) -
                        arma::as_scalar((x_center_ / x_scale_) *
                                        coef0_.elem(non_int_ind));
                    for (size_t j {1}; j < p_; ++j) {
                        coef_[j] = coef0_[j] / x_scale_[j - 1];
                    }
                } else {
                    for (size_t j {0}; j < p0_; ++j) {
                        coef_[j] = coef0_[j] / x_scale_[j];
                    }
                }
            }
        }
        // transform coef for original data to the one for standardized data
        inline arma::vec rev_rescale_coef(const arma::vec& beta) const
        {
            if (control_.standardize_) {
                arma::vec beta0 { beta };
                double tmp {0};
                for (size_t j {1}; j < beta.n_elem; ++j) {
                    beta0(j) *= x_scale_(j - 1);
                    tmp += beta(j) * x_center_(j - 1);
                }
                beta0(0) += tmp;
                return beta0;
            }
            return beta;
        }

        inline arma::vec predict() const
        {
            return predict0(coef0_);
        }
        //! @param beta coef vector for original x
        inline arma::vec predict(const arma::vec& beta) const
        {
            arma::vec beta0 { beta };
            if (control_.standardize_) {
                beta0 = rev_rescale_coef(beta0);
            }
            return predict0(beta0);
        }
        // define objective function (negative log-likehood function)
        inline double objective() const
        {
            return objective0(coef0_);
        }
        //! @param beta coef estimates for the original x
        inline double objective(const arma::vec& beta) const
        {
            arma::vec beta0 { beta };
            if (control_.standardize_) {
                beta0 = rev_rescale_coef(beta0);
            }
            return objective0(beta0);
        }
        inline double net_penalty(const arma::vec& beta,
                                  const double l1_lambda,
                                  const double l2_lambda,
                                  const arma::vec& penalty_factor) const
        {
            if (control_.intercept_) {
                arma::mat beta0int { beta.tail_rows(p0_) };
                return l1_lambda * l1_norm(beta % penalty_factor) +
                    l2_lambda * sum_of_square(beta0int);
            }
            return l1_lambda * l1_norm(beta) +
                l2_lambda * sum_of_square(beta);
        }
        inline double net_penalty() const
        {
            return net_penalty(coef0_,
                               control_.l1_lambda_,
                               control_.l2_lambda_,
                               control_.penalty_factor_);
        }
        inline double get_l1_lambda_max(const arma::vec& penalty_factor) const
        {
            arma::uvec active_penalty { arma::find(penalty_factor > 0.0) };
            arma::uvec penalty_free { arma::find(penalty_factor == 0.0) };
            arma::vec beta { arma::zeros(p_) };
            arma::vec grad_beta { arma::abs(gradient0(beta)) };
            double l1_lambda_max { 0.0 };
            for (arma::uvec::iterator it { active_penalty.begin() };
                 it != active_penalty.end(); ++it) {
                double tmp { grad_beta(*it) };
                tmp /= control_.penalty_factor_(*it);
                if (l1_lambda_max < tmp) {
                    l1_lambda_max = tmp;
                }
            }
            return l1_lambda_max;
        }
        inline LogisticReg* set_l1_lambda_max()
        {
            l1_lambda_max_ = get_l1_lambda_max(control_.penalty_factor_);
            return this;
        }

        // fit regular logistic regression model
        inline void fit();

        // fit regularized logistic regression model
        // for a perticular lambda
        inline void net_fit();

        // for a sequence of lambda's
        inline void net_path();

        // function that helps update y_
        inline void update_y(const arma::vec& y) { y_ = y; }

        // getters
        inline arma::mat get_x(const bool rescale,
                               const bool with_intercept) const
        {
            arma::mat out {x_};
            if (rescale && control_.standardize_) {
                for (size_t j { 0 }; j < p0_; ++j) {
                    out.col(j + int_intercept_) =
                        x_scale_(j) * out.col(j + int_intercept_) +
                        x_center_(j);
                }
            }
            if (! with_intercept && control_.intercept_) {
                out.shed_col(0);
            }
            return out;
        }
        inline arma::vec get_xbeta() const
        {
            return mat2vec(x_ * coef0_);
        }

    };

    // fitting regular logistic model by monotonic quadratic approximation
    // algorithm non-integer y vector is allowed
    // reference: Bohning and Lindsay (1988) SIAM
    inline void LogisticReg::fit()
    {
        set_bl_iter_mat();
        check_start()->check_offset();
        arma::vec beta0 { control_.start_ };
        double ell { arma::datum::inf };
        if (control_.verbose_ > 1) {
            Rcpp::Rcout << "\n" << std::string(40, '=')
                        << "\nStarting from\n"
                        << arma2rvec(beta0)
                        << "\n";
        }
        if (control_.verbose_ > 0) {
            ell = objective(beta0);
        }
        arma::vec beta { beta0 };
        arma::vec y_hat;
        arma::mat iter_mat { bl_iter_mat_ * x_.t() };
        // main loop
        for (size_t i {0}; i < control_.max_iter_; ++i) {
            y_hat = predict0(beta0);
            beta = beta0 + iter_mat * (y_ - y_hat);
            if (control_.verbose_ > 1) {
                Rcpp::Rcout << "\n"
                            << std::string(40, '=')
                            << "\nitartion: "
                            << i + 1
                            << "\n  coef estimates: "
                            << arma2rvec(beta)
                            << "\n";
            }
            if (control_.verbose_ > 0) {
                double ell_old { ell };
                ell = objective0(beta);
                Rcpp::Rcout << "\n  The negative log-likelihood changed\n";
                Rprintf("  from %15.15f\n", ell_old);
                Rprintf("    to %15.15f\n", ell);
                if (ell_old < ell) {
                    Rcpp::Rcout << "Warning: The negative log-likelihood"
                                << " somehow increased\n";
                }
            }
            // if relative tolerance is statisfied
            if (rel_l1_norm(beta, beta0) < control_.epsilon_) {
                if (control_.verbose_ > 0) {
                    Rcpp::Rcout << "\nReached convergence criterion\n";
                }
                break;
            }
            // update beta
            beta0 = beta;
        }
        coef0_ = std::move(beta);
        rescale_coef();
    }

    // run one cycle of coordinate descent over a given active set
    inline void LogisticReg::net_one_update(
        arma::vec& beta,
        arma::uvec& is_active,
        const double l1_lambda,
        const double l2_lambda,
        const arma::vec& penalty_factor,
        const bool update_active,
        const unsigned int verbose
        )
    {
        double ell_verbose { 0.0 }, obj_verbose { 0.0 }, reg_verbose { 0.0 };
        if (verbose > 2) {
            Rcpp::Rcout << "\nStarting values of beta:\n"
                        << beta
                        << "\nThe active set of beta:\n"
                        << arma2rvec(is_active)
                        << "\n";
        }
        if (verbose > 1) {
            obj_verbose = objective0(beta);
            reg_verbose = net_penalty(beta, l1_lambda, l2_lambda,
                                      penalty_factor);
            ell_verbose = obj_verbose + reg_verbose;
        }
        arma::vec beta_old { beta };
        for (size_t j {0}; j < p_; ++j) {
            if (is_active[j] == 0) {
                continue;
            }
            double dlj { gradient0(beta, j) };
            // if cmd_lowerbound_ = 0 and l1_lambda > 0, numer will be 0
            double numer {
                soft_threshold(cmd_lowerbound_[j] * beta[j] - dlj,
                               penalty_factor[j] * l1_lambda)
            };
            if (isAlmostEqual(numer, 0)) {
                beta[j] = 0;
            } else {
                double denom {
                    cmd_lowerbound_[j] + 2 * l2_lambda *
                    static_cast<double>(j >= int_intercept_)
                };
                beta[j] = numer / denom;
            }
            if (update_active) {
                // check if it has been shrinkaged to zero
                if (isAlmostEqual(beta[j], 0)) {
                    is_active[j] = 0;
                } else {
                    is_active[j] = 1;
                }
            }
        }
        if (verbose > 1) {
            double ell_old { ell_verbose };
            Rcpp::Rcout << "The objective function changed\n";
            Rprintf("  from %7.7f (obj. %7.7f + reg. %7.7f)\n",
                    ell_verbose, obj_verbose, reg_verbose);
            obj_verbose = objective0(beta);
            reg_verbose = net_penalty(beta, l1_lambda, l2_lambda,
                                      penalty_factor);
            ell_verbose = obj_verbose + reg_verbose;
            Rprintf("    to %7.7f (obj. %7.7f + reg. %7.7f)\n",
                    ell_verbose, obj_verbose, reg_verbose);
            if (ell_verbose > ell_old) {
                Rcpp::Rcout << "Warning: "
                            << "the objective function somehow increased\n";
            }
        }
    }

    // run a complete cycle of CMD for a given active set and lambda
    inline void LogisticReg::net_active_update(
        arma::vec& beta,
        arma::uvec& is_active,
        const double l1_lambda,
        const double l2_lambda,
        const arma::vec& penalty_factor,
        const bool varying_active,
        const unsigned int max_iter,
        const double epsilon,
        const unsigned int verbose
        )
    {
        unsigned int num_iter {0};
        arma::vec beta0 { beta };
        if (varying_active) {
            arma::uvec is_active_strong { is_active },
                is_active_varying { is_active };
            if (verbose > 1) {
                Rcpp::Rcout << "The size of active set from strong rule: "
                            << l1_norm(is_active_strong)
                            << "\n";
            }
            for (size_t i {0}; i < max_iter; ++i) {
                // cycles over the active set
                size_t ii {0};
                while (ii < max_iter) {
                    net_one_update(beta, is_active_varying, l1_lambda,
                                   l2_lambda, penalty_factor, true, verbose);
                    if (rel_l1_norm(beta, beta0) < epsilon) {
                        num_iter = ii + 1;
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                // run a full cycle over the converged beta
                net_one_update(beta, is_active, l1_lambda,
                               l2_lambda, penalty_factor, true, verbose);
                // check if two active sets coincide
                if (l1_norm(is_active_varying - is_active) > 0) {
                    // if different, repeat this process
                    if (verbose > 1) {
                        Rcpp::Rcout << "Changed the active set from "
                                    << l1_norm(is_active_varying)
                                    << " to "
                                    << l1_norm(is_active)
                                    << " after "
                                    << num_iter + 1
                                    << " iteration(s)\n";
                    }
                    is_active_varying = is_active;
                    // recover the active set
                    is_active = is_active_strong;
                } else {
                    if (verbose > 1) {
                        Rcpp::Rcout << "Converged over the active set after "
                                    << num_iter + 1
                                    << " iteration(s)\n";
                        Rcpp::Rcout << "The size of active set is "
                                    << l1_norm(is_active) << "\n";
                    }
                    num_iter = i + 1;
                    break;
                }
            }
        } else {
            // regular coordinate descent
            for (size_t i {0}; i < max_iter; ++i) {
                net_one_update(beta,
                               is_active,
                               l1_lambda,
                               l2_lambda,
                               penalty_factor,
                               false,
                               verbose);
                if (rel_l1_norm(beta, beta0) < epsilon) {
                    num_iter = i + 1;
                    break;
                }
                beta0 = beta;
            }
        }
        if (verbose > 0) {
            if (num_iter < max_iter) {
                Rcpp::Rcout << "Converged after "
                            << num_iter
                            << " iteration(s)\n";
            } else {
                msg("Reached the maximum number of iteratons.");
            }
        }
    }

    // regularized logistic model by coordinate-majorization-descent algorithm
    // for particular lambda's for lasso penalty and ridge penalty
    // lambda_1 * factor * lasso + lambda_2 * ridge
    inline void LogisticReg::net_fit()
    {
        set_cmd_lowerbound();
        check_start()->check_offset()->check_penalty_factor();
        // use the given starting values
        arma::vec beta { gen_start(control_.start_) };
        arma::uvec is_active { arma::ones<arma::uvec>(p_) };
        net_active_update(beta,
                          is_active,
                          control_.l1_lambda_,
                          control_.l2_lambda_,
                          control_.penalty_factor_,
                          control_.varying_active_,
                          control_.max_iter_,
                          control_.epsilon_,
                          control_.verbose_);
        coef0_ = beta;
        rescale_coef();
    }


    // for a sequence of lambda's
    // lambda * (penalty_factor * alpha_ * lasso + (1 - alpha_) / 2 * ridge)
    inline void LogisticReg::net_path()
    {
        set_cmd_lowerbound();
        check_start()->check_offset()->check_penalty_factor();
        const bool is_ridge_only { isAlmostEqual(control_.alpha_, 0.0) };
        arma::uvec active_penalty {
            arma::find(control_.penalty_factor_ > 0.0)
        };
        arma::uvec penalty_free { arma::find(control_.penalty_factor_ == 0.0) };
        // construct lambda sequence
        arma::vec beta { arma::zeros(x_.n_cols) };
        arma::vec grad_beta, strong_rhs;
        // if alpha = 0 and lambda is specified
        if (is_ridge_only && ! control_.lambda_.empty()) {
            control_.lambda_ = arma::reverse(arma::unique(control_.lambda_));
            l1_lambda_max_ = 0.0;    // not well defined
            lambda_max_ = 0.0;       // not well defined
        } else {
            // need to determine l1_lambda_max
            set_l1_lambda_max();
            lambda_max_ = l1_lambda_max_ / std::max(control_.alpha_, 1e-2);
            // set up lambda sequence
            if (control_.lambda_.empty()) {
                double log_lambda_max { std::log(lambda_max_) };
                control_.lambda_ = arma::exp(
                    arma::linspace(log_lambda_max,
                                   log_lambda_max +
                                   std::log(control_.lambda_min_ratio_),
                                   control_.nlambda_)
                    );
            } else {
                control_.lambda_ = arma::reverse(
                    arma::unique(control_.lambda_));
            }
        }
        // initialize the estimate matrix
        coef_mat_ = arma::zeros(p_, control_.lambda_.n_elem);
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(p_) };
        // for ridge penalty
        if (is_ridge_only) {
            is_active_strong = arma::ones<arma::uvec>(x_.n_cols);
            for (size_t li { 0 }; li < control_.lambda_.n_elem; ++li) {
                net_active_update(beta,
                                  is_active_strong,
                                  0.0,
                                  0.5 * control_.lambda_(li),
                                  control_.penalty_factor_,
                                  false,
                                  control_.max_iter_,
                                  control_.epsilon_,
                                  control_.verbose_);
                coef0_ = beta;
                rescale_coef();
                coef_mat_.col(li) = coef_;
            }
            return;             // early exit
        }
        // get solution of lambda_max for a warm start
        for (arma::uvec::iterator it { penalty_free.begin() };
             it != penalty_free.end(); ++it) {
            is_active_strong(*it) = 1;
        }
        double l1_lambda { lambda_max_ * control_.alpha_ };
        double l2_lambda { 0.5 * lambda_max_ * (1 - control_.alpha_) };
        net_active_update(beta,
                          is_active_strong,
                          l1_lambda,
                          l2_lambda,
                          control_.penalty_factor_,
                          false,
                          control_.max_iter_,
                          control_.epsilon_,
                          control_.verbose_);
        double old_l1_lambda { l1_lambda_max_ }; // for strong rule
        // outer loop for the lambda sequence
        for (size_t k {0}; k < control_.lambda_.n_elem; ++k) {
            double lambda_k { control_.lambda_(k) };
            l1_lambda = lambda_k * control_.alpha_;
            l2_lambda = 0.5 * lambda_k * (1 - control_.alpha_);
            // early exit for large lambda greater than lambda_max
            if (l1_lambda >= l1_lambda_max_) {
                coef0_ = beta;
                rescale_coef();
                coef_mat_.col(k) = coef_;
                continue;
            }
            // update acitve set by strong rule (for lambda < lamda_max)
            grad_beta = arma::abs(gradient0(beta));
            strong_rhs = (2 * l1_lambda - old_l1_lambda) *
                control_.penalty_factor_;
            for (arma::uvec::iterator it { active_penalty.begin() };
                 it != active_penalty.end(); ++it) {
                if (is_active_strong(*it) > 0) {
                    continue;
                }
                if (grad_beta(*it) >= strong_rhs(*it)) {
                    is_active_strong(*it) = 1;
                }
            }
            arma::uvec is_active_strong_old { is_active_strong };
            strong_rhs = l1_lambda * control_.penalty_factor_;
            bool kkt_failed { true };
            // eventually, strong rule will guess correctly
            while (kkt_failed) {
                arma::uvec is_strong_rule_failed {
                    arma::zeros<arma::uvec>(p_)
                };
                // update beta
                net_active_update(beta,
                                  is_active_strong,
                                  l1_lambda,
                                  l2_lambda,
                                  control_.penalty_factor_,
                                  control_.varying_active_,
                                  control_.max_iter_,
                                  control_.epsilon_,
                                  control_.verbose_);
                // check kkt condition
                if (control_.verbose_ > 0) {
                    msg("Checking the KKT condition for the null set.");
                }
                for (arma::uvec::iterator it { active_penalty.begin() };
                     it != active_penalty.end(); ++it) {
                    if (is_active_strong_old(*it)) {
                        continue;
                    }
                    if (std::abs(gradient0(beta, *it)) > strong_rhs(*it)) {
                        is_strong_rule_failed(*it) = 1;
                    }
                }
                if (arma::accu(is_strong_rule_failed) > 0) {
                    is_active_strong = is_active_strong_old ||
                        is_strong_rule_failed;
                    if (control_.verbose_ > 0) {
                        Rcpp::Rcout << "The strong rule failed for "
                                    << arma::accu(is_strong_rule_failed)
                                    << " group(s)\nThe size of old active set: "
                                    << l1_norm(is_active_strong_old)
                                    << "\nThe size of new active set: "
                                    << l1_norm(is_active_strong)
                                    << "\n";
                    }
                } else {
                    if (control_.verbose_ > 0) {
                        msg("The strong rule worked.\n");
                    }
                    kkt_failed = false;
                }
            }
            coef0_ = beta;
            rescale_coef();
            coef_mat_.col(k) = coef_;
        }
    }

}


#endif
