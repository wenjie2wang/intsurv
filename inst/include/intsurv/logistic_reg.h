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

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {

    // define class for inputs and outputs
    class LogisticReg
    {
    public:
        // model data ========================================================
        arma::mat x_;           // (standardized) x
        arma::vec y_;
        bool intercept_;
        arma::vec offset_;      // offset term

        // internals =========================================================
        unsigned int n_obs_;    // number of observations
        double dn_obs_;         // double version of n_obs
        arma::rowvec x_center_; // the column center of x
        arma::rowvec x_scale_;  // the scale of x
        unsigned int int_intercept_;
        unsigned int p0_;      // number of covariates without intercept
        unsigned int p1_;      // number of covariates with possible intercept
        // iteration matrix in Bohning and Lindsay (1988)
        arma::mat bl_iter_mat_;
        // for regularized coordinate majorization descent
        arma::rowvec cmd_lowerbound_;

        // controls ==========================================================
        bool standardize_;            // is x standardized
        arma::vec l1_penalty_factor_; // adaptive weights for lasso penalty
        double l1_lambda_max_;        // the "big enough" lambda => zero coef
        // for a sinle l1_lambda and l2_lambda
        double l1_lambda_;      // tuning parameter for lasso penalty
        double l2_lambda_;      // tuning parameter for ridge penalty
        double alpha_;          // tuning parameter
        arma::vec lambda_vec_;  // lambda sequence

        // outputs ===========================================================
        arma::vec coef0_;       // coef before rescaling
        // for one solution
        arma::vec coef_;        // coef (rescaled for origin x)
        arma::vec en_coef_;     // (rescaled) elastic net estimates
        arma::vec xbeta_;       // sorted x * coef
        arma::vec prob_vec_;    // sorted linkinv response
        double neg_ll_;         // negative log-likelihood
        unsigned int coef_df_;  // number of non-zero coef estimates
        // for a lambda sequence
        arma::mat coef_mat_;     // coef matrix (rescaled for origin x)
        arma::mat en_coef_mat_;  // elastic net estimates
        arma::vec neg_ll_vec_;   // negative log-likelihood vector
        arma::uvec coef_df_vec_; // coef df vector

        // default constructor
        LogisticReg() {}

        // constructors
        LogisticReg(const arma::mat& x,
                    const arma::vec& y,
                    const bool intercept = true,
                    const bool standardize = true) :
            x_ (x),
            y_ (y),
            intercept_ (intercept),
            standardize_ (standardize)
        {
            int_intercept_ = static_cast<unsigned int>(intercept);
            n_obs_ = x_.n_rows;
            p0_ = x_.n_cols;
            p1_ = p0_ + int_intercept_;
            dn_obs_ = static_cast<double>(n_obs_);
            if (standardize_) {
                if (intercept_) {
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
            if (intercept_) {
                x_ = arma::join_horiz(arma::ones(n_obs_), x_);
            }
            // initialize offset_
            reset_offset();
        }

        // function members
        // set offset
        inline void set_offset(const arma::vec& offset)
        {
            if (offset.n_elem == n_obs_) {
                offset_ = offset;
            } else {
                throw std::length_error(
                    "The length of the specified offset must match sample size."
                    );
            }
        }
        // reset offset to zeros
        inline void reset_offset()
        {
            offset_ = arma::zeros(n_obs_);
        }

        // transfer coef for standardized data to coef for non-standardized data
        inline void rescale_coef()
        {
            coef_ = coef0_;
            if (standardize_) {
                if (intercept_) {
                    arma::uvec non_int_ind {
                        arma::regspace<arma::uvec>(1, p0_)
                    };
                    coef_[0] = coef0_(0) -
                        arma::as_scalar((x_center_ / x_scale_) *
                                        coef0_.elem(non_int_ind));
                    for (size_t j {1}; j < p1_; ++j) {
                        coef_[j] = coef0_[j] / x_scale_[j - 1];
                    }
                } else {
                    for (size_t j {0}; j < p0_; ++j) {
                        coef_[j] = coef0_[j] / x_scale_[j];
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
                beta0(j) *= x_scale_(j - 1);
                tmp += beta(j) * x_center_(j - 1);
            }
            beta0(0) += tmp;
            return beta0;
        }

        // additional methods for coxph_cure
        // update coef0, en_coef, and coef_df from a new coef
        inline void update_from_coef(const double l2_lambda = 0)
        {
            // update coef0_
            coef0_ = rev_rescale_coef(coef_);
            arma::vec beta { coef0_ };
            // update en_coef_
            if (l2_lambda > 0) {
                coef0_ *= (1 + l2_lambda);
                rescale_coef();
                en_coef_ = coef_;
                // overwrite the naive elastic net estimate
                coef0_ = beta;
                rescale_coef();
            } else {
                en_coef_ = coef_;
            }
            coef_df_ = get_coef_df(beta);
        }

        inline arma::vec linkinv(const arma::vec& eta,
                                 const double pmin) const;

        // here beta is coef_ vector for non-standardized data
        inline arma::vec predict(const arma::vec& beta,
                                 const double pmin) const;

        inline double objective() const;
        inline double objective(const arma::vec& beta) const;

        inline arma::vec gradient(const arma::vec& beta,
                                  const double pmin) const;
        inline double gradient(const arma::vec& beta,
                               const unsigned int k,
                               const double pmin) const;

        inline double objective(const arma::vec& beta, arma::vec& grad) const;

        // Firth-type score function
        inline arma::vec firth_score(const arma::vec& beta) const;
        inline arma::vec firth_score(const arma::vec& beta,
                                     const unsigned int k) const;

        // compute iteration matrix in Bohning and Lindsay (1988)
        inline void compute_bl_iter_mat(const bool force_update);

        // fit regular logistic regression model
        inline void fit(const arma::vec& start,
                        const unsigned int max_iter,
                        const double rel_tol,
                        const double pmin,
                        const bool early_stop,
                        const bool verbose);

        // fit firth logistic regression model
        inline void firth_fit(const arma::vec& start,
                              const unsigned int max_iter,
                              const double rel_tol,
                              const double pmin);

        // compute cov lowerbound used in regularied model
        inline void set_cmd_lowerbound();

        // update step for regularized logistic regression model
        inline void regularized_fit_update(arma::vec& beta,
                                           arma::uvec& is_active,
                                           const double l1_lambda,
                                           const double l2_lambda,
                                           const arma::vec& l1_penalty_factor,
                                           const bool update_active,
                                           const double pmin,
                                           const bool early_stop,
                                           const bool verbose);

        inline void reg_active_fit(arma::vec& beta,
                                   const arma::uvec& is_active,
                                   const double l1_lambda,
                                   const double l2_lambda,
                                   const arma::vec& l1_penalty_factor,
                                   const bool varying_active_set,
                                   const unsigned int max_iter,
                                   const double rel_tol,
                                   const double pmin,
                                   const bool early_stop,
                                   const bool verbose);

        // fit regularized logistic regression model
        // for a perticular lambda
        inline void regularized_fit(const double l1_lambda,
                                    const double l2_lambda,
                                    const arma::vec& l1_penalty_factor,
                                    const arma::vec& start,
                                    const unsigned int max_iter,
                                    const double rel_tol,
                                    const double pmin,
                                    const bool early_stop,
                                    const bool verbose);

        // overload for a sequence of lambda's
        inline void regularized_fit(arma::vec lambda,
                                    const double alpha_,
                                    const unsigned int nlambda,
                                    double lambda_min_ratio,
                                    const arma::vec& l1_penalty_factor,
                                    const unsigned int max_iter,
                                    const double rel_tol,
                                    const double pmin,
                                    const bool early_stop,
                                    const bool verbose);

        // function that helps update y_
        inline void update_y(const arma::vec& y) { y_ = y; }

        // some simple functions
        inline unsigned int sample_size() const
        {
            return n_obs_;
        }

        // getters
        inline arma::mat get_x(const bool& include_intercept = true) const
        {
            arma::mat out {x_};
            if (include_intercept && intercept_) {
                out.shed_col(0);
            }
            if (standardize_) {
                for (size_t j {0}; j < out.n_cols; ++j) {
                    out.col(j) = x_scale_(j) * out.col(j) + x_center_(j);
                }
            }
            return out;
        }
        inline arma::vec get_y() const { return y_; }

    };

    // define inverse link function
    inline arma::vec LogisticReg::linkinv(const arma::vec& beta,
                                          const double pmin = 1e-5) const
    {
        arma::vec p_vec { 1 / (1 + arma::exp(- mat2vec(x_ * beta) - offset_)) };
        // special care prevents coef_ diverging
        // reference: Friedman, J., Hastie, T., & Tibshirani, R. (2010)
        set_pmin_bound(p_vec, pmin);
        return p_vec;
    }
    inline arma::vec LogisticReg::predict(const arma::vec& beta,
                                          const double pmin = 1e-5) const
    {
        arma::vec beta0 { beta };
        if (standardize_) {
            beta0 = rev_rescale_coef(beta0);
        }
        return linkinv(beta0, pmin);
    }

    // define objective function (negative log-likehood function)
    inline double LogisticReg::objective(const arma::vec& beta) const
    {
        double res { 0 };
        arma::vec tmp { arma::zeros(2) };
        for (size_t i { 0 }; i < x_.n_rows; ++i) {
            double x_beta { arma::as_scalar(x_.row(i) * beta + offset_(i)) };
            tmp[1] = x_beta;
            res += log_sum_exp(tmp) - y_(i) * x_beta;
        }
        return res / dn_obs_;
    }
    inline double LogisticReg::objective() const
    {
        return objective(coef0_);
    }

    // define gradient function
    inline arma::vec LogisticReg::gradient(const arma::vec& beta,
                                           const double pmin) const
    {
        return x_.t() * (linkinv(beta, pmin) - y_) / dn_obs_;
    }
    // define gradient function at k-th dimension
    inline double LogisticReg::gradient(const arma::vec& beta,
                                        const unsigned int k,
                                        const double pmin) const
    {
        return arma::mean((linkinv(beta, pmin) - y_) % x_.col(k));
    }

    // Firth-type score function
    inline arma::vec LogisticReg::firth_score(const arma::vec& pi_vec) const
    {
        arma::mat b_mat {x_};
        for (size_t i {0}; i < x_.n_rows; ++i) {
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
        arma::vec res { arma::zeros(x_.n_cols) };
        for (size_t i {0}; i < x_.n_rows; ++i) {
            arma::rowvec tmp {
                (y_(i) - pi_vec(i) + hat_vec(i) * (0.5 - pi_vec(i))) * x_.row(i)
            };
            res += tmp.t();
        }
        return res / dn_obs_;
    }

    // define objective function and overwrites graidient
    inline double LogisticReg::objective(const arma::vec& beta,
                                         arma::vec& grad) const
    {
        arma::vec x_beta {x_ * beta + offset_};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        grad = x_.t() * (exp_x_beta / (1 + exp_x_beta) - y_);
        arma::vec y_x_beta {y_ % x_beta};
        double res {
            arma::as_scalar(arma::mean(arma::log(1 + exp_x_beta) - y_x_beta))
        };
        return res;
    }

    // compute iteration matrix in Bohning and Lindsay (1988)
    inline void LogisticReg::compute_bl_iter_mat(
        const bool force_update = false
        )
    {
        if (force_update || bl_iter_mat_.is_empty()) {
            bl_iter_mat_ = 4 * arma::inv_sympd(x_.t() * x_);
        }
    }

    // fitting regular logistic model by monotonic quadratic approximation
    // algorithm non-integer y_ vector is allowed
    // reference: Bohning and Lindsay (1988) SIAM
    inline void LogisticReg::fit(const arma::vec& start = 0,
                                 const unsigned int max_iter = 1000,
                                 const double rel_tol = 1e-6,
                                 const double pmin = 1e-5,
                                 const bool early_stop = false,
                                 const bool verbose = false)
    {
        arma::vec beta0 { arma::zeros(x_.n_cols) };
        if (start.n_elem == x_.n_cols) {
            beta0 = start;
        }
        double ell { arma::datum::inf };
        if (verbose) {
            Rcpp::Rcout << "\n" << std::string(40, '=')
                        << "\nStarting from\n"
                        << arma2rvec(beta0)
                        << "\n.";
            ell = objective(beta0);
            Rcpp::Rcout << "\nThe negative log-likelihood is "
                        << ell
                        << "\n.";
        }
        arma::vec beta { beta0 };
        arma::vec eta { arma::zeros(n_obs_) };
        arma::vec y_hat { eta };
        compute_bl_iter_mat();
        arma::mat iter_mat { bl_iter_mat_ * x_.t() };
        size_t i {0};

        while (i < max_iter) {
            y_hat = linkinv(beta0, pmin);
            beta = beta0 + iter_mat * (y_ - y_hat);
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(40, '=')
                            << "\nitartion: " << i + 1
                            << "\n  estimated coef_: "
                            << arma2rvec(beta)
                            << "\n.";
            }
            // if relative tolerance is statisfied
            if (rel_l1_norm(beta, beta0) < rel_tol) {
                double ell_old { ell };
                ell = objective(beta);
                if (verbose) {
                    Rcpp::Rcout << "\n  The negative log-likelihood changed\n";
                    Rprintf("  from %15.15f\n", ell_old);
                    Rprintf("    to %15.15f\n", ell);
                    Rcpp::Rcout << "\nReached convergence criterion.\n";
                }
                break;
            }
            // if early stop
            if (early_stop) {
                double ell_old { ell };
                ell = objective(beta);
                if (verbose) {
                    Rcpp::Rcout << "\n  The negative log-likelihood changed\n";
                    Rprintf("  from %15.15f\n", ell_old);
                    Rprintf("    to %15.15f\n", ell);
                }
                if (ell_old < ell) {
                    if (verbose) {
                        Rcpp::Rcout
                            << "\nEarly stopped the algorithm"
                            << " with estimates from"
                            << " iteration " << i << ".\n";
                    }
                    beta = beta0;
                    ell = ell_old;
                    break;
                }
            }
            // update beta
            beta0 = beta;
            i++;
        }
        coef0_ = beta;
        // rescale coef_ back
        rescale_coef();
        // compute score and prob
        xbeta_ = x_ * beta;
        prob_vec_ = linkinv(beta);
        // compute negative log-likelihood
        neg_ll_ = ell;
        coef_df_ = beta.n_elem;
    }

    // fitting Firth logistic model with monotonic quadratic approximation
    // algorithm;  non-integer y_ vector is allowed.
    inline void LogisticReg::firth_fit(const arma::vec& start = 0,
                                       const unsigned int max_iter = 300,
                                       const double rel_tol = 1e-5,
                                       const double pmin = 1e-5
        )
    {
        arma::vec beta0 { arma::zeros(x_.n_cols) };
        if (start.n_elem == x_.n_cols) {
            beta0 = start;
        }
        arma::vec beta { beta0 }, score_vec { beta0 };
        arma::vec eta { arma::zeros(y_.n_elem) };
        arma::vec y_hat { eta };
        compute_bl_iter_mat();
        size_t i {0};
        while (i < max_iter) {
            y_hat = linkinv(beta0, pmin);
            score_vec = firth_score(y_hat);
            beta = beta0 + bl_iter_mat_ * score_vec;
            if (rel_l1_norm(beta, beta0) < rel_tol) {
                break;
            }
            // update beta
            beta0 = beta;
            i++;
        }
        coef0_ = beta;
        // rescale coef_ back
        rescale_coef();
        // compute score and prob
        xbeta_ = x_ * beta;
        prob_vec_ = linkinv(beta);
        // compute negative log-likelihood
        neg_ll_ = objective();
        coef_df_ = beta.n_elem;
    }

    // compute CMD lowerbound vector
    inline void LogisticReg::set_cmd_lowerbound()
    {
        if (standardize_) {
            cmd_lowerbound_ = arma::ones<arma::rowvec>(p1_) * 0.25;
        } else {
            cmd_lowerbound_ = arma::sum(arma::square(x_), 0) / (4 * dn_obs_);
        }
    }

    // run one cycle of coordinate descent over a given active set
    inline void LogisticReg::regularized_fit_update(
        arma::vec& beta,
        arma::uvec& is_active,
        const double l1_lambda,
        const double l2_lambda,
        const arma::vec& l1_penalty_factor,
        const bool update_active = false,
        const double pmin = 1e-5,
        const bool early_stop = false,
        const bool verbose = false
        )
    {
        // compute lowerbound vector used in CMD algorithm
        double dlj { 0 };
        arma::vec beta_old = beta;
        for (size_t j {0}; j < p1_; ++j) {
            if (is_active[j] > 0) {
                dlj = gradient(beta, j, pmin);
                // update beta
                double numer {
                    soft_threshold(cmd_lowerbound_[j] * beta[j] - dlj,
                                   l1_penalty_factor[j] * l1_lambda)
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
        }
        // if early stop
        if (early_stop) {
            double ell_old { objective(beta_old) };
            ell_old = ell_old +
                l1_lambda * l1_norm(beta_old % l1_penalty_factor) +
                l2_lambda * sum_of_square(beta_old.tail(p0_));
            double ell_new { objective(beta) };
            ell_new = ell_new +
                l1_lambda * l1_norm(beta % l1_penalty_factor) +
                l2_lambda * sum_of_square(beta.tail(p0_));
            if (verbose) {
                Rcpp::Rcout << "The objective function changed\n";
                Rprintf("  from %15.15f\n", ell_old);
                Rprintf("    to %15.15f\n", ell_new);
            }
            if (ell_new > ell_old) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "the objective function somehow increased\n";
                    Rcpp::Rcout << "\nEarly stopped the CMD iterations "
                                << "with estimates from the last step.\n";
                }
                beta = beta_old;
            }
        }
    }

    // run a complete cycle of CMD for a given active set and lambda
    inline void LogisticReg::reg_active_fit(
        arma::vec& beta,
        const arma::uvec& is_active,
        const double l1_lambda,
        const double l2_lambda,
        const arma::vec& l1_penalty_factor,
        const bool varying_active_set = false,
        const unsigned int max_iter = 200,
        const double rel_tol = 1e-5,
        const double pmin = 1e-5,
        const bool early_stop = false,
        const bool verbose = false
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
                                           l2_lambda, l1_penalty_factor, true,
                                           pmin, early_stop, verbose);
                    if (rel_l1_norm(beta, beta0) < rel_tol) {
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                // run a full cycle over the converged beta
                regularized_fit_update(beta, is_active_new, l1_lambda,
                                       l2_lambda, l1_penalty_factor, true,
                                       pmin, early_stop, verbose);
                // check two active sets coincide
                if (l1_norm(is_active_new - is_active_stored) > 0) {
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
                                       l2_lambda, l1_penalty_factor, false,
                                       pmin, early_stop, verbose);
                if (rel_l1_norm(beta, beta0) < rel_tol) {
                    break;
                }
                beta0 = beta;
                i++;
            }
        }
    }

    // regularized logistic model by coordinate-majorization-descent algorithm
    // for particular lambda's for lasso penalty and ridge penalty
    // lambda_1 * factor * lasso + lambda_2 * ridge
    inline void LogisticReg::regularized_fit(
        const double l1_lambda = 0,
        const double l2_lambda = 0,
        const arma::vec& l1_penalty_factor = 0,
        const arma::vec& start = 0,
        const unsigned int max_iter = 300,
        const double rel_tol = 1e-5,
        const double pmin = 1e-5,
        const bool early_stop = false,
        const bool verbose = false
        )
    {
        set_cmd_lowerbound();
        // set penalty terms
        arma::vec l1_penalty { arma::ones(p0_) };
        if (l1_penalty_factor.n_elem == p0_) {
            l1_penalty = l1_penalty_factor * x_.n_cols /
                arma::sum(l1_penalty_factor);
        }
        l1_penalty_factor_ = l1_penalty;

        arma::vec beta { arma::zeros(x_.n_cols) };
        arma::vec grad_zero { arma::abs(gradient(beta, pmin)) };
        arma::vec grad_beta { grad_zero }, strong_rhs { grad_beta };

        // large enough lambda for all-zero coef_ (except intercept_)
        // excluding variable with zero penalty factor
        arma::uvec active_l1_penalty { arma::find(l1_penalty > 0) };
        grad_zero = grad_zero.tail(p0_);
        l1_lambda_max_ = arma::max(grad_zero.elem(active_l1_penalty) /
                                   l1_penalty.elem(active_l1_penalty));

        if (intercept_) {
            l1_penalty = arma::join_vert(arma::zeros(1), l1_penalty);
        }

        // get solution of l1_lambda_max_ for a warm start
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x_.n_cols) };
        if (intercept_) {
            // only needs to estimate intercept_
            is_active_strong(0) = 1;
            regularized_fit_update(beta, is_active_strong, l1_lambda_max_,
                                   l2_lambda, l1_penalty, false,
                                   pmin, early_stop, verbose);
        }
        coef0_ = beta;
        // rescale coef_ back
        rescale_coef();

        // early exit for lambda greater than lambda_max
        if (l1_lambda >= l1_lambda_max_) {
            en_coef_ = coef_;
            coef_df_ = 0;
            // compute score and prob
            xbeta_ = x_ * beta;
            prob_vec_ = linkinv(beta);
            // compute negative log-likelihood
            neg_ll_ = objective();
            // compute degree of freedom
            coef_df_ = get_coef_df(beta);
            // record other inputs
            l1_lambda_ = l1_lambda;
            l2_lambda_ = l2_lambda;
            return;
        }

        // use the input start if correctly specified
        if (start.n_elem == x_.n_cols) {
            beta = start;
        }

        // update active set by strong rule
        grad_beta = arma::abs(gradient(coef0_, pmin));
        strong_rhs = (2 * l1_lambda - l1_lambda_max_) * l1_penalty;

        for (size_t j { int_intercept_ }; j < p1_; ++j) {
            if (grad_beta(j) >= strong_rhs(j)) {
                is_active_strong(j) = 1;
            } else {
                beta(j) = 0;
            }
        }
        arma::uvec is_active_strong_new { is_active_strong };

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x_.n_cols > x_.n_rows || x_.n_cols > 50) {
            varying_active_set = true;
        }

        bool kkt_failed { true };
        strong_rhs = l1_lambda * l1_penalty;
        // eventually, strong rule will guess correctly
        while (kkt_failed) {
            // update beta
            reg_active_fit(beta, is_active_strong, l1_lambda, l2_lambda,
                           l1_penalty, varying_active_set, max_iter, rel_tol,
                           pmin, early_stop, verbose);
            // check kkt condition
            for (size_t j { int_intercept_ }; j < p1_; ++j) {
                if (is_active_strong(j)) {
                    continue;
                }
                if (std::abs(gradient(beta, j, pmin)) > strong_rhs(j)) {
                    // update active set
                    is_active_strong_new(j) = 1;
                }
            }
            if (l1_norm(is_active_strong - is_active_strong_new)) {
                is_active_strong = is_active_strong_new;
            } else {
                kkt_failed = false;
            }
        }
        // compute elastic net estimates, then rescale them back
        coef0_ = (1 + l2_lambda) * beta;
        rescale_coef();
        en_coef_ = coef_;
        // overwrite the naive elastic net estimate
        coef0_ = beta;
        rescale_coef();
        // compute score and prob
        xbeta_ = x_ * beta;
        prob_vec_ = linkinv(beta);
        // compute negative log-likelihood
        neg_ll_ = objective();
        // compute degree of freedom
        coef_df_ = get_coef_df(beta);
        // record other inputs
        l1_lambda_ = l1_lambda;
        l2_lambda_ = l2_lambda;
    }


    // for a sequence of lambda's
    // lambda * (penalty_factor * alpha_ * lasso + (1 - alpha_) / 2 * ridge)
    inline void LogisticReg::regularized_fit(
        arma::vec lambda = 0,
        const double alpha = 1,
        const unsigned int nlambda = 1,
        double lambda_min_ratio = 1e-4,
        const arma::vec& l1_penalty_factor = 0,
        const unsigned int max_iter = 300,
        const double rel_tol = 1e-5,
        const double pmin = 1e-5,
        const bool early_stop = false,
        const bool verbose = false
        )
    {
        set_cmd_lowerbound();
        // check alpha_
        if ((alpha < 0) || (alpha > 1)) {
            throw std::range_error("Alpha must be between 0 and 1.");
        }
        alpha_ = alpha;
        arma::vec l1_penalty { arma::ones(p0_) };
        if (l1_penalty_factor.n_elem == p0_) {
            l1_penalty = l1_penalty_factor * x_.n_cols /
                arma::sum(l1_penalty_factor);
        }
        l1_penalty_factor_ = l1_penalty;

        // construct lambda sequence
        arma::vec beta { arma::zeros(x_.n_cols) };
        arma::vec grad_beta { arma::abs(gradient(beta, pmin)) };
        arma::vec strong_rhs;
        l1_lambda_max_ =
            arma::max(grad_beta.tail(l1_penalty.n_elem) / l1_penalty) /
            std::max(alpha_, 1e-2);
        // take unique lambda and sort descendingly
        lambda = arma::reverse(arma::unique(lambda));
        // construct lambda sequence
        arma::vec lambda_seq;
        if (nlambda > 1) {
            double log_lambda_max { std::log(l1_lambda_max_) };
            if (x_.n_cols > n_obs_) {
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
        lambda_vec_ = lambda_seq;

        // update penalty for intercept_
        if (intercept_) {
            l1_penalty = arma::join_vert(arma::zeros(1), l1_penalty);
        }

        // initialize the estimate matrix
        coef_mat_ = arma::zeros(x_.n_cols, lambda_seq.n_elem);
        en_coef_mat_ = coef_mat_;
        neg_ll_vec_ = arma::zeros(lambda_seq.n_elem);
        coef_df_vec_ = arma::zeros<arma::uvec>(lambda_seq.n_elem);

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x_.n_cols) };

        // get solution of lambda_max for a warm start
        if (intercept_) {
            // only needs to estimate intercept_
            is_active_strong(0) = 1;
            regularized_fit_update(
                beta, is_active_strong, l1_lambda_max_ * alpha_,
                l1_lambda_max_ * (1 - alpha_) / 2, l1_penalty, false,
                pmin, early_stop, verbose
                );
        }
        coef0_ = beta;
        // rescale coef_ back
        rescale_coef();

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x_.n_cols > n_obs_) {
            varying_active_set = true;
        }

        // outer loop for the lambda sequence
        for (size_t k {0}; k < lambda_seq.n_elem; ++k) {
            // early exit for large lambda greater than lambda_max
            if (alpha_ * lambda_seq(k) >= l1_lambda_max_) {
                coef_mat_.col(k) = coef_;
                continue;
            }
            // update acitve set by strong rule (for lambda < lamda_max)
            grad_beta = arma::abs(gradient(beta, pmin));
            if (k == 0) {
                // use lambda_max
                strong_rhs = alpha_ *
                    (2 * lambda_seq(k) - l1_lambda_max_) * l1_penalty;
            } else {
                // use the last lambda
                strong_rhs = alpha_ *
                    (2 * lambda_seq(k) - lambda_seq(k - 1)) * l1_penalty;
            }
            for (size_t j { int_intercept_ }; j < p1_; ++j) {
                if (grad_beta(j) > strong_rhs(j)) {
                    is_active_strong(j) = 1;
                } else {
                    beta(j) = 0;
                }
            }
            arma::uvec is_active_strong_new { is_active_strong };
            strong_rhs = alpha_ * lambda_seq(k) * l1_penalty;
            bool kkt_failed { true };
            // eventually, strong rule will guess correctly
            while (kkt_failed) {
                // update beta
                reg_active_fit(beta, is_active_strong, lambda_seq(k) * alpha_,
                               lambda_seq(k) * (1 - alpha_) / 2, l1_penalty,
                               varying_active_set, max_iter, rel_tol,
                               pmin, early_stop, verbose);
                // check kkt condition
                for (size_t j { int_intercept_ }; j < p1_; ++j) {
                    if (is_active_strong(j) > 0) {
                        continue;
                    }
                    if (std::abs(gradient(beta, j, pmin)) > strong_rhs(j)) {
                        // update active set
                        is_active_strong_new(j) = 1;
                    }
                }
                if (l1_norm(is_active_strong - is_active_strong_new)) {
                    is_active_strong = is_active_strong_new;
                } else {
                    kkt_failed = false;
                }
            }
            // compute elastic net estimates
            coef0_ = (1 + (1 - alpha_) * lambda_seq(k) / 2) * beta;
            rescale_coef();
            en_coef_mat_.col(k) = coef_;
            // compute naive elastic net estimates
            coef0_ = beta;
            rescale_coef();
            coef_mat_.col(k) = coef_;
            // compute negative log-likelihood
            neg_ll_vec_(k) = objective();
            // compute degree of freedom
            coef_df_vec_(k) = get_coef_df(beta);
        }
        // prepare outputs

    }

}


#endif
