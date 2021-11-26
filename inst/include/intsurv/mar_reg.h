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

#ifndef INTSURV_MAR_REG_H
#define INTSURV_MAR_REG_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {

    class MarReg
    {
    public:
        arma::mat x_;           // standardized x with possible intercept
        bool standardize_;      // is x standardized
        arma::rowvec x_center_;
        arma::rowvec x_scale_;
        arma::vec a_;
        arma::vec b_;
        arma::vec a_bar_;
        arma::vec b_bar_;
        bool intercept_;

        // outputs
        unsigned int n_obs_;    // number of observation
        unsigned int p0_;       // number of predictors excluding the intercepts
        double l1_lambda_max_;
        arma::vec eta_;
        double alpha0_;
        arma::vec coef_;        // (eta_, alpha0)

        // cache variables
        arma::rowvec cmd_lowerbound_;
        unsigned int int_intercept_;
        unsigned int p1_;       // p0 + int_intercept_
        unsigned int p2_;       // p1_ + 1
        double dn_obs_;         // double version of n_obs_
        arma::vec xeta_;

        // control variables
        double l1_lambda_;
        double l2_lambda_;
        arma::vec l1_penalty_factor_;
        double pmin_ = 1e-5;

        // default
        MarReg() {}

        MarReg(const arma::mat& x,
               const bool intercept = true,
               const bool standardize = true) :
            x_ (x),
            standardize_ (standardize),
            intercept_ (intercept)
        {
            int_intercept_ = static_cast<unsigned int>(intercept_);
            n_obs_ = x_.n_rows;
            dn_obs_ = static_cast<double>(n_obs_);
            p0_ = x_.n_cols;
            p1_ = p0_ + int_intercept_;
            p2_ = p1_ + 1;
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
                        throw std::range_error(
                            "The design 'x' contains constant column.");
                    }
                }
            }
            if (intercept_) {
                x_ = arma::join_horiz(arma::ones(n_obs_), x_);
            }
            cmd_lowerbound_ = arma::mean(arma::square(x_), 0) / 2.0;
        }

        // setter
        inline MarReg* set_a(const arma::vec& a)
        {
            if (a.n_elem == n_obs_) {
                a_ = a;
            }
            return this;
        }
        inline MarReg* set_b(const arma::vec& b)
        {
            if (b.n_elem == n_obs_) {
                b_ = b;
            }
            return this;
        }
        inline MarReg* set_a_bar(const arma::vec& a_bar)
        {
            if (a_bar.n_elem == n_obs_) {
                a_bar_ = a_bar;
            }
            return this;
        }
        inline MarReg* set_b_bar(const arma::vec& b_bar)
        {
            if (b_bar.n_elem == n_obs_) {
                b_bar_ = b_bar;
            }
            return this;
        }

        // transfer coef for standardized data to coef for non-standardized data
        inline void rescale_eta()
        {
            if (standardize_) {
                arma::vec eta0 { eta_ };
                if (intercept_) {
                    eta_(0) = eta0(0) - arma::as_scalar((x_center_ / x_scale_) *
                                                        eta0.tail_rows(p0_));
                    for (size_t j {1}; j < p1_; ++j) {
                        eta_(j) = eta0(j) / x_scale_(j - 1);
                    }
                } else {
                    for (size_t j {0}; j < p0_; ++j) {
                        eta_(j) = eta0(j) / x_scale_(j);
                    }
                }
            }
        }

        // objective function
        inline double objective(const arma::vec& eta,
                                const double alpha0) const
        {
            arma::vec xeta { mat2vec(x_ * eta) };
            return - arma::mean(
                alpha0 * a_ + (a_ + a_bar_) % xeta -
                (a_ + b_) % arma::log(1 + arma::exp(xeta + alpha0)) -
                (a_bar_ + b_bar_) % arma::log(1 + arma::exp(xeta)));
        }
        inline double objective() const
        {
            return - arma::mean(
                alpha0_ * a_ + (a_ + a_bar_) % xeta_ -
                (a_ + b_) % arma::log(1 + arma::exp(xeta_ + alpha0_)) -
                (a_bar_ + b_bar_) % arma::log(1 + arma::exp(xeta_)));
        }

        // inverse link function
        inline arma::vec linkinv(const arma::vec& xeta) const
        {
            arma::vec p_vec {
                1.0 / (1.0 + arma::exp(- xeta))
            };
            // special care prevents coef diverging
            // reference: Friedman, J., Hastie, T., & Tibshirani, R. (2010)
            arma::vec::iterator it_end { p_vec.end() };
            for (arma::vec::iterator it { p_vec.begin() }; it != it_end; ++it) {
                if (*it < pmin_) {
                    *it = pmin_;
                } else if (*it > 1 - pmin_) {
                    *it = 1 - pmin_;
                }
            }
            return p_vec;
        }

        // return fitted q_{j,1} and q_{j,0}
        inline arma::vec q_vec(const double z = 0.0) const
        {
            return linkinv(xeta_ + z * alpha0_);
        }

        // compute gradient
        inline double cmd_grad_eta(const arma::vec& xeta,
                                   const double alpha0,
                                   const unsigned int l) const
        {
            return - arma::mean(
                x_.col(l) % (
                    (a_ + a_bar_) -
                    (a_ + b_)  % linkinv(xeta + alpha0) -
                    (a_bar_ + b_bar_) % linkinv(xeta)
                    )
                );
        }
        inline arma::vec cmd_grad_eta(const arma::vec& xeta,
                                      const double alpha0) const
        {
            arma::vec p1 { linkinv(xeta + alpha0) };
            arma::vec p0 { linkinv(xeta) };
            arma::vec tmp { ((a_ + a_bar_) -
                             (a_ + b_) % std::move(p1) -
                             (a_bar_ + b_bar_) % std::move(p0)) };
            return - x_.t() * tmp / dn_obs_;
        }

        inline double cmd_grad_alpha0(const arma::vec& xeta,
                                      const double alpha0) const
        {
            return - arma::mean(a_ - (a_ + b_) % linkinv(xeta + alpha0));
        }

        // one full cycle of coordinate descent
        inline void run_one_full_cycle(
            arma::vec& eta,
            double& alpha0,
            arma::vec& xeta,
            const double l1_lambda,
            const double l2_lambda,
            const arma::vec& l1_penalty_factor
            ) const
        {
            // update eta
            for (size_t l {0}; l < p1_; l++) {
                double dl { cmd_grad_eta(xeta, alpha0, l) };
                double eta_l_old { eta(l) };
                double numer {
                    soft_threshold(cmd_lowerbound_(l) * eta_l_old - dl,
                                   l1_penalty_factor(l) * l1_lambda)
                };
                if (numer == 0) {
                    eta(l) = 0.0;
                } else {
                    double denom { cmd_lowerbound_(l) + 2 * l2_lambda };
                    eta(l) = numer / denom;
                }
                xeta += (eta(l) - eta_l_old) * x_.col(l);
            }
            // update alpha0
            double dl { cmd_grad_alpha0(xeta, alpha0) };
            double numer { soft_threshold(0.5 * alpha0 - dl, l1_lambda) };
            if (numer == 0) {
                alpha0 = 0.0;
            } else {
                double denom { 0.5 + 2.0 * l2_lambda };
                alpha0 = numer / denom;
            }
        }

        inline void run_cmd_full_cycles(
            arma::vec& eta,
            double& alpha0,
            arma::vec& xeta,
            const double l1_lambda,
            const double l2_lambda,
            const arma::vec& l1_penalty_factor,
            const unsigned int max_iter,
            const double rel_tol
            ) const
        {
            arma::vec eta_old { eta };
            double alpha0_old { alpha0 };
            for (size_t i {0}; i < max_iter; ++i) {
                run_one_full_cycle(eta, alpha0, xeta, l1_lambda, l2_lambda,
                                   l1_penalty_factor);
                if (std::abs(alpha0 - alpha0_old) +
                    rel_l1_norm(eta, eta_old) < rel_tol) {
                    break;
                }
                eta_old = eta;
                alpha0_old = alpha0;
            }
        }

        // fit with optional regularization
        inline void regularized_fit(
            const arma::vec& a,
            const arma::vec& b,
            const arma::vec& a_bar,
            const arma::vec& b_bar,
            const arma::vec& start = arma::vec(),
            const double l1_lambda = 0,
            const double l2_lambda = 0,
            const arma::vec& l1_penalty_factor = arma::vec(),
            const unsigned int max_iter = 200,
            const double rel_tol = 1e-5,
            const double pmin = 1e-5
            )
        {
            a_ = a;
            b_ = b;
            a_bar_ = a_bar;
            b_bar_ = b_bar;
            l1_lambda_ = l1_lambda;
            l2_lambda_ = l2_lambda;
            pmin_ = pmin;
            // set initial values
            arma::vec eta, xeta, grad_eta, strong_rhs;
            double alpha0, grad_alpha0;
            // l1 penalty factor
            l1_penalty_factor_ = arma::ones(p2_);
            if (l1_penalty_factor.n_elem == p2_) {
                l1_penalty_factor_ = l1_penalty_factor * p2_ /
                    arma::sum(l1_penalty_factor);
            } else {
                // by default, do not penalize the intercept
                l1_penalty_factor_(0) = 0.0;
            }
            // starting values
            const bool has_start { start.n_elem == p2_ };
            if (has_start) {
                eta = start.head_rows(p1_);
                alpha0 = start(p1_);
                xeta = mat2vec(x_ * eta);
                grad_eta = cmd_grad_eta(xeta, alpha0);
                grad_alpha0 = cmd_grad_alpha0(xeta, alpha0);
            } else {
                eta = arma::zeros(p1_);
                alpha0 = 0.0;
                xeta = arma::zeros(n_obs_);
                grad_eta = arma::zeros(p1_);
                grad_alpha0 = 0.0;
            }
            l1_lambda_max_ = 0;
            if (l1_penalty_factor_(p1_) > 0) {
                // if alpha0 has positive l1 penalty
                l1_lambda_max_ = std::abs(grad_alpha0) /
                    l1_penalty_factor_(p1_);
            }
            for (size_t k {0}; k < p1_; ++k) {
                if (l1_penalty_factor_(k) > 0) {
                    l1_lambda_max_ = std::max(l1_lambda_max_,
                                              std::abs(grad_eta(k)));
                }
            }
            // TODO: add strong rule and active set updating
            run_cmd_full_cycles(eta, alpha0, xeta,
                                l1_lambda_, l2_lambda_, l1_penalty_factor_,
                                max_iter, rel_tol);
            eta_ = eta;
            alpha0_ = alpha0;
            rescale_eta();
            coef_ = eta_;
            coef_.resize(eta_.n_elem + 1);
            coef_(eta.n_elem) = alpha0;
            xeta_ = xeta;
            return;
        }

    };

}  // Intsurv


#endif /* INTSURV_MAR_REG_H */
