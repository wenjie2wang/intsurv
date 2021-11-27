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

#ifndef INTSURV_COXPH_REG_H
#define INTSURV_COXPH_REG_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {

    class CoxphReg {
    public:
        // model data ========================================================
        // all the following data are sorted accordingly
        arma::mat x_;           // (sorted and standardized) design matrix
        arma::vec time_;        // (sorted) observed times
        arma::vec event_;       // (sorted) event indicators

        // internals that depend on the data =================================
        unsigned int n_obs_;    // number of observations
        double dn_obs_;         // double version of n_obs
        unsigned int p_;        // number of given predictors
        arma::uvec ord_;        // index sorting the rows by time and event
        arma::uvec rev_ord_;    // index reverting the sorting
        bool has_ties_ {false}; // if there exists ties on event_ times
        arma::rowvec x_center_; // the column center of x
        arma::rowvec x_scale_;  // the scale of x
        // at each unique time point
        arma::uvec event_ind_;  // indices of event times
        // index indicating the first record on each distinct event time
        arma::uvec uni_event_ind_;
        // index indicating the first record on each time if event/censored
        arma::uvec uni_time_ind_;
        // size of risk-set at each time point
        arma::vec riskset_size_time_;
        // at each unique event time point
        arma::vec unique_time_;
        arma::vec d_time0_;      // event times
        arma::vec d_time_;       // distinct event times
        arma::mat d_x_;          // design matrix aggregated at d_time
        arma::vec d_offset_;     // offset aggregated at d_time
        arma::vec d_offset_haz_; // offset_haz aggregated at d_time
        arma::vec delta_n_;      // event counts at d_time
        arma::vec riskset_size_; // size of risk-set at d_time

        // inverse matrix of formula (5.9) in Bohning and Lindsay (1988) SIAM
        arma::mat inv_bl_cov_lowerbound_;
        // for regularized regression
        arma::vec cmd_lowerbound_;

        // control parameters that do not depend on data =====================
        bool standardize_;      // is x standardized
        arma::vec l1_penalty_factor_; // adaptive weights for lasso penalty
        arma::vec offset_;      // offset term
        arma::vec offset_haz_;  // offset term for baseline hazard only
        double l1_lambda_;      // tuning parameter for lasso penalty
        double l2_lambda_;      // tuning parameter for ridge penalty
        double alpha_;          // tuning parameter
        arma::vec lambda_vec_;   // lambda sequence

        // outputs ===========================================================
        arma::vec coef0_;       // coef for the standardized x
        double l1_lambda_max_;  // the "big enough" lambda => zero coef

        // for a single l1_lambda_ and l2_lambda_
        arma::vec coef_;        // covariate coefficient estimates
        arma::vec en_coef_;     // (rescaled) elastic net estimates
        double neg_ll_;         // partial negative log-likelihood
        unsigned int coef_df_;  // number of non-zero coef estimates
        arma::vec xbeta_;       // sorted x * coef_
        double bic_;            // log(num_event) * coef_df + 2 * neg_ll

        // for a lambda sequence
        arma::mat coef_mat_;     // coef_ matrix (rescaled for origin x_)
        arma::mat en_coef_mat_;  // elastic net estimates
        arma::vec neg_ll_vec_;   // negative log-likelihood vector
        arma::uvec coef_df_vec_; // coef df vector
        arma::vec bic_vec_;      // log(num_event) * coef_df_ + 2 * neg_ll_

        // hazard and survival estimates at every time point (unique or not)
        arma::vec h0_time_;
        arma::vec S0_time_;
        arma::vec H0_time_;
        arma::vec h_time_;
        arma::vec H_time_;
        arma::vec S_time_;
        arma::vec hc_time_;
        arma::vec Hc_time_;
        arma::vec Sc_time_;
        // baseline estimates at unique times only
        arma::vec h0_est_;
        arma::vec H0_est_;
        arma::vec S0_est_;
        arma::vec hc_est_;
        arma::vec Hc_est_;
        arma::vec Sc_est_;

        // default constructor
        CoxphReg() {}

        // constructors
        CoxphReg(const arma::vec& time,
                 const arma::vec& event,
                 const arma::mat& x,
                 const bool standardize = true)
        {
            // sort based on time and event
            // time: ascending order
            // event: events first, then censoring at the same time point
            arma::uvec des_event_ind { arma::sort_index(event, "descend") };
            arma::uvec asc_time_ind {
                arma::stable_sort_index(time.elem(des_event_ind), "ascend")
            };
            ord_ = des_event_ind.elem(asc_time_ind);
            rev_ord_ = arma::sort_index(ord_);
            // do actual sorting
            time_ = time.elem(ord_);
            event_ = event.elem(ord_);
            x_ = x.rows(ord_);
            n_obs_ = x_.n_rows;
            dn_obs_ = static_cast<double>(n_obs_);
            p_ = x_.n_cols;
            // standardize covariates
            standardize_ = standardize;
            if (standardize_) {
                x_center_ = arma::mean(x_);
                x_scale_ = arma::stddev(x_, 1);
                for (size_t j {0}; j < x_.n_cols; ++j) {
                    if (x_scale_(j) > 0) {
                        x_.col(j) = (x_.col(j) - x_center_(j)) / x_scale_(j);
                    } else {
                        // coef will be zero, set non-zero for rescaling
                        x_.col(j) = arma::zeros(n_obs_);
                        x_scale_(j) = - 1.0;
                    }
                }
            }
            // binary event indicator
            event_ind_ = arma::find(event_ > 0);
            // check if there exists ties on event times
            d_time0_ = time_.elem(event_ind_);
            has_ties_ = any_duplicated(d_time0_);
            // default value when no tied event times
            uni_event_ind_ = event_ind_;
            d_time_ = d_time0_;
            delta_n_ = event_.elem(event_ind_);
            d_x_ = x_.rows(event_ind_);
            for (size_t j {0}; j < x_.n_cols; ++j) {
                d_x_.col(j) %= delta_n_;
            }
            uni_time_ind_ = find_first_unique(time_);
            if (has_ties_) {
                // update uni_event_ind_
                uni_event_ind_ = vec_intersection(uni_time_ind_, event_ind_);
                d_time_ = time_.elem(uni_event_ind_);
                // aggregate at distinct event times
                delta_n_ = aggregate_sum(delta_n_, d_time0_);
                d_x_ = aggregate_sum(d_x_, d_time0_);
            }
            // initialize offset
            reset_offset();
            reset_offset_haz();
            riskset_size_ = arma::ones(time_.n_elem);
            riskset_size_ = cum_sum(riskset_size_, true).elem(uni_event_ind_);
        }

        // re-scale coef_ back
        inline void rescale_coef()
        {
            coef_ = coef0_;
            if (standardize_) {
                for (size_t j {0}; j < coef0_.n_elem; ++j) {
                    coef_[j] = coef0_[j] / x_scale_[j];
                }
            }
        }

        // partially update event non-zero-one indicators
        // assuming only old event weights strictly > 0 need update
        inline void update_event_weight(const arma::vec& event,
                                        const bool is_sorted = true)
        {
            event_ = event;
            if (! is_sorted) {
                event_ = event_.elem(ord_);
            }
            // need to update delta_n_ and d_x_
            delta_n_ = event_.elem(event_ind_);
            d_x_ = x_.rows(event_ind_);
            for (size_t j {0}; j < x_.n_cols; ++j) {
                d_x_.col(j) %= delta_n_;
            }
            if (has_ties_) {
                // aggregate at distinct event_ times
                delta_n_ = aggregate_sum(delta_n_, d_time0_);
                d_x_ = aggregate_sum(d_x_, d_time0_);
            }
        }

        // set offset
        inline void set_offset(const arma::vec& offset,
                               const bool is_sorted = true)
        {
            if (offset.n_elem == n_obs_) {
                offset_ = offset;
            } else if (offset.n_elem == 1 || offset.empty()) {
                reset_offset();
                return;
            } else {
                throw std::length_error(
                    "The length of the specified offset must match sample size."
                    );
            }
            if (! is_sorted) {
                // update offset_ for appropriate input
                offset_ = offset_.elem(ord_);
            }
            // update d_offset_ as well
            d_offset_ = offset_.elem(event_ind_) % event_.elem(event_ind_);
            if (has_ties_) {
                d_offset_ = aggregate_sum(d_offset_, d_time0_);
            }
        }
        inline void reset_offset()
        {
            offset_ = arma::zeros(n_obs_);
            if (has_ties_) {
                d_offset_ = arma::zeros(d_time_.n_elem);
            } else {
                d_offset_ = arma::zeros(event_ind_.n_elem);
            }
        }
        // set offset for denominator in baseline hazard function
        // for cure rate model
        inline void set_offset_haz(const arma::vec& offset,
                                   const bool is_sorted = true)
        {
            if (offset.n_elem == n_obs_) {
                offset_haz_ = offset;
            } else if (offset.n_elem == 1 || offset.empty()) {
                reset_offset_haz();
                return;
            } else {
                throw std::length_error(
                    "The length of offset must match sample size.");
            }
            if (! is_sorted) {
                // update offset_ for appropriate input
                offset_haz_ = offset_haz_.elem(ord_);
            }
            // update d_offset_haz as well
            d_offset_haz_ = offset_haz_.elem(event_ind_) %
                event_.elem(event_ind_);
            if (has_ties_) {
                d_offset_haz_ = aggregate_sum(d_offset_haz_, d_time0_);
            }
        }
        inline void reset_offset_haz()
        {
            offset_haz_ = arma::zeros(n_obs_);
            if (has_ties_) {
                d_offset_haz_ = arma::zeros(d_time_.n_elem);
            } else {
                d_offset_haz_ = arma::zeros(event_ind_.n_elem);
            }
        }

        // function that computes baseline estimates
        inline void compute_haz_surv_time(const arma::vec& beta);
        inline void compute_haz_surv_time();
        inline void compute_censor_haz_surv_time();

        // prepare hazard and survival estimates on unique time_ points
        inline void est_haz_surv();

        // function that computes objective function only
        inline double objective() const;
        inline double objective(const arma::vec& beta) const;

        // function that computes gradients only
        inline arma::vec gradient(const arma::vec& beta) const;
        inline double gradient(const arma::vec& beta,
                               const unsigned int k) const;

        // function that computes objective and overwrite gradients
        inline double objective(const arma::vec& beta,
                                arma::vec& grad) const;

        // function computing B&L covariance lower bound matrix
        inline void compute_inv_bl_cov_lowerbound(const bool force_update);

        // function computing B&L lower bound for step size
        inline double bl_step_lowerbound(const arma::mat& x,
                                         const arma::vec& h) const;

        // get the lower bound of second derivative in CMD algorithm
        inline void compute_cmd_lowerbound(const bool force_update);

        // some simple functions
        inline unsigned int sample_size() const { return n_obs_; }
        inline void compute_bic() {
            bic_ = std::log(arma::sum(event_)) * coef_df_ + 2 * neg_ll_;
        }

        // additional methods for coxph_cure
        // revserse the rescale process to get coef0_ from a new coef_
        inline void rev_rescale_coef()
        {
            coef0_ = coef_;
            if (standardize_) {
                for (size_t j {0}; j < coef_.n_elem; ++j) {
                    coef0_[j] = coef_[j] * x_scale_[j];
                }
            }
        }
        // update coef0_, en_coef_, and coef_df_ from a new coef_
        inline void set_en_coef(const double l2_lambda = 0) {
            // update coef0_
            rev_rescale_coef();
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
            coef_df_ = compute_coef_df(beta);
        }

        // fit regular Cox model
        inline void fit(const arma::vec& start,
                        const unsigned int max_iter,
                        const double rel_tol,
                        const bool early_stop,
                        const bool verbose);

        // one-full cycle of coordinate-majorization-descent
        inline void regularized_fit_update(arma::vec& beta,
                                           arma::uvec& is_active,
                                           const double l1_lambda,
                                           const double l2_lambda,
                                           const arma::vec& l1_penalty_factor,
                                           const bool update_active,
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
                                   const bool early_stop,
                                   const bool verbose);

        // fit regularized Cox model with adaptive lasso penalty
        // for a perticular lambda
        inline void regularized_fit(const double l1_lambda,
                                    const double l2_lambda,
                                    const arma::vec& l1_penalty_factor,
                                    const arma::vec& start,
                                    const unsigned int max_iter,
                                    const double rel_tol,
                                    const bool early_stop,
                                    const bool verbose);

        // overload for a sequence of lambda's
        inline void regularized_fit(arma::vec lambda,
                                    const double alpha,
                                    const unsigned int nlambda,
                                    double lambda_min_ratio,
                                    const arma::vec& l1_penalty_factor,
                                    const unsigned int max_iter,
                                    const double rel_tol,
                                    const bool early_stop,
                                    const bool verbose);

    };
    // end of class definition

    // compute baseline hazard function and its friends
    // here beta is coef_ estimate for non-standardized data
    inline void CoxphReg::compute_haz_surv_time(const arma::vec& beta)
    {
        if (standardize_) {
            if (! beta.is_empty() && beta.n_elem == x_.n_cols) {
                // re-scale the input beta
                arma::vec beta0 { beta % x_scale_.t() };
                xbeta_ = x_ * beta0 + arma::as_scalar(x_center_ * beta);
            } else {
                // use estimated coef_
                xbeta_ = x_ * coef0_ + arma::as_scalar(x_center_ * coef_);
            }
        } else {
            if (! beta.is_empty() && beta.n_elem == x_.n_cols) {
                xbeta_ = x_ * beta;
            } else {
                // use estimated coef_
                xbeta_ = x_ * coef_;
            }
        }
        arma::vec exp_risk_score { arma::exp(xbeta_ + offset_) };
        // 1. hazard rate function
        arma::vec h0_numer { aggregate_sum(event_, time_, false) };
        arma::vec h0_denom { exp_risk_score % arma::exp(offset_haz_) };
        h0_denom = aggregate_sum(h0_denom, time_, false, true, true);
        h0_time_ = h0_numer / h0_denom;
        h_time_ = h0_time_ % exp_risk_score;

        // 2. baseline cumulative hazard function
        H0_time_ = arma::zeros(h0_time_.n_elem);
        for (size_t i: uni_event_ind_) {
            H0_time_(i) = h0_time_(i);
        }
        H0_time_ = cum_sum(H0_time_);
        H_time_ = H0_time_ % exp_risk_score;

        // 3. baseline survival function
        S0_time_ = arma::exp(- H0_time_);
        S_time_ = arma::exp(- H_time_);
    }
    inline void CoxphReg::compute_haz_surv_time()
    {
        compute_haz_surv_time(coef_);
    }
    inline void CoxphReg::compute_censor_haz_surv_time()
    {
        arma::vec censor_ind { 1 - event_ };
        arma::vec delta_c { aggregate_sum(censor_ind, time_, false) };
        if (riskset_size_time_.is_empty()) {
            arma::vec tmp { arma::ones(time_.n_elem) };
            riskset_size_time_ =
                aggregate_sum(tmp, time_, false, true, true);
        }
        hc_time_ = delta_c / riskset_size_time_;
        Hc_time_ = arma::zeros(hc_time_.n_elem);
        for (size_t i: uni_time_ind_) {
            Hc_time_(i) = hc_time_(i);
        }
        Hc_time_ = cum_sum(Hc_time_);
        Sc_time_ = arma::exp(- Hc_time_);
    }

    // prepare hazard and survival estimates on unique time_ points
    // should be used after compute_haz_surv_time() and
    // compute_censor_haz_surv_time()
    inline void CoxphReg::est_haz_surv()
    {
        unique_time_ = time_.elem(uni_time_ind_);
        if (h0_time_.is_empty()) {
            compute_haz_surv_time();
        }
        if (hc_time_.is_empty()) {
            compute_censor_haz_surv_time();
        }
        h0_est_ = h0_time_.elem(uni_time_ind_);
        H0_est_ = H0_time_.elem(uni_time_ind_);
        S0_est_ = S0_time_.elem(uni_time_ind_);
        hc_est_ = hc_time_.elem(uni_time_ind_);
        Hc_est_ = Hc_time_.elem(uni_time_ind_);
        Sc_est_ = Sc_time_.elem(uni_time_ind_);
    }

    // the negative log-likelihood function based on the broslow's formula
    inline double CoxphReg::objective(const arma::vec& beta) const
    {
        const arma::vec dx_beta {
            d_x_ * beta + d_offset_ + d_offset_haz_
        };
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + offset_ + offset_haz_)
        };
        const arma::vec h0_denom {
            cum_sum(exp_x_beta, true)
        };
        arma::vec log_h0_denom_event {
            arma::log(h0_denom.elem(uni_event_ind_))
        };
        return - arma::sum(dx_beta - delta_n_ % log_h0_denom_event);
    }
    inline double CoxphReg::objective() const
    {
        return objective(coef0_);
    }

    // the gradient of negative loglikelihood function
    inline arma::vec CoxphReg::gradient(const arma::vec& beta) const
    {
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + offset_ + offset_haz_)
        };
        const arma::vec h0_denom { cum_sum(exp_x_beta, true) };
        arma::mat numer_mat { arma::zeros(x_.n_rows, x_.n_cols) };
        for (size_t i {0}; i < x_.n_rows; ++i) {
            numer_mat.row(i) = exp_x_beta(i) * x_.row(i);
        }
        numer_mat = cum_sum(numer_mat, true);
        arma::vec h0_denom_event {h0_denom.elem(uni_event_ind_)};
        numer_mat = numer_mat.rows(uni_event_ind_);
        for (size_t j {0}; j < x_.n_cols; ++j) {
            numer_mat.col(j) = numer_mat.col(j) % delta_n_ / h0_denom_event;
        }
        return - (arma::sum(d_x_, 0) - arma::sum(numer_mat, 0)).t();
    }
    // the gradient of negative loglikelihood function at k-th direction
    inline double CoxphReg::gradient(const arma::vec& beta,
                                     const unsigned int k) const
    {
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + offset_ + offset_haz_)
        };
        arma::vec h0_denom { cum_sum(exp_x_beta, true) };
        arma::vec numer { cum_sum(mat2vec(x_.col(k) % exp_x_beta), true) };
        h0_denom = h0_denom.elem(uni_event_ind_);
        numer = numer.elem(uni_event_ind_);
        return - arma::sum(d_x_.col(k) - delta_n_ % numer / h0_denom);
    }

    // the negative log-likelihood function based on the broslow's formula
    inline double CoxphReg::objective(const arma::vec& beta,
                                      arma::vec& grad) const
    {
        const arma::vec dx_beta {
            d_x_ * beta + d_offset_ + d_offset_haz_
        };
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + offset_ + offset_haz_)
        };
        const arma::vec h0_denom { cum_sum(exp_x_beta, true) };
        arma::mat numer_mat { arma::zeros(x_.n_rows, x_.n_cols) };
        for (size_t i {0}; i < x_.n_rows; ++i) {
            numer_mat.row(i) = exp_x_beta(i) * x_.row(i);
        }
        numer_mat = cum_sum(numer_mat, true);
        arma::vec h0_denom_event {h0_denom.elem(uni_event_ind_)};
        numer_mat = numer_mat.rows(uni_event_ind_);
        for (size_t j {0}; j < x_.n_cols; ++j) {
            numer_mat.col(j) = numer_mat.col(j) % delta_n_ / h0_denom_event;
        }
        // overwrite grad
        grad = - (arma::sum(d_x_, 0) - arma::sum(numer_mat, 0)).t();
        return - arma::sum(dx_beta - delta_n_ % arma::log(h0_denom_event));
    }

    // one lower bound of covariance matrix
    // reference: formula (5.9) in Bohning and Lindsay (1988) SIAM
    inline void CoxphReg::compute_inv_bl_cov_lowerbound(
        const bool force_update = false
        )
    {
        if (inv_bl_cov_lowerbound_.is_empty() || force_update) {
            arma::mat res { arma::zeros(x_.n_cols, x_.n_cols) };
            // compute x_ * x_^t and store them in a cube
            arma::cube sum_risk_xxt {
                arma::cube(x_.n_cols, x_.n_cols, time_.n_elem,
                           arma::fill::zeros)
            };
            for (size_t i {0}; i < time_.n_elem; ++i) {
                sum_risk_xxt.slice(i) = crossprod(x_.row(i));
            }
            sum_risk_xxt = cum_sum(sum_risk_xxt, true);
            sum_risk_xxt = sum_risk_xxt.slices(uni_event_ind_);
            arma::mat sum_risk_x { cum_sum(x_, true) };
            sum_risk_x = sum_risk_x.rows(uni_event_ind_);
            for (unsigned int i {0}; i < d_time_.n_elem; ++i) {
                res += (sum_risk_xxt.slice(i) - crossprod(sum_risk_x.row(i)) /
                        riskset_size_(i)) * (delta_n_(i) / 2);
            }
            inv_bl_cov_lowerbound_ = arma::inv_sympd(res);
        }
    }
    inline double CoxphReg::bl_step_lowerbound(const arma::mat& x,
                                               const arma::vec& h) const
    {
        // compute x * h and store them in a vector
        arma::vec htx { x * h };
        arma::vec max_risk_htx { cum_max(htx, true) };
        max_risk_htx = max_risk_htx.elem(uni_event_ind_);
        arma::vec min_risk_htx { cum_min(htx, true) };
        min_risk_htx = min_risk_htx.elem(uni_event_ind_);
        double Mi {0}, mi {0}, res {0};
        for (unsigned int i {0}; i < d_time_.n_elem; ++i) {
            Mi = max_risk_htx(i);
            mi = min_risk_htx(i);
            res += ((std::pow(Mi, 2) + std::pow(mi, 2)) / 4 - Mi * mi / 2) *
                delta_n_(i);
        }
        return res;
    }

    // compute lower bound of second derivative
    // in coordinate-majorization-descent algorithm (cmd)
    inline void CoxphReg::compute_cmd_lowerbound(
        const bool force_update = false
        )
    {
        if (cmd_lowerbound_.is_empty() || force_update) {
            // compute D_k at each event_ time_ point
            arma::mat max_risk_x { cum_max(x_, true) };
            arma::mat min_risk_x { cum_min(x_, true) };
            max_risk_x = max_risk_x.rows(uni_event_ind_);
            min_risk_x = min_risk_x.rows(uni_event_ind_);
            arma::vec res { arma::zeros(x_.n_cols) };
            for (size_t k {0}; k < x_.n_cols; ++k) {
                res(k) = arma::sum(
                    pow(max_risk_x.col(k) - min_risk_x.col(k), 2) % delta_n_
                    ) / (4 * dn_obs_);
            }
            cmd_lowerbound_ = res;
        }
    }

    // fit regular Cox model by monotonic quadratic approximation algorithm
    // that allows non-integer "event_" and tied events
    // reference: Bohning and Lindsay (1988) SIAM
    inline void CoxphReg::fit(const arma::vec& start = 0,
                              const unsigned int max_iter = 100,
                              const double rel_tol = 1e-6,
                              const bool early_stop = true,
                              const bool verbose = false)
    {
        compute_inv_bl_cov_lowerbound();
        arma::vec beta0 { arma::zeros(x_.n_cols) };
        if (start.n_elem == x_.n_cols) {
            beta0 = start;
        }
        arma::vec beta { beta0 }, h_vec { beta0 }, grad_vec { beta0 };
        size_t i {0};
        double b_new {0}, alpha_ {0}, ell {0}, ell_old {arma::datum::inf};
        while (i <= max_iter) {
            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // compute negative log-likelihood function and update gradient
            ell = objective(beta, grad_vec);
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(40, '=')
                            << "\niteration: " << i
                            << "\n  estimated coef: "
                            << arma2rvec(beta)
                            << "\n  negative log-likelihood: "
                            << ell
                            << "\n";
            }
            h_vec = - inv_bl_cov_lowerbound_ * grad_vec;
            b_new = bl_step_lowerbound(x_, h_vec);
            alpha_ = - arma::as_scalar(crossprod(h_vec, grad_vec)) / b_new;
            beta = beta0 + alpha_ * h_vec;

            // if tolerance is reached
            if (rel_l1_norm(beta, beta0) < rel_tol) {
                break;
            }
            // if the objective function somehow increased,
            // which is theorically impossible
            // but possible in practice due to numerical issue
            if (ell_old < ell) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "The negative log-likelihood increased\n";
                    Rprintf("  from %15.15f\n", ell_old);
                    Rprintf("    to %15.15f\n", ell);
                }
                if (early_stop) {
                    if (verbose) {
                        Rcpp::Rcout << "\nEarly stopped the algorithm"
                                    << " with estimates from"
                                    << " iteration " << i - 1 << "\n";
                    }
                    beta = beta0;
                    ell = ell_old;
                    break;
                }
            }
            // update-steps
            beta0 = beta;
            ell_old = ell;
            i++;
        }
        coef0_ = beta;
        rescale_coef();
        neg_ll_ = ell;
        coef_df_ = beta.n_elem;
        compute_bic();
    }

    // run one cycle of coordinate descent over a given active set
    inline void CoxphReg::regularized_fit_update(
        arma::vec& beta,
        arma::uvec& is_active,
        const double l1_lambda,
        const double l2_lambda,
        const arma::vec& l1_penalty_factor,
        const bool update_active = false,
        const bool early_stop = false,
        const bool verbose = false
        )
    {
        double dlj { 0 };
        arma::vec beta_old = beta;
        for (size_t j {0}; j < beta.n_elem; ++j) {
            if (is_active[j]) {
                dlj = gradient(beta, j) / dn_obs_;
                // update beta
                beta[j] = soft_threshold(
                    cmd_lowerbound_[j] * beta[j] - dlj,
                    l1_penalty_factor[j] * l1_lambda
                    ) /
                    (cmd_lowerbound_[j] + 2 * l2_lambda);
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
            ell_old = ell_old / dn_obs_ +
                l1_lambda * l1_norm(beta_old % l1_penalty_factor) +
                l2_lambda * sum_of_square(beta_old);
            double ell_new { objective(beta) };
            ell_new = ell_new / dn_obs_ +
                l1_lambda * l1_norm(beta % l1_penalty_factor) +
                l2_lambda * sum_of_square(beta);
            if (verbose) {
                Rcpp::Rcout << "The objective function changed\n";
                Rprintf("  from %15.15f\n", ell_old);
                Rprintf("    to %15.15f\n", ell_new);
            }
            if (ell_new > ell_old) {
                if (verbose) {
                    Rcpp::Rcout << "Warning: "
                                << "the objective function somehow increased.\n"
                                << "Early stopped the CMD iterations "
                                << "with estimates from the last step.\n";
                }
                beta = beta_old;
            }
        }
    }

    // run a complete cycle of CMD for a given active set and lambda
    inline void CoxphReg::reg_active_fit(
        arma::vec& beta,
        const arma::uvec& is_active,
        const double l1_lambda,
        const double l2_lambda,
        const arma::vec& l1_penalty_factor,
        const bool varying_active_set = false,
        const unsigned int max_iter = 300,
        const double rel_tol = 1e-5,
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
                                           early_stop, verbose);
                    if (rel_l1_norm(beta, beta0) < rel_tol) {
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                // run a full cycle over the converged beta
                regularized_fit_update(beta, is_active_new, l1_lambda,
                                       l2_lambda, l1_penalty_factor, true,
                                       early_stop, verbose);
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
                                       l2_lambda, l1_penalty_factor, false,
                                       early_stop, verbose);
                if (rel_l1_norm(beta, beta0) < rel_tol) {
                    break;
                }
                beta0 = beta;
                i++;
            }
        }
    }

    // fitting regularized Cox model with coordinate-majorizatio-descent
    // algorithm that allows non-integer "event_" and tied events
    // for a perticular l1_lambda_ and l2_lambda_
    // lambda_1 * lasso + lambda_2 * ridge
    inline void CoxphReg::regularized_fit(
        const double l1_lambda = 0,
        const double l2_lambda = 0,
        const arma::vec& l1_penalty_factor = 0,
        const arma::vec& start = 0,
        const unsigned int max_iter = 300,
        const double rel_tol = 1e-5,
        const bool early_stop = false,
        const bool verbose = false
        )
    {
        // compute lowerbound vector used in CMD algorithm
        compute_cmd_lowerbound();
        // set penalty terms
        arma::vec l1_penalty { arma::ones(x_.n_cols) };
        if (l1_penalty_factor.n_elem == x_.n_cols) {
            // re-scale so that sum(factor) = number of predictors
            l1_penalty = l1_penalty_factor * x_.n_cols /
                arma::sum(l1_penalty_factor);
        }
        l1_penalty_factor_ = l1_penalty;
        l1_lambda_ = l1_lambda;
        l2_lambda_ = l2_lambda;

        // the maximum (large enough) lambda that results in all-zero estimates
        arma::vec beta { arma::zeros(x_.n_cols) };
        arma::vec grad_beta { beta }, strong_rhs { beta }, grad_zero { beta };
        // excluding variable with zero penalty factor
        arma::uvec active_l1_penalty { arma::find(l1_penalty > 0) };
        grad_zero = arma::abs(gradient(beta));
        l1_lambda_max_ = arma::max(
            grad_zero.elem(active_l1_penalty) /
            l1_penalty.elem(active_l1_penalty)
            ) / x_.n_rows;

        // early exit for large lambda greater than lambda_max
        coef0_ = beta;
        coef_ = beta;
        if (l1_lambda > l1_lambda_max_) {
            // no need to rescale all-zero coef_
            en_coef_ = coef_;
            coef_df_ = 0;
            // compute negative log-likelihood
            neg_ll_ = objective();
            coef_df_ = compute_coef_df(beta);
            compute_bic();
            return;
        }

        // use the input starting value
        if (start.n_elem == x_.n_cols) {
            beta = start;
        }

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x_.n_cols) };
        // update active set by strong rule
        grad_beta = arma::abs(gradient(coef0_)) / x_.n_rows;
        strong_rhs = (2 * l1_lambda_ - l1_lambda_max_) * l1_penalty;
        for (size_t j {0}; j < x_.n_cols; ++j) {
            if (grad_beta(j) > strong_rhs(j)) {
                is_active_strong(j) = 1;
            } else {
                beta(j) = 0;
            }
        }
        arma::uvec is_active_strong_new { is_active_strong };

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x_.n_cols > x_.n_rows) {
            varying_active_set = true;
        }

        bool kkt_failed { true };
        strong_rhs = l1_lambda_ * l1_penalty;
        // eventually, strong rule will guess correctly
        while (kkt_failed) {
            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();

            // update beta
            reg_active_fit(beta, is_active_strong, l1_lambda_, l2_lambda_,
                           l1_penalty, varying_active_set, max_iter, rel_tol,
                           early_stop, verbose);
            // check kkt condition
            for (size_t j {0}; j < x_.n_cols; ++j) {
                if (is_active_strong(j)) {
                    continue;
                }
                if (std::abs(gradient(beta, j)) / x_.n_rows >
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
        coef0_ = (1 + l2_lambda_) * beta;
        rescale_coef();
        en_coef_ = coef_;
        // overwrite the naive elastic net estimate
        coef0_ = beta;
        rescale_coef();
        // compute negative log-likelihood
        neg_ll_ = objective();
        coef_df_ = compute_coef_df(beta);
        compute_bic();
    }


    // for a sequence of lambda's
    // lambda * (penalty_factor * alpha_ * lasso + (1 - alpha_) / 2 * ridge)
    inline void CoxphReg::regularized_fit(
        arma::vec lambda = 0,
        const double alpha = 1,
        const unsigned int nlambda = 1,
        double lambda_min_ratio = 1e-4,
        const arma::vec& l1_penalty_factor = 0,
        const unsigned int max_iter = 300,
        const double rel_tol = 1e-5,
        const bool early_stop = false,
        const bool verbose = false
        )
    {
        // check alpha_
        if ((alpha < 0) || (alpha > 1)) {
            throw std::range_error("Alpha must be between 0 and 1.");
        }
        alpha_ = alpha;

        // set penalty terms
        arma::vec l1_penalty { arma::ones(x_.n_cols) };
        if (l1_penalty_factor.n_elem == x_.n_cols) {
            // re-scale so that sum(factor) = number of predictors
            l1_penalty = l1_penalty_factor * x_.n_cols /
                arma::sum(l1_penalty_factor);
        }
        l1_penalty_factor_ = l1_penalty;

        // the maximum (large enough) lambda that results in all-zero estimates
        arma::vec beta { arma::zeros(x_.n_cols) };
        arma::vec grad_beta { beta }, strong_rhs { beta };
        l1_lambda_max_ =
            arma::max(arma::abs(gradient(beta)) / l1_penalty_factor_) /
            x_.n_rows / std::max(alpha_, 1e-3);

        // take unique lambda and sort descendingly
        lambda = arma::reverse(arma::unique(lambda));
        // construct lambda sequence
        arma::vec lambda_seq {
            arma::zeros(std::max(nlambda, lambda.n_elem))
        };
        if (nlambda > 1) {
            double log_lambda_max { std::log(l1_lambda_max_) };
            if (x_.n_cols > x_.n_rows) {
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

        // initialize the estimate matrix
        coef_mat_ = arma::zeros(x_.n_cols, lambda_seq.n_elem);
        en_coef_mat_ = coef_mat_;
        neg_ll_vec_ = arma::zeros(lambda_seq.n_elem);
        coef_df_vec_ = arma::zeros<arma::uvec>(lambda_seq.n_elem);
        bic_vec_ = neg_ll_vec_;

        // start values
        coef0_ = beta;
        coef_ = beta;

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x_.n_cols) };

        // optim with varying active set when p > n
        bool varying_active_set { false };
        if (x_.n_cols > x_.n_rows) {
            varying_active_set = true;
        }

        // outer loop for the lambda sequence
        for (size_t k {0}; k < lambda_seq.n_elem; ++k) {
            // early exit for large lambda greater than lambda_max
            if (alpha_ * lambda_seq(k) >= l1_lambda_max_) {
                // no re-scale is needed
                coef_mat_.col(k) = coef_;
                en_coef_mat_.col(k) = coef_;
                // compute negative log-likelihood
                neg_ll_vec_(k) = objective();
                coef_df_vec_(k) = compute_coef_df(beta);
                neg_ll_ = neg_ll_vec_(k);
                coef_df_ = coef_df_vec_(k);
                compute_bic();
                bic_vec_(k) = bic_;
                continue;
            }
            // allow users to stop the loop
            Rcpp::checkUserInterrupt();

            // update acitve set by strong rule (for lambda < lamda_max)
            grad_beta = arma::abs(gradient(beta)) / x_.n_rows;
            if (k == 0) {
                // use lambda_max
                strong_rhs = alpha_ *
                    (2 * lambda_seq(k) - l1_lambda_max_) * l1_penalty;
            } else {
                // use the last lambda
                strong_rhs = alpha_ *
                    (2 * lambda_seq(k) - lambda_seq(k - 1)) * l1_penalty;
            }
            for (size_t j {0}; j < x_.n_cols; ++j) {
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
                // allow users to stop the loop
                Rcpp::checkUserInterrupt();

                reg_active_fit(beta, is_active_strong, lambda_seq(k) * alpha_,
                               lambda_seq(k) * (1 - alpha_) / 2, l1_penalty,
                               varying_active_set, max_iter, rel_tol,
                               early_stop, verbose);
                // check kkt condition
                for (size_t j {0}; j < x_.n_cols; ++j) {
                    if (is_active_strong(j)) {
                        continue;
                    }
                    if (std::abs(gradient(beta, j)) / x_.n_rows >
                        strong_rhs(j)) {
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
            coef0_ = (1 + (1 - alpha_) * lambda_seq(k) / 2) * beta;
            rescale_coef();
            en_coef_mat_.col(k) = coef_;

            // compute naive elastic net estimates
            coef0_ = beta;
            rescale_coef();
            coef_mat_.col(k) = coef_;

            // compute negative log-likelihood
            neg_ll_vec_(k) = objective();
            coef_df_vec_(k) = compute_coef_df(beta);
            neg_ll_ = neg_ll_vec_(k);
            coef_df_ = coef_df_vec_(k);
            compute_bic();
            bic_vec_(k) = bic_;
        }
        // end of regularized fit
    }

}


#endif
