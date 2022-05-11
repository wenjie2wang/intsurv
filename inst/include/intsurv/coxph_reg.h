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
#include "control.h"
#include "utils.h"

namespace Intsurv {

    class CoxphReg
    {
    protected:
        // internals that depend on the data =================================
        double dn_obs_;         // double version of n_obs
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

        // function that computes objective function only
        inline double objective0(const arma::vec& beta) const;

        // function that computes gradients only
        inline arma::vec gradient0(const arma::vec& beta) const;
        inline double gradient0(const arma::vec& beta,
                                const unsigned int k) const;

        // function that computes objective and overwrite gradients
        inline double objective0(const arma::vec& beta,
                                 arma::vec& grad) const;

        // function computing B&L covariance lower bound matrix
        inline void compute_inv_bl_cov_lowerbound(const bool force_update);

        // function computing B&L lower bound for step size
        inline double bl_step_lowerbound(const arma::mat& x,
                                         const arma::vec& h) const;

        // get the lower bound of second derivative in CMD algorithm
        inline void set_cmd_lowerbound(const bool force_update);

        // starting values
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
            if (offset.n_elem == 1 || offset.empty()) {
                return arma::zeros(n_obs_);
            }
            throw std::length_error(
                "The length of the specified offset must match sample size."
                );
        }
        // generate and set penalty factor
        inline arma::vec gen_penalty_factor(
            const arma::vec& penalty_factor = arma::vec()
            ) const
        {
            if (penalty_factor.is_empty()) {
                return arma::ones(p_);
            }
            if (penalty_factor.n_elem == p_) {
                if (arma::any(penalty_factor < 0.0)) {
                    throw std::range_error(
                        "The 'penalty_factor' cannot be negative.");
                }
                return penalty_factor;
            }
            // else
            throw std::range_error("Incorrect length of the 'penalty_factor'.");
        }

        // make sure the following are well set
        inline CoxphReg* check_start()
        {
            if (control_.start_.n_elem != p_) {
                control_.start_ = gen_start(control_.start_);
            }
            return this;
        }
        inline CoxphReg* check_offset()
        {
            if (control_.offset_.n_elem != n_obs_) {
                control_.offset_ = gen_offset(control_.offset_);
            }
            return this;
        }
        inline CoxphReg* check_penalty_factor()
        {
            if (control_.penalty_factor_.n_elem != p_) {
                control_.penalty_factor_ = gen_penalty_factor(
                    control_.penalty_factor_);
            }
            return this;
        }

        // one-full cycle of coordinate-majorization-descent
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
        unsigned int n_obs_;    // number of observations
        unsigned int p_;        // number of given predictors
        arma::uvec ord_;        // index sorting the rows by time and event
        arma::uvec rev_ord_;    // index reverting the sorting
        // all the following data are sorted accordingly
        arma::mat x_;           // (sorted and standardized) design matrix
        arma::vec time_;        // (sorted) observed times
        arma::vec event_;       // (sorted) event indicators

        // controls
        Control control_;

        // more data
        arma::vec offset_haz_;     // offset term for baseline hazard only
        double l1_lambda_max_;     // the "big enough" lambda => zero coef
        double lambda_max_;

        // outputs ===========================================================
        arma::vec coef0_;       // coef for the standardized x
        // for a single control_.l1_lambda_ and control_.l2_lambda_
        arma::vec coef_;        // covariate coefficient estimates
        // for a lambda sequence
        arma::mat coef_mat_;     // coef_ matrix (rescaled for origin x_)
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
        arma::vec unique_time_;
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
                 const Control& control = Control()) :
            control_ (control)
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
            if (control_.standardize_) {
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
            if (control_.standardize_) {
                for (size_t j {0}; j < coef0_.n_elem; ++j) {
                    coef_[j] = coef0_[j] / x_scale_[j];
                }
            }
        }
        // re-scale all estimates for the original data
        inline void rescale_estimates()
        {
            rescale_coef();
            if (control_.standardize_) {
                double tmp { arma::as_scalar(x_center_ * coef_) };
                tmp = std::exp(tmp);
                h0_time_ /= tmp;
                H0_time_ /= tmp;
                S0_time_ = arma::exp(- H0_time_);
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

        inline CoxphReg* set_start(const arma::vec& start)
        {
            control_.start_ = gen_start(start);
            return this;
        }
        inline CoxphReg* set_start()
        {
            return set_start(control_.start_);
        }

        inline CoxphReg* set_offset(const arma::vec& offset,
                                    const bool is_sorted = true)
        {
            control_.offset_ = gen_offset(offset);
            if (! is_sorted) {
                // update control_.offset_ for appropriate input
                control_.offset_ = control_.offset_.elem(ord_);
            }
            // update d_offset_ as well
            d_offset_ = control_.offset_.elem(event_ind_) %
                event_.elem(event_ind_);
            if (has_ties_) {
                d_offset_ = aggregate_sum(d_offset_, d_time0_);
            }
            return this;
        }
        inline CoxphReg* set_offset()
        {
            return set_offset(control_.offset_, false);
        }

        inline CoxphReg* reset_offset()
        {
            control_.offset_ = arma::zeros(n_obs_);
            if (has_ties_) {
                d_offset_ = arma::zeros(d_time_.n_elem);
            } else {
                d_offset_ = arma::zeros(event_ind_.n_elem);
            }
            return this;
        }

        inline CoxphReg* set_penalty_factor(const arma::vec& penalty_factor)
        {
            control_.penalty_factor_ = gen_penalty_factor(penalty_factor);
            return this;
        }
        inline CoxphReg* set_penalty_factor()
        {
            return set_penalty_factor(control_.penalty_factor_);
        }

        // set offset for denominator in baseline hazard function
        // for cure rate model
        inline CoxphReg* set_offset_haz(const arma::vec& offset,
                                        const bool is_sorted = true)
        {
            if (offset.n_elem == n_obs_) {
                offset_haz_ = offset;
            } else if (offset.n_elem == 1 || offset.empty()) {
                return reset_offset_haz();
            } else {
                throw std::length_error(
                    "The length of offset must match sample size.");
            }
            if (! is_sorted) {
                // update control_.offset_ for appropriate input
                offset_haz_ = offset_haz_.elem(ord_);
            }
            // update d_offset_haz as well
            d_offset_haz_ = offset_haz_.elem(event_ind_) %
                event_.elem(event_ind_);
            if (has_ties_) {
                d_offset_haz_ = aggregate_sum(d_offset_haz_, d_time0_);
            }
            return this;
        }

        inline CoxphReg* reset_offset_haz()
        {
            offset_haz_ = arma::zeros(n_obs_);
            if (has_ties_) {
                d_offset_haz_ = arma::zeros(d_time_.n_elem);
            } else {
                d_offset_haz_ = arma::zeros(event_ind_.n_elem);
            }
            return this;
        }

        // function that computes baseline estimates
        inline void compute_haz_surv_time(const arma::vec& beta);
        inline void compute_haz_surv_time();
        inline void compute_censor_haz_surv_time();

        // prepare hazard and survival estimates on unique time points
        inline void est_haz_surv();

        // additional methods for coxph_cure
        // revserse the rescale process to get coef0_ from a new coef_
        inline void rev_rescale_coef()
        {
            coef0_ = coef_;
            if (control_.standardize_) {
                for (size_t j {0}; j < coef_.n_elem; ++j) {
                    coef0_[j] = coef_[j] * x_scale_[j];
                }
            }
        }

        inline arma::vec rev_rescale_coef(const arma::vec& beta) const
        {
            if (control_.standardize_) {
                arma::vec out { beta };
                for (size_t j {0}; j < out.n_elem; ++j) {
                    out[j] = out[j] * x_scale_[j];
                }
                return out;
            }
            return beta;
        }

        inline double objective() const
        {
            return objective0(coef0_);
        }

        inline double objective(const arma::vec& beta) const
        {
            arma::vec beta0 { rev_rescale_coef(beta) };
            return objective0(beta0);
        }

        inline double net_penalty(const arma::vec& beta,
                                  const double l1_lambda,
                                  const double l2_lambda,
                                  const arma::vec& penalty_factor) const
        {
            return l1_lambda * l1_norm(beta % penalty_factor) +
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

        inline void set_l1_lambda_max()
        {
            l1_lambda_max_ = get_l1_lambda_max(control_.penalty_factor_);
        }

        inline arma::vec get_xbeta(const arma::vec& beta) const
        {
            if (control_.standardize_) {
                arma::vec xbeta;
                // re-scale the input beta
                arma::vec beta0 { beta % x_scale_.t() };
                xbeta = x_ * beta0 + arma::as_scalar(x_center_ * beta);
                return xbeta;
            }
            return mat2vec(x_ * beta);
        }

        inline arma::vec get_xbeta() const
        {
            return get_xbeta(coef_);
        }

        // fit regular Cox model
        inline void fit();

        // fit regularized Cox model with adaptive lasso penalty
        // for a perticular lambda
        inline void net_fit();

        // for a sequence of lambda's
        inline void net_path();

    };
    // end of class definition

    // compute baseline hazard function and its friends
    // here beta is the estimate for non-standardized data
    inline void CoxphReg::compute_haz_surv_time(const arma::vec& beta)
    {
        arma::vec xbeta { get_xbeta(beta) };
        arma::vec exp_risk_score { arma::exp(xbeta + control_.offset_) };
        // 1. hazard rate function
        arma::vec h0_numer { aggregate_sum(event_, time_, false) };
        arma::vec h0_denom { exp_risk_score % arma::exp(offset_haz_) };
        h0_denom = aggregate_sum(h0_denom, time_, false, true, true);
        h0_time_ = h0_numer / h0_denom;
        // set h0_time to be zero if h0_denom is zero
        h0_time_.replace(arma::datum::nan, 0.0);
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
    inline double CoxphReg::objective0(const arma::vec& beta) const
    {
        const arma::vec dx_beta {
            d_x_ * beta + d_offset_ + d_offset_haz_
        };
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + control_.offset_ + offset_haz_)
        };
        const arma::vec h0_denom {
            cum_sum(exp_x_beta, true)
        };
        arma::vec log_h0_denom_event {
            arma::log(h0_denom.elem(uni_event_ind_))
        };
        return - arma::sum(dx_beta - delta_n_ % log_h0_denom_event) / dn_obs_;
    }
    // the gradient of negative loglikelihood function
    inline arma::vec CoxphReg::gradient0(const arma::vec& beta) const
    {
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + control_.offset_ + offset_haz_)
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
        return - (arma::sum(d_x_, 0) - arma::sum(numer_mat, 0)).t() / dn_obs_;
    }
    // the gradient of negative loglikelihood function at k-th direction
    inline double CoxphReg::gradient0(const arma::vec& beta,
                                      const unsigned int k) const
    {
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + control_.offset_ + offset_haz_)
        };
        arma::vec h0_denom { cum_sum(exp_x_beta, true) };
        arma::vec numer { cum_sum(mat2vec(x_.col(k) % exp_x_beta), true) };
        h0_denom = h0_denom.elem(uni_event_ind_);
        numer = numer.elem(uni_event_ind_);
        return - arma::sum(d_x_.col(k) - delta_n_ % numer / h0_denom) / dn_obs_;
    }

    // the negative log-likelihood function based on the broslow's formula
    inline double CoxphReg::objective0(const arma::vec& beta,
                                       arma::vec& grad) const
    {
        const arma::vec dx_beta {
            d_x_ * beta + d_offset_ + d_offset_haz_
        };
        const arma::vec exp_x_beta {
            arma::exp(x_ * beta + control_.offset_ + offset_haz_)
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
        grad = - (arma::sum(d_x_, 0) - arma::sum(numer_mat, 0)).t() / dn_obs_;
        return - arma::sum(dx_beta - delta_n_ % arma::log(h0_denom_event)) /
            dn_obs_;
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
    inline void CoxphReg::set_cmd_lowerbound(
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
    inline void CoxphReg::fit()
    {
        compute_inv_bl_cov_lowerbound();
        check_start()->check_offset();
        arma::vec beta0 { control_.start_ };
        arma::vec beta { beta0 }, h_vec { beta0 }, grad_vec { beta0 };
        double b_new {0.0}, alpha_ {0.0};
        double ell { arma::datum::inf };
        if (control_.verbose_ > 1) {
            Rcpp::Rcout << "\n" << std::string(40, '=')
                        << "\nStarting from\n"
                        << arma2rvec(beta0)
                        << "\n";
        }
        if (control_.verbose_ > 0) {
            ell = objective0(beta0);
        }
        for (size_t i {0}; i < control_.max_iter_; ++i) {
            // allow users to stop the main loop
            Rcpp::checkUserInterrupt();
            // compute negative log-likelihood function and update gradient
            grad_vec = gradient0(beta0) * dn_obs_;
            h_vec = - inv_bl_cov_lowerbound_ * grad_vec;
            b_new = bl_step_lowerbound(x_, h_vec);
            alpha_ = - arma::as_scalar(crossprod(h_vec, grad_vec)) / b_new;
            beta = beta0 + alpha_ * h_vec;
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
            // if tolerance is reached
            if (rel_l1_norm(beta, beta0) < control_.epsilon_) {
                if (control_.verbose_ > 0) {
                    Rcpp::Rcout << "\nReached convergence criterion\n";
                }
                break;
            }
            beta0 = beta;
        }
        coef0_ = beta;
        rescale_coef();
    }

    // run one cycle of coordinate descent over a given active set
    inline void CoxphReg::net_one_update(
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
                        << arma2rvec(beta)
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
        for (size_t j {0}; j < beta.n_elem; ++j) {
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
                    cmd_lowerbound_[j] + 2 * l2_lambda
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
    inline void CoxphReg::net_active_update(
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
                size_t ii {0};
                while (ii < max_iter) {
                    net_one_update(beta, is_active_varying,
                                   l1_lambda, l2_lambda, penalty_factor,
                                   true, verbose);
                    if (rel_l1_norm(beta, beta0) <epsilon) {
                        num_iter = ii + 1;
                        break;
                    }
                    beta0 = beta;
                    ii++;
                }
                // run a full cycle over the converged beta
                net_one_update(beta, is_active,
                               l1_lambda, l2_lambda, penalty_factor,
                               true, verbose);
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
            for (size_t i {0}; i < control_.max_iter_; ++i) {
                net_one_update(beta, is_active,
                               l1_lambda, l2_lambda, penalty_factor,
                               false, verbose);
                if (rel_l1_norm(beta, beta0) < epsilon) {
                    num_iter = i + 1;
                    break;
                }
                beta0 = beta;
            }
        }
        if (control_.verbose_ > 0) {
            if (num_iter < control_.max_iter_) {
                Rcpp::Rcout << "Converged after "
                            << num_iter
                            << " iteration(s)\n";
            } else {
                msg("Reached the maximum number of iteratons.");
            }
        }
    }

    // fitting regularized Cox model with coordinate-majorizatio-descent
    // algorithm that allows non-integer "event_" and tied events
    // for a perticular control_.l1_lambda_ and control_.l2_lambda_
    // lambda_1 * lasso + lambda_2 * ridge
    inline void CoxphReg::net_fit()
    {
        set_cmd_lowerbound();
        check_start()->check_offset()->check_penalty_factor();
        // use the given starting values
        arma::vec beta { control_.start_ };
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
    inline void CoxphReg::net_path()
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
                control_.lambda_ = arma::reverse(arma::unique(control_.lambda_));
            }
        }
        // initialize the estimate matrix
        coef_mat_ = arma::zeros(p_, control_.lambda_.n_elem);
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(p_) };
        // for ridge penalty
        if (is_ridge_only) {
            is_active_strong = arma::ones<arma::uvec>(x_.n_cols);
            for (size_t li { 0 }; li < control_.lambda_.n_elem; ++li) {
                net_active_update(beta, is_active_strong,
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
