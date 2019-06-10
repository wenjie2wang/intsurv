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

#ifndef COXPH_REG_H
#define COXPH_REG_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {

    // define class for inputs and outputs
    class CoxphReg {
    private:
        arma::uvec ord;            // index sorting the rows
        arma::uvec rev_ord;        // index reverting the sorting

        // all the following results are sorted accordingly
        arma::vec time;            // (sorted) observed times
        arma::vec event;           // (sorted) event indicators
        arma::mat x;               // (sorted and standardized) design matrix

        bool standardize;          // is x standardized
        arma::rowvec x_center;     // the column center of x
        arma::rowvec x_scale;      // the scale of x
        arma::vec coef0;           // coef before scaling

        bool hasTies {false};      // if there exists ties on event times

        // at each unique time point
        arma::uvec event_ind;      // indices of event times
        arma::uvec uni_event_ind;  // the index indicating the first record
                                   // on each distinct event time
        arma::uvec uni_time_ind;   // the index indicating the first record on
                                   // each time whether event or censoring
        arma::vec offset;          // offset term
        // size of risk-set at each time point
        arma::vec riskset_size_time;

        // at each unique event time point
        arma::vec d_time;          // distinct event times
        arma::mat d_x;             // design matrix aggregated at d_time
        arma::vec d_offset;        // offset terms aggregated at d_time
        arma::vec delta_n;         // event counts at d_time
        arma::vec riskset_size;    // size of risk-set at d_time

        // inverse matrix of formula (5.9) in Bohning and Lindsay (1988) SIAM
        arma::mat inv_bl_cov_lowerbound_1;
        // for regularized regression
        arma::vec cmd_lowerbound;

    public:
        unsigned int nObs;           // number of observations
        arma::vec l1_penalty_factor; // adaptive weights for lasso penalty
        double l1_lambda_max;        // the "big enough" lambda => zero coef

        // for a single l1_lambda and l2_lambda
        double l1_lambda;            // tuning parameter for lasso penalty
        double l2_lambda;            // tuning parameter for ridge penalty
        arma::vec coef;              // covariate coefficient estimates
        arma::vec en_coef;           // (rescaled) elastic net estimates
        double negLogL;              // partial negative log-likelihood
        unsigned int coef_df;        // number of non-zero coef estimates
        arma::vec xBeta;             // x * coef

        // for a lambda sequence
        double alpha;           // tuning parameter
        arma::vec lambda_vec;   // lambda sequence
        arma::mat coef_mat;     // coef matrix (rescaled for origin x)
        arma::mat en_coef_mat;  // elastic net estimates
        arma::vec negLogL_vec;  // negative log-likelihood vector
        arma::uvec coef_df_vec; // coef df vector

        // hazard and survival estimates at every time point (unique or not)
        arma::vec h0_time;
        arma::vec S0_time;
        arma::vec H0_time;
        arma::vec h_time;
        arma::vec H_time;
        arma::vec S_time;
        arma::vec hc_time;
        arma::vec Hc_time;
        arma::vec Sc_time;

        // baseline estimates at unique times only
        arma::vec unique_time;
        arma::vec h0_est;
        arma::vec H0_est;
        arma::vec S0_est;
        arma::vec hc_est;
        arma::vec Hc_est;
        arma::vec Sc_est;

        // default constructor
        CoxphReg() {}

        // constructors
        CoxphReg(const arma::vec& time_,
                 const arma::vec& event_,
                 const arma::mat& x_,
                 const bool& standardize_ = true)
        {
            // sort based on time and event
            // time: ascending order
            // event: events first, then censoring at the same time point
            arma::uvec des_event_ind { arma::sort_index(event_, "descend") };
            arma::uvec asc_time_ind {
                arma::stable_sort_index(time_.elem(des_event_ind), "ascend")
            };
            ord = des_event_ind.elem(asc_time_ind);
            rev_ord = arma::sort_index(ord);

            // do actual sorting
            time = time_.elem(ord);
            event = event_.elem(ord);
            x = x_.rows(ord);
            this->nObs = x.n_rows;

            // standardize covariates
            standardize = standardize_;
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

            // binary event indicator
            event_ind = arma::find(event > 0);
            // check if there exists ties on event times
            arma::vec d_time0 {time.elem(event_ind)};
            hasTies = any_duplicated(d_time0);
            // default value when no tied event times
            uni_event_ind = event_ind;
            d_time = d_time0;
            delta_n = event.elem(event_ind);
            d_x = x.rows(event_ind);
            for (size_t j {0}; j < x.n_cols; ++j) {
                d_x.col(j) = d_x.col(j) % delta_n;
            }
            // initialize offset
            this->set_offset(arma::zeros(1));
            this->d_offset = this->offset.elem(event_ind) % delta_n;

            uni_time_ind = find_first_unique(time);
            if (hasTies) {
                // re-define uni_event_ind
                uni_event_ind = vec_intersection(uni_time_ind, event_ind);
                d_time = time.elem(uni_event_ind);
                // aggregate at distinct event times
                delta_n = aggregate_sum(delta_n, d_time0);
                d_x = aggregate_sum(d_x, d_time0);
                d_offset = aggregate_sum(d_offset, d_time0);
            }
            riskset_size = arma::ones(time.n_elem);
            riskset_size = cum_sum(riskset_size, true).elem(uni_event_ind);
        }

        // re-scale coef back
        inline void rescale_coef()
        {
            this->coef = coef0;
            if (this->standardize) {
                for (size_t j {0}; j < coef0.n_elem; ++j) {
                    this->coef[j] = coef0[j] / x_scale[j];
                }
            }
        }

        // partially update event non-zero-one indicators
        // assuming only old event weights strictly between 0 and 1 need update
        inline void update_event_weight(const arma::vec& event_,
                                        const bool& is_sorted = true)
        {
            event = event_;
            if (! is_sorted) {
                event = event.elem(ord);
            }
            // need to update delta_n and d_x
            delta_n = event.elem(event_ind);
            d_x = x.rows(event_ind);
            for (size_t j {0}; j < x.n_cols; ++j) {
                d_x.col(j) = d_x.col(j) % delta_n;
            }
            if (hasTies) {
                arma::vec d_time0 {time.elem(event_ind)};
                // aggregate at distinct event times
                delta_n = aggregate_sum(delta_n, d_time0);
                d_x = aggregate_sum(d_x, d_time0);
            }
        }

        // set offset
        inline void set_offset(const arma::vec& offset_,
                               const bool& is_sorted = true)
        {
            if (offset_.n_elem == x.n_rows) {
                offset = offset_;
                if (! is_sorted) {
                    // update offset for appropriate input
                    offset = offset.elem(ord);
                }
            } else {
                offset = arma::zeros(time.n_elem);
            }
        }

        // function that computes baseline estimates
        inline void compute_haz_surv_time(const arma::vec& beta);
        inline void compute_haz_surv_time();
        inline void compute_censor_haz_surv_time();

        // prepare hazard and survival estimates on unique time points
        inline void est_haz_surv();

        // function that computes objective function only
        inline double objective() const;
        inline double objective(const arma::vec& beta) const;

        // function that computes gradients only
        inline arma::vec gradient(const arma::vec& beta) const;
        inline double gradient(const arma::vec& beta,
                               const unsigned int& k) const;

        // function that computes objective and overwrite gradients
        inline double objective(const arma::vec& beta,
                                arma::vec& grad) const;

        // function computing B&L covariance lower bound matrix
        inline void compute_inv_bl_cov_lowerbound_1(const bool& force_update);

        // function computing B&L lower bound for step size
        inline double bl_step_lowerbound(const arma::mat& x,
                                         const arma::vec& h) const;

        // get the lower bound of second derivative in CMD algorithm
        inline void compute_cmd_lowerbound(const bool& force_update);

        // some simple functions
        inline unsigned int sample_size() const { return time.n_elem; }

        // helper function to access some private members
        inline arma::vec get_time() const { return time; }
        inline arma::vec get_event() const { return event; }
        inline arma::mat get_x() const { return x; }
        inline arma::uvec get_sort_index() { return this->ord; }

        // fit regular Cox model
        inline void fit(const arma::vec& start,
                        const unsigned int& max_iter,
                        const double& rel_tol,
                        const bool& early_stop,
                        const bool& verbose);

        // one-full cycle of coordinate-majorization-descent
        inline void regularized_fit_update(arma::vec& beta,
                                           arma::uvec& is_active,
                                           const double& l1_lambda,
                                           const double& l2_lambda,
                                           const arma::vec& l1_penalty_factor,
                                           const bool& update_active,
                                           const bool& early_stop,
                                           const bool& verbose
            );

        inline void reg_active_fit(arma::vec& beta,
                                   const arma::uvec& is_active,
                                   const double& l1_lambda,
                                   const double& l2_lambda,
                                   const arma::vec& l1_penalty_factor,
                                   const bool& varying_active_set,
                                   const unsigned int& max_iter,
                                   const double& rel_tol,
                                   const bool& early_stop,
                                   const bool& verbose
            );

        // fit regularized Cox model with adaptive lasso penalty
        // for a perticular lambda
        inline void regularized_fit(const double& l1_lambda,
                                    const double& l2_lambda,
                                    const arma::vec& l1_penalty_factor,
                                    const arma::vec& start,
                                    const unsigned int& max_iter,
                                    const double& rel_tol,
                                    const bool& early_stop,
                                    const bool& verbose);

        // overload for a sequence of lambda's
        inline void regularized_fit(arma::vec lambda,
                                    const double& alpha,
                                    const unsigned int& nlambda,
                                    double lambda_min_ratio,
                                    const arma::vec& l1_penalty_factor,
                                    const unsigned int& max_iter,
                                    const double& rel_tol,
                                    const bool& early_stop,
                                    const bool& verbose);

    };
    // end of class definition

    // compute baseline hazard function and its friends
    // here beta is coef estimate for non-standardized data
    inline void CoxphReg::compute_haz_surv_time(const arma::vec& beta)
    {
        if (this->standardize) {
            if (! beta.is_empty() && beta.n_elem == x.n_cols) {
                // re-scale the input beta
                arma::vec beta0 { beta % x_scale.t() };
                this->xBeta = x * beta0 + arma::as_scalar(x_center * beta);
            } else {
                // use estimated coef
                this->xBeta = x * this->coef0 +
                    arma::as_scalar(x_center * this->coef);
            }
        } else {
            if (! beta.is_empty() && beta.n_elem == x.n_cols) {
                this->xBeta = x * beta;
            } else {
                // use estimated coef
                this->xBeta = x * this->coef;
            }
        }
        arma::vec exp_x_beta { arma::exp(this->xBeta) };

        // 1. hazard rate function
        arma::vec h0_numer { aggregate_sum(event, time, false) };
        arma::vec h0_denom { exp_x_beta % arma::exp(this->offset) };
        h0_denom = aggregate_sum(h0_denom, time, false, true, true);
        this->h0_time = h0_numer / h0_denom;
        this->h_time = this->h0_time % exp_x_beta;

        // 2. baseline cumulative hazard function
        this->H0_time = arma::zeros(h0_time.n_elem);
        for (size_t i: uni_event_ind) {
            this->H0_time(i) = h0_time(i);
        }
        this->H0_time = cum_sum(this->H0_time);
        this->H_time = this->H0_time % exp_x_beta;

        // 3. baseline survival function
        this->S0_time = arma::exp(- this->H0_time);
        this->S_time = arma::exp(- this->H_time);
    }
    inline void CoxphReg::compute_haz_surv_time()
    {
        compute_haz_surv_time(this->coef);
    }
    inline void CoxphReg::compute_censor_haz_surv_time()
    {
        arma::vec censor_ind { 1 - event };
        arma::vec delta_c { aggregate_sum(censor_ind, time, false) };
        if (this->riskset_size_time.is_empty()) {
            arma::vec tmp { arma::ones(this->time.n_elem) };
            this->riskset_size_time =
                aggregate_sum(tmp, time, false, true, true);
        }
        this->hc_time = delta_c / this->riskset_size_time;
        this->Hc_time = arma::zeros(this->hc_time.n_elem);
        for (size_t i: uni_time_ind) {
            this->Hc_time(i) = this->hc_time(i);
        }
        this->Hc_time = cum_sum(this->Hc_time);
        this->Sc_time = arma::exp(- this->Hc_time);
    }

    // prepare hazard and survival estimates on unique time points
    // should be used after compute_haz_surv_time() and
    // compute_censor_haz_surv_time()
    inline void CoxphReg::est_haz_surv()
    {
        this->unique_time = this->time.elem(this->uni_time_ind);
        if (this->h0_time.is_empty()) {
            this->compute_haz_surv_time();
        }
        if (this->hc_time.is_empty()) {
            this->compute_censor_haz_surv_time();
        }
        this->h0_est = this->h0_time.elem(this->uni_time_ind);
        this->H0_est = this->H0_time.elem(this->uni_time_ind);
        this->S0_est = this->S0_time.elem(this->uni_time_ind);
        this->hc_est = this->hc_time.elem(this->uni_time_ind);
        this->Hc_est = this->Hc_time.elem(this->uni_time_ind);
        this->Sc_est = this->Sc_time.elem(this->uni_time_ind);
    }

    // the negative log-likelihood function based on the broslow's formula
    inline double CoxphReg::objective(const arma::vec& beta) const
    {
        const arma::vec dx_beta {d_x * beta + d_offset};
        const arma::vec exp_x_beta {arma::exp(x * beta + offset)};
        const arma::vec h0_denom {cum_sum(exp_x_beta, true)};
        arma::vec log_h0_denom_event {arma::log(h0_denom.elem(uni_event_ind))};
        return - arma::sum(dx_beta - delta_n % log_h0_denom_event);
    }
    inline double CoxphReg::objective() const
    {
        return this->objective(this->coef0);
    }

    // the gradient of negative loglikelihood function
    inline arma::vec CoxphReg::gradient(const arma::vec& beta) const
    {
        const arma::vec exp_x_beta {arma::exp(x * beta + offset)};
        const arma::vec h0_denom {cum_sum(exp_x_beta, true)};
        arma::mat numer_mat {arma::zeros(x.n_rows, x.n_cols)};
        for (size_t i {0}; i < x.n_rows; ++i) {
            numer_mat.row(i) = exp_x_beta(i) * x.row(i);
        }
        numer_mat = cum_sum(numer_mat, true);
        arma::vec h0_denom_event {h0_denom.elem(uni_event_ind)};
        numer_mat = numer_mat.rows(uni_event_ind);
        for (size_t j {0}; j < x.n_cols; ++j) {
            numer_mat.col(j) = numer_mat.col(j) % delta_n / h0_denom_event;
        }
        return - (arma::sum(d_x, 0) - arma::sum(numer_mat, 0)).t();
    }
    // the gradient of negative loglikelihood function at k-th direction
    inline double CoxphReg::gradient(const arma::vec& beta,
                                     const unsigned int& k) const
    {
        const arma::vec exp_x_beta { arma::exp(x * beta + offset) };
        arma::vec h0_denom { cum_sum(exp_x_beta, true) };
        arma::vec numer { cum_sum(mat2vec(x.col(k) % exp_x_beta), true) };
        h0_denom = h0_denom.elem(uni_event_ind);
        numer = numer.elem(uni_event_ind);
        return - arma::sum(d_x.col(k) - delta_n % numer / h0_denom);
    }

    // the negative log-likelihood function based on the broslow's formula
    inline double CoxphReg::objective(const arma::vec& beta,
                                      arma::vec& grad) const
    {
        const arma::vec dx_beta {d_x * beta + d_offset};
        const arma::vec exp_x_beta {arma::exp(x * beta + offset)};
        const arma::vec h0_denom {cum_sum(exp_x_beta, true)};
        arma::mat numer_mat {arma::zeros(x.n_rows, x.n_cols)};
        for (size_t i {0}; i < x.n_rows; ++i) {
            numer_mat.row(i) = exp_x_beta(i) * x.row(i);
        }
        numer_mat = cum_sum(numer_mat, true);
        arma::vec h0_denom_event {h0_denom.elem(uni_event_ind)};
        numer_mat = numer_mat.rows(uni_event_ind);
        for (size_t j {0}; j < x.n_cols; ++j) {
            numer_mat.col(j) = numer_mat.col(j) % delta_n / h0_denom_event;
        }
        // overwrite grad
        grad = - (arma::sum(d_x, 0) - arma::sum(numer_mat, 0)).t();
        return - arma::sum(dx_beta - delta_n % arma::log(h0_denom_event));
    }

    // one lower bound of covariance matrix
    // reference: formula (5.9) in Bohning and Lindsay (1988) SIAM
    inline void CoxphReg::compute_inv_bl_cov_lowerbound_1(
        const bool& force_update = false
        )
    {
        if (inv_bl_cov_lowerbound_1.is_empty() || force_update) {
            arma::mat res { arma::zeros(x.n_cols, x.n_cols) };
            // compute x * x^t and store them in a cube
            arma::cube sum_risk_xxt {
                arma::cube(x.n_cols, x.n_cols, time.n_elem, arma::fill::zeros)
            };
            for (size_t i {0}; i < time.n_elem; ++i) {
                sum_risk_xxt.slice(i) = crossprod(x.row(i));
            }
            sum_risk_xxt = cum_sum(sum_risk_xxt, true);
            sum_risk_xxt = sum_risk_xxt.slices(uni_event_ind);
            arma::mat sum_risk_x { cum_sum(x, true) };
            sum_risk_x = sum_risk_x.rows(uni_event_ind);
            for (unsigned int i {0}; i < d_time.n_elem; ++i) {
                res += (sum_risk_xxt.slice(i) - crossprod(sum_risk_x.row(i)) /
                        riskset_size(i)) * (delta_n(i) / 2);
            }
            this->inv_bl_cov_lowerbound_1 = arma::inv_sympd(res);
        }
    }
    inline double CoxphReg::bl_step_lowerbound(const arma::mat& x,
                                               const arma::vec& h) const
    {
        // compute x * h and store them in a vector
        arma::vec htx { x * h };
        arma::vec max_risk_htx { cum_max(htx, true) };
        max_risk_htx = max_risk_htx.elem(uni_event_ind);
        arma::vec min_risk_htx { cum_min(htx, true) };
        min_risk_htx = min_risk_htx.elem(uni_event_ind);
        double Mi {0}, mi {0}, res {0};
        for (unsigned int i {0}; i < d_time.n_elem; ++i) {
            Mi = max_risk_htx(i);
            mi = min_risk_htx(i);
            res += ((std::pow(Mi, 2) + std::pow(mi, 2)) / 4 - Mi * mi / 2) *
                delta_n(i);
        }
        return res;
    }

    // compute lower bound of second derivative
    // in coordinate-majorization-descent algorithm (cmd)
    inline void CoxphReg::compute_cmd_lowerbound(
        const bool& force_update = false
        )
    {
        if (this->cmd_lowerbound.is_empty() || force_update) {
            // compute D_k at each event time point
            arma::mat max_risk_x { cum_max(x, true) };
            arma::mat min_risk_x { cum_min(x, true) };
            max_risk_x = max_risk_x.rows(uni_event_ind);
            min_risk_x = min_risk_x.rows(uni_event_ind);
            arma::vec res { arma::zeros(x.n_cols) };
            for (size_t k {0}; k < x.n_cols; ++k) {
                res(k) = arma::sum(
                    pow(max_risk_x.col(k) - min_risk_x.col(k), 2) % delta_n
                    ) / (4 * x.n_rows);
            }
            this->cmd_lowerbound = res;
        }
    }

    // fit regular Cox model by monotonic quadratic approximation algorithm
    // that allows non-integer "event" and tied events
    // reference: Bohning and Lindsay (1988) SIAM
    inline void CoxphReg::fit(const arma::vec& start = 0,
                              const unsigned int& max_iter = 100,
                              const double& rel_tol = 1e-6,
                              const bool& early_stop = true,
                              const bool& verbose = false)
    {
        this->compute_inv_bl_cov_lowerbound_1();
        arma::vec beta0 { arma::zeros(x.n_cols) };
        if (start.n_elem == x.n_cols) {
            beta0 = start;
        }
        arma::vec beta { beta0 }, h_vec { beta0 }, grad_vec { beta0 };
        size_t i {0};
        double b_new {0}, alpha {0}, ell {0}, ell_old {arma::datum::inf};
        while (i <= max_iter) {
            // compute negative log-likelihood function and update gradient
            ell = this->objective(beta, grad_vec);
            if (verbose) {
                Rcpp::Rcout << "\n" << std::string(40, '=')
                            << "\niteration: " << i
                            << "\n  estimated coef: "
                            << arma2rvec(beta)
                            << "\n  negative log-likelihood: "
                            << ell
                            << std::endl;
            }
            h_vec = - this->inv_bl_cov_lowerbound_1 * grad_vec;
            b_new = this->bl_step_lowerbound(x, h_vec);
            alpha = - arma::as_scalar(crossprod(h_vec, grad_vec)) / b_new;
            beta = beta0 + alpha * h_vec;

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
                                << "The negative log-likelihood increased"
                                << std::endl;
                    Rprintf("  from %15.15f\n", ell_old);
                    Rprintf("    to %15.15f\n", ell);
                }
                if (early_stop) {
                    if (verbose) {
                        Rcpp::Rcout << "\nEarly stopped the algorithm"
                                    << " with estimates from"
                                    << " iteration " << i - 1 << "."
                                    << std::endl;
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
        this->coef0 = beta;
        this->rescale_coef();
        this->negLogL = ell;
        this->coef_df = beta.n_elem;
    }

    // run one cycle of coordinate descent over a given active set
    inline void CoxphReg::regularized_fit_update(
        arma::vec& beta,
        arma::uvec& is_active,
        const double& l1_lambda,
        const double& l2_lambda,
        const arma::vec& l1_penalty_factor,
        const bool& update_active = false,
        const bool& early_stop = false,
        const bool& verbose = false
        )
    {
        // compute lowerbound vector used in CMD algorithm
        compute_cmd_lowerbound();
        double dlj { 0 };
        arma::vec beta_old = beta;
        for (size_t j {0}; j < beta.n_elem; ++j) {
            if (is_active[j]) {
                dlj = gradient(beta, j) / time.n_elem;
                // update beta
                beta[j] = soft_threshold(
                    this->cmd_lowerbound[j] * beta[j] - dlj,
                    l1_penalty_factor[j] * l1_lambda
                    ) /
                    (this->cmd_lowerbound[j] + 2 * l2_lambda);
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
            double ell_old { this->objective(beta_old) };
            ell_old = ell_old / this->nObs +
                l1_lambda * l1_norm(beta_old % l1_penalty_factor) +
                l2_lambda * sum_of_square(beta_old);
            double ell_new { this->objective(beta) };
            ell_new = ell_new / this->nObs +
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
                                << "the objective function somehow increased\n";
                    Rcpp::Rcout << "\nEarly stopped the CMD iterations "
                                << "with estimates from the last step"
                                << std::endl;
                }
                beta = beta_old;
            }
        }
    }

    // run a complete cycle of CMD for a given active set and lambda
    inline void CoxphReg::reg_active_fit(
        arma::vec& beta,
        const arma::uvec& is_active,
        const double& l1_lambda,
        const double& l2_lambda,
        const arma::vec& l1_penalty_factor,
        const bool& varying_active_set = false,
        const unsigned int& max_iter = 300,
        const double& rel_tol = 1e-5,
        const bool& early_stop = false,
        const bool& verbose = false
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
    // algorithm that allows non-integer "event" and tied events
    // for a perticular l1_lambda and l2_lambda
    // lambda_1 * lasso + lambda_2 * ridge
    inline void CoxphReg::regularized_fit(
        const double& l1_lambda = 0,
        const double& l2_lambda = 0,
        const arma::vec& l1_penalty_factor = 0,
        const arma::vec& start = 0,
        const unsigned int& max_iter = 300,
        const double& rel_tol = 1e-5,
        const bool& early_stop = false,
        const bool& verbose = false
        )
    {
        // set penalty terms
        arma::vec l1_penalty { arma::ones(x.n_cols) };
        if (l1_penalty_factor.n_elem == x.n_cols) {
            // re-scale so that sum(factor) = number of predictors
            l1_penalty = l1_penalty_factor * x.n_cols /
                arma::sum(l1_penalty_factor);
        }
        this->l1_penalty_factor = l1_penalty;

        // the maximum (large enough) lambda that results in all-zero estimates
        arma::vec beta { arma::zeros(x.n_cols) };
        arma::vec grad_beta { beta }, strong_rhs { beta };
        this->l1_lambda_max =
            arma::max(arma::abs(this->gradient(beta)) /
                      l1_penalty) / this->x.n_rows;

        // early exit for large lambda greater than lambda_max
        this->coef0 = beta;
        this->coef = beta;
        if (l1_lambda > this->l1_lambda_max) {
            // no need to rescale all-zero coef
            this->en_coef = this->coef;
            this->coef_df = 0;
            return;
        }

        // use the input starting value
        if (start.n_elem == x.n_cols) {
            beta = start;
        }

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x.n_cols) };
        // update active set by strong rule
        grad_beta = arma::abs(this->gradient(this->coef0)) / this->x.n_rows;
        strong_rhs = (2 * l1_lambda - this->l1_lambda_max) * l1_penalty;
        for (size_t j {0}; j < x.n_cols; ++j) {
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

        bool kkt_failed { true };
        strong_rhs = l1_lambda * l1_penalty;
        // eventually, strong rule will guess correctly
        while (kkt_failed) {
            // update beta
            reg_active_fit(beta, is_active_strong, l1_lambda, l2_lambda,
                           l1_penalty, varying_active_set, max_iter, rel_tol,
                           early_stop, verbose);
            // check kkt condition
            for (size_t j {0}; j < x.n_cols; ++j) {
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
        // record other inputs
        this->l1_lambda = l1_lambda;
        this->l2_lambda = l2_lambda;
        this->coef_df = get_coef_df(beta);
    }


    // for a sequence of lambda's
    // lambda * (penalty_factor * alpha * lasso + (1 - alpha) / 2 * ridge)
    inline void CoxphReg::regularized_fit(
        arma::vec lambda = 0,
        const double& alpha = 1,
        const unsigned int& nlambda = 1,
        double lambda_min_ratio = 1e-4,
        const arma::vec& l1_penalty_factor = 0,
        const unsigned int& max_iter = 300,
        const double& rel_tol = 1e-5,
        const bool& early_stop = false,
        const bool& verbose = false
        )
    {
        // check alpha
        if ((alpha < 0) || (alpha > 1)) {
            throw std::range_error("Alpha must be between 0 and 1.");
        }
        this->alpha = alpha;

        // set penalty terms
        arma::vec l1_penalty { arma::ones(x.n_cols) };
        if (l1_penalty_factor.n_elem == x.n_cols) {
            // re-scale so that sum(factor) = number of predictors
            l1_penalty = l1_penalty_factor * x.n_cols /
                arma::sum(l1_penalty_factor);
        }
        this->l1_penalty_factor = l1_penalty;

        // the maximum (large enough) lambda that results in all-zero estimates
        arma::vec beta { arma::zeros(x.n_cols) };
        arma::vec grad_beta { beta }, strong_rhs { beta };
        this->l1_lambda_max =
            arma::max(arma::abs(this->gradient(beta)) / l1_penalty) /
            this->x.n_rows / std::max(alpha, 1e-10);

        // take unique lambda and sort descendingly
        lambda = arma::reverse(arma::unique(lambda));
        // construct lambda sequence
        arma::vec lambda_seq {
            arma::zeros(std::max(nlambda, lambda.n_elem))
        };
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

        // initialize the estimate matrix
        this->coef_mat = arma::zeros(x.n_cols, lambda_seq.n_elem);
        this->en_coef_mat = this->coef_mat;
        this->negLogL_vec = arma::zeros(lambda_seq.n_elem);
        this->coef_df_vec = arma::zeros<arma::uvec>(lambda_seq.n_elem);

        // start values
        this->coef0 = beta;
        this->coef = beta;

        // for active set
        arma::uvec is_active_strong { arma::zeros<arma::uvec>(x.n_cols) };

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
            for (size_t j {0}; j < this->x.n_cols; ++j) {
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
                reg_active_fit(beta, is_active_strong, lambda_seq(k) * alpha,
                               lambda_seq(k) * (1 - alpha) / 2, l1_penalty,
                               varying_active_set, max_iter, rel_tol,
                               early_stop, verbose);
                // check kkt condition
                for (size_t j {0}; j < x.n_cols; ++j) {
                    if (is_active_strong(j)) {
                        continue;
                    }
                    if (std::abs(this->gradient(beta, j)) / this->x.n_rows >
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
            this->coef0 = (1 + (1 - alpha) * lambda_seq(k) / 2) * beta;
            this->rescale_coef();
            this->en_coef_mat.col(k) = this->coef;

            // compute naive elastic net estimates
            this->coef0 = beta;
            this->rescale_coef();
            this->coef_mat.col(k) = this->coef;

            // compute negative log-likelihood
            this->negLogL_vec(k) = this->objective();
            this->coef_df_vec(k) = get_coef_df(beta);
        }
        // end of regularized fit
    }

}


#endif
