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
    class RcppCoxph {
    private:
        arma::vec time;           // sorted observed times
        arma::vec event;          // sorted event indicators
        arma::mat x;              // sorted design matrix
        bool hasTies {false};     // if there exists ties on event times
        arma::uvec uni_event_ind; // the index indicating the first record
                                  // on each distinct event time

        arma::uvec event_ind;     // indices of event times
        arma::vec d_time;         // distinct event times
        arma::mat d_x;            // design matrix aggregated at d_time
        arma::vec delta_n;        // event counts at d_time
        arma::vec riskset_size;   // size of risk-set at d_time

    public:
        double partial_logL {0}; // partial log-likelihood
        arma::vec coef;          // covariate coefficient estimates
        arma::vec h0;            // baseline hazard estimates

        // constructors
        RcppCoxph(const arma::vec time_,
                  const arma::vec event_,
                  const arma::mat x_)
        {
            // sort based on time and event
            // time: ascending order
            // event: events first, then censoring at the same time point
            arma::uvec s_event_ind {arma::sort_index(event_, "descend")};
            time = time_.elem(s_event_ind);
            event = event_.elem(s_event_ind);
            x = x_.rows(s_event_ind);
            arma::uvec s_time_ind {arma::stable_sort_index(time, "ascend")};
            time = time.elem(s_time_ind);
            event = event.elem(s_time_ind);
            x = x.rows(s_time_ind);
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
            if (hasTies) {
                arma::uvec uni_time_ind {find_first_unique(time)};
                // re-define uni_event_ind
                uni_event_ind = vec_intersection(uni_time_ind, event_ind);
                d_time = time.elem(uni_event_ind);
                // aggregate at distinct event times
                delta_n = aggregate_sum(delta_n, d_time0);
                d_x = aggregate_sum(d_x, d_time0);
            }
            riskset_size = arma::ones(time.n_elem);
            riskset_size = cum_sum(riskset_size, true).elem(uni_event_ind);
        }

        // function that computes objective function only
        inline double objective(const arma::vec& beta) const;

        // function that computes gradients only
        inline arma::vec gradient(const arma::vec& beta) const;
        inline double gradient(const arma::vec& beta,
                               const unsigned int k) const;

        // function that computes objective and overwrite gradients
        inline double objective(const arma::vec& beta, arma::vec& grad) const;

        // function computing B&L covariance lower bound matrix
        inline arma::mat bl_cov_lowerbound_mat() const;

        // function computing B&L lower bound for step size
        inline double bl_step_lowerbound(const arma::mat& x,
                                         const arma::vec& h) const;

        // get the lower bound of second derivative in CMD algorithm
        inline arma::vec cmd_lowerbound() const;

        // some simple functions
        unsigned int sample_size() const
        {
            return time.n_elem;
        }

        // helper function to access some private members
        arma::vec get_time() const{ return time; }
        arma::vec get_event() const { return event; }
        arma::vec get_x() const { return x; }

    };


    // the negative log-likelihood function based on the broslow's formula
    inline double RcppCoxph::objective(const arma::vec& beta) const
    {
        const arma::vec dx_beta {d_x * beta};
        const arma::vec exp_x_beta {arma::exp(x * beta)};
        const arma::vec h0_denom {cum_sum(exp_x_beta, true)};
        arma::vec log_h0_denom_event {arma::log(h0_denom.elem(uni_event_ind))};
        return - arma::sum(dx_beta - delta_n % log_h0_denom_event);
    }

    // the gradient of negative loglikelihood function
    inline arma::vec RcppCoxph::gradient(const arma::vec& beta) const
    {
        const arma::vec exp_x_beta {arma::exp(x * beta)};
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
    inline double RcppCoxph::gradient(const arma::vec& beta,
                                      const unsigned int k) const
    {
        const arma::vec exp_x_beta { arma::exp(x * beta) };
        arma::vec h0_denom { cum_sum(exp_x_beta, true) };
        arma::vec numer { cum_sum(mat2vec(x.col(k) % exp_x_beta), true) };
        h0_denom = h0_denom.elem(uni_event_ind);
        numer = numer.elem(uni_event_ind);
        return - arma::sum(d_x.col(k) - delta_n % numer / h0_denom);
    }

    // the negative log-likelihood function based on the broslow's formula
    inline double RcppCoxph::objective(const arma::vec& beta,
                                       arma::vec& grad) const
    {
        const arma::vec dx_beta {d_x * beta};
        const arma::vec exp_x_beta {arma::exp(x * beta)};
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
    inline arma::mat RcppCoxph::bl_cov_lowerbound_mat() const
    {
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
        return res;
    }
    inline double RcppCoxph::bl_step_lowerbound(const arma::mat& x,
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
    inline arma::vec RcppCoxph::cmd_lowerbound() const
    {
        // compute D_k at each event time point
        arma::mat max_risk_x { Intsurv::cum_max(x, true) };
        arma::mat min_risk_x { Intsurv::cum_min(x, true) };
        max_risk_x = max_risk_x.rows(uni_event_ind);
        min_risk_x = min_risk_x.rows(uni_event_ind);
        arma::vec res { arma::zeros(x.n_cols) };
        for (size_t k {0}; k < x.n_cols; ++k) {
            res(k) = arma::sum(
                pow(max_risk_x.col(k) - min_risk_x.col(k), 2) % delta_n
                ) / (4 * x.n_rows);
        }
        return res;
    }

}


#endif
