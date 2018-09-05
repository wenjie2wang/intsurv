//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
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

#ifndef RCPP_COXPH_H
#define RCPP_COXPH_H

#include <RcppArmadillo.h>
#include "utils.hpp"

namespace Intsurv {

    // define class for inputs and outputs
    class RcppCoxph {
    private:
        arma::vec time;           // observed times
        arma::vec event;          // event indicators
        arma::mat z;              // design matrix
        bool hasTies {false};     // if there exists ties on event times
        arma::uvec uni_event_ind; // the first unique index of event times
        arma::uvec event_ind;     // indices of event times
        arma::vec d_time;         // distinct event times
        arma::mat d_z;            // design matrix aggregated at d_time
        arma::vec delta_n;        // event counts at d_time

    public:
        double partial_logL {0}; // partial log-likelihood
        arma::vec coef;          // covariate coefficient estimates
        arma::mat h0;            // baseline hazard estimates

        // constructors
        RcppCoxph(const arma::vec time_,
                  const arma::vec event_,
                  const arma::mat z_)
        {
            // sort event and z based on time
            arma::uvec s_time_ind {arma::sort_index(time_)};
            time = time_.elem(s_time_ind);
            event = event_.elem(s_time_ind);
            z = z_.rows(s_time_ind);

            event_ind = arma::find(event > 0);
            // check if there exists ties on time
            hasTies = any_duplicated(time.elem(event_ind));
            if (hasTies) {
                arma::vec d_time0 {time.elem(event_ind)};
                arma::uvec uni_time_ind {find_unique(time)};
                uni_event_ind = vec_intersection(uni_time_ind, event_ind);
                d_time = time.elem(uni_event_ind);
                // aggregate at distinct event times
                delta_n = event.elem(event_ind);
                delta_n = aggregate_sum(delta_n, d_time0);
                d_z = z.rows(event_ind);
                d_z = aggregate_sum_cols(d_z, d_time0);
            } else {
                d_time = time.elem(event_ind);
                delta_n = arma::ones(event_ind.n_rows);
                d_z = z.rows(event_ind);
            }
        }

        // function members
        double objective(const arma::vec& beta) const;
        arma::vec gradient(const arma::vec& beta) const;

    };

    // the negative loglikelihood function based on the broslow's formula
    double RcppCoxph::objective(const arma::vec& beta) const
    {
        const arma::vec dz_beta {d_z * beta};
        const arma::vec exp_z_beta {arma::exp(z * beta)};
        const arma::vec h0_denom {cum_sum(exp_z_beta, true)};
        if (hasTies) {
            arma::vec log_h0_denom_event {
                arma::log(h0_denom.elem(uni_event_ind))
            };
            return - arma::sum(dz_beta - delta_n % log_h0_denom_event);
        } else {
            arma::vec log_h0_denom_event {
                arma::log(h0_denom.elem(event_ind))
            };
            return - (arma::sum(dz_beta) - arma::sum(log_h0_denom_event));
        }
    }

    // the negative gradient of negative loglikelihood function
    arma::vec RcppCoxph::gradient(const arma::vec& beta) const
    {
        const arma::vec exp_z_beta {arma::exp(z * beta)};
        const arma::vec h0_denom {cum_sum(exp_z_beta, true)};
        arma::mat numer_mat {arma::zeros(z.n_rows, z.n_cols)};
        for (size_t i {0}; i < z.n_rows; ++i) {
            numer_mat.row(i) = exp_z_beta(i) * z.row(i);
        }
        numer_mat = cum_sum_cols(numer_mat, true);
        if (hasTies) {
            arma::vec h0_denom_event {h0_denom.elem(uni_event_ind)};
            numer_mat = numer_mat.rows(uni_event_ind);
            for (size_t j {0}; j < z.n_cols; ++j) {
                numer_mat.col(j) = numer_mat.col(j) % delta_n / h0_denom_event;
            }
            return - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
        } else {
            arma::vec h0_denom_event {h0_denom.elem(event_ind)};
            numer_mat = numer_mat.rows(event_ind);
            for (size_t j {0}; j < z.n_cols; ++j) {
                numer_mat.col(j) = numer_mat.col(j) / h0_denom_event;
            }
            return - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
        }
    }


}


#endif
