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

#ifndef INTSURV_NONPARAMETRIC_H
#define INTSURV_NONPARAMETRIC_H

#include <RcppArmadillo.h>
#include <map>
#include "utils.h"


namespace Intsurv {

    // define class for inputs and outputs
    class NelsonAalen {
    protected:
        arma::uvec des_event_ind_; // index sorting events descendingly
        arma::uvec asc_time_ind_;  // index sorting times ascendingly
        arma::vec time_;           // (sorted) observed times
        arma::vec event_;          // (sorted) event indicators

    public:
        arma::vec uni_event_time_; // sorted unique event_ times
        // at each unique event time
        arma::vec delta_event_;    // number of events
        arma::vec riskset_size_;   // size of riskset
        arma::vec inst_rate_;      // instantaneous hazard rate function
        arma::vec cum_rate_;       // cumulative hazard rate function

        // constructors
        NelsonAalen(const arma::vec& time,
                    const arma::vec& event) :
            time_ (time),
            event_ (event)
        {
            // sort based on time and event
            // time: ascending order
            // event: events first, then censoring at the same time point
            des_event_ind_ = arma::sort_index(event_, "descend");
            time_ = time_.elem(des_event_ind_);
            event_ = event_.elem(des_event_ind_);
            asc_time_ind_ = arma::stable_sort_index(time_, "ascend");
            time_ = time_.elem(asc_time_ind_);
            event_ = event_.elem(asc_time_ind_);

            // check for any tied events
            arma::uvec event_ind { arma::find(event_ > 0) };
            arma::vec d_time0 { time_.elem(event_ind) };
            bool has_ties { any_duplicated(d_time0) };

            // default for no ties
            delta_event_ = event_.elem(event_ind);
            arma::uvec uni_event_ind { event_ind };
            uni_event_time_ = d_time0;

            if (has_ties) {
                arma::uvec uni_time_ind { find_first_unique(time_) };
                uni_event_ind = vec_intersection(uni_time_ind, event_ind);
                uni_event_time_ = time_.elem(uni_event_ind);
                delta_event_ = aggregate_sum(delta_event_, d_time0);
            }
            riskset_size_ = arma::ones(time_.n_elem);
            riskset_size_ = cum_sum(riskset_size_, true).elem(uni_event_ind);
            inst_rate_ = delta_event_ / riskset_size_;
            cum_rate_ = cum_sum(inst_rate_);
        }

        // function members
        // estimate instantaneous hazard rate by nearest left neighbor
        // the nearest right neighbor is used if the left neighbor doesn't exist
        inline arma::vec step_inst_rate(const arma::vec& new_time) const
        {
            // create a map for fast comparison
            std::map<double, double> step_map;
            for (size_t i {0}; i < uni_event_time_.n_elem; ++i) {
                step_map.insert(std::make_pair(uni_event_time_(i),
                                               inst_rate_(i)));
            }
            arma::vec res { arma::zeros(new_time.n_elem) };
            std::map<double, double>::iterator it;
            for (size_t i {0}; i < new_time.n_elem; ++i) {
                it = step_map.upper_bound(new_time(i));
                if (it == step_map.begin()) {
                    res(i) = inst_rate_(0);
                } else {
                    --it;
                    res(i) = it->second;
                }
            }
            return res;
        }

        // estimate cumulative hazard function at new time_ points
        inline arma::vec step_cum_rate(const arma::vec& new_time) const
        {
            // create a map for fast comparison
            std::map<double, double> step_map;
            for (size_t i {0}; i < uni_event_time_.n_elem; ++i) {
                step_map.insert(std::make_pair(uni_event_time_(i),
                                               cum_rate_(i)));
            }
            arma::vec res { arma::zeros(new_time.n_elem) };
            std::map<double, double>::iterator it;
            for (size_t i {0}; i < new_time.n_elem; ++i) {
                it = step_map.upper_bound(new_time(i));
                if (it != step_map.begin()) {
                    --it;
                    res(i) = it->second;
                }
            }
            return res;
        }



    };


}


#endif
