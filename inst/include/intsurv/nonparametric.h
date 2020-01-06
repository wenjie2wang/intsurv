//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2020  Wenjie Wang <wjwang.stat@gmail.com>
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

#ifndef NONPARAMETRIC_H
#define NONPARAMETRIC_H

#include <RcppArmadillo.h>
#include <map>
#include "utils.h"


namespace Intsurv {

    // define class for inputs and outputs
    class NelsonAalen {
    private:
        arma::uvec des_event_ind; // index sorting events descendingly
        arma::uvec asc_time_ind;  // index sorting times ascendingly
        arma::vec time;         // (sorted) observed times
        arma::vec event;        // (sorted) event indicators

    public:
        arma::vec uni_event_time; // sorted unique event times
        // at each unique event time
        arma::vec delta_event;    // number of events
        arma::vec riskset_size;   // size of riskset
        arma::vec inst_rate;      // instantaneous hazard rate function
        arma::vec cum_rate;       // cumulative hazard rate function

        // constructors
        NelsonAalen(const arma::vec& time_,
                    const arma::vec& event_)
        {
            // sort based on time and event
            // time: ascending order
            // event: events first, then censoring at the same time point
            des_event_ind = arma::sort_index(event_, "descend");
            time = time_.elem(des_event_ind);
            event = event_.elem(des_event_ind);
            asc_time_ind = arma::stable_sort_index(time, "ascend");
            time = time.elem(asc_time_ind);
            event = event.elem(asc_time_ind);

            // check for any tied events
            arma::uvec event_ind { arma::find(event > 0) };
            arma::vec d_time0 { time.elem(event_ind) };
            bool hasTies { any_duplicated(d_time0) };

            // default for no ties
            delta_event = event.elem(event_ind);
            arma::uvec uni_event_ind { event_ind };
            uni_event_time = d_time0;

            if (hasTies) {
                arma::uvec uni_time_ind { find_first_unique(time) };
                uni_event_ind = vec_intersection(uni_time_ind, event_ind);
                uni_event_time = time.elem(uni_event_ind);
                delta_event = aggregate_sum(delta_event, d_time0);
            }
            riskset_size = arma::ones(time.n_elem);
            riskset_size = cum_sum(riskset_size, true).elem(uni_event_ind);
            inst_rate = delta_event / riskset_size;
            cum_rate = cum_sum(inst_rate);
        }

        // function members
        // estimate instantaneous hazard rate by nearest left neighbor
        // the nearest right neighbor is used if the left neighbor doesn't exist
        inline arma::vec step_inst_rate(const arma::vec& new_time) const
        {
            // create a map for fast comparison
            std::map<double, double> step_map;
            for (size_t i {0}; i < this->uni_event_time.n_elem; ++i) {
                step_map.insert(std::make_pair(this->uni_event_time(i),
                                               this->inst_rate(i)));
            }
            arma::vec res { arma::zeros(new_time.n_elem) };
            std::map<double, double>::iterator it;
            for (size_t i {0}; i < new_time.n_elem; ++i) {
                it = step_map.upper_bound(new_time(i));
                if (it == step_map.begin()) {
                    res(i) = this->inst_rate(0);
                } else {
                    --it;
                    res(i) = it->second;
                }
            }
            return res;
        }

        // estimate cumulative hazard function at new time points
        inline arma::vec step_cum_rate(const arma::vec& new_time) const
        {
            // create a map for fast comparison
            std::map<double, double> step_map;
            for (size_t i {0}; i < this->uni_event_time.n_elem; ++i) {
                step_map.insert(std::make_pair(this->uni_event_time(i),
                                               this->cum_rate(i)));
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
