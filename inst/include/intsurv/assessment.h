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

#ifndef INTSURV_ASSESSMENT_H
#define INTSURV_ASSESSMENT_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {
// A straightforward implementation of Harrel's C-index_ that allows weights
    class Concordance {
    public:
        double index_ {0.0};      // C-index_ result
        double comparable_ {0.0}; // number of comparable_ pairs
        double concordant_ {0.0}; // number of concordant_ pairs
        // number of comparable_ pairs with tied risk scores
        double tied_risk_ {0.0};

        // constructors with weight
        Concordance(arma::vec time,
                    arma::vec event,
                    arma::vec risk_score,
                    arma::vec weight)
        {
            // sort based on time and event
            // time: ascending order
            // event: events first, then censoring at the same time point
            arma::uvec des_event_ind { arma::sort_index(event, "descend") };
            arma::uvec asc_time_ind {
                arma::stable_sort_index(time.elem(des_event_ind), "ascend")
            };
            arma::uvec ord { des_event_ind.elem(asc_time_ind) };

            // do the actual sorting
            time = time.elem(ord);
            event = event.elem(ord);
            risk_score = risk_score.elem(ord);
            weight = weight.elem(ord);

            // do the actual computation
            unsigned int nObs { time.n_elem };
            for (size_t i { 0 }; i < nObs - 1; ++i) {
                // only comparable_ when event(i) > 0
                if (event(i) > 0) {
                    for (size_t j { i + 1 }; j < nObs; ++j) {
                        // not comparable_ if time(i) = time(j) and event(j) = 1
                        if (isAlmostEqual(time(i), time(j)) && event(j) > 0) {
                            continue;
                        }
                        // otherwise, comparable_
                        comparable_ += weight(j);
                        // determine the concordance
                        if (isAlmostEqual(time(i), time(j))) {
                            // case 1. tied times but event(j) = 0
                            if (is_gt(risk_score(i), risk_score(j))) {
                                concordant_ += weight(j);
                            }
                        } else {
                            // case 2. distinct times
                            if (isAlmostEqual(risk_score(i), risk_score(j))) {
                                tied_risk_ += weight(j);
                            } else if (risk_score(i) > risk_score(j)) {
                                concordant_ += weight(j);
                            }
                        }
                    }
                }
            }
            index_ = (concordant_ + tied_risk_ / 2) /
                comparable_;
        }

        // constructor without weight
        Concordance(arma::vec time,
                    arma::vec event,
                    arma::vec risk_score)
        {
            Concordance(time, event, risk_score, arma::ones(time.n_elem));
        }

    };
}

#endif
