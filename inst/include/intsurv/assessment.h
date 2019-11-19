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

#ifndef ASSESSMENT_H
#define ASSESSMENT_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {
// A straightforward implementation of Harrel's C-index that allows weights
    class Concordance {
    public:
        double index {0.0};      // C-index result
        double comparable {0.0}; // number of comparable pairs
        double concordant {0.0}; // number of concordant pairs
        // number of comparable pairs with tied risk scores
        double tied_risk {0.0};

        // constructors
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
                // only comparable when event(i) > 0
                if (event(i) > 0) {
                    for (size_t j { i + 1 }; j < nObs; ++j) {
                        // not comparable if time(i) = time(j) and event(j) = 1
                        if (isAlmostEqual(time(i), time(j)) && event(j) > 0) {
                            continue;
                        }
                        // otherwise, comparable
                        this->comparable += weight(i);
                        // determine the concordance
                        if (isAlmostEqual(time(i), time(j))) {
                            // case 1. tied times but event(j) = 0
                            if (is_gt(risk_score(i), risk_score(j))) {
                                this->concordant += weight(i);
                            }
                        } else {
                            // case 2. distinct times
                            if (isAlmostEqual(risk_score(i), risk_score(j))) {
                                this->tied_risk += weight(i);
                            } else if (risk_score(i) > risk_score(j)) {
                                this->concordant += weight(i);
                            }
                        }
                    }
                }
            }
            this->index = (this->concordant + this->tied_risk / 2) /
                this->comparable;
        }
        // constructor with weights
        Concordance(arma::vec time,
                    arma::vec event,
                    arma::vec risk_score)
        {
            Concordance(time, event, risk_score, arma::ones(time.n_elem));
        }

    };
}

#endif
