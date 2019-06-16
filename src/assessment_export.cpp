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


#include <RcppArmadillo.h>
#include <stdexcept>
#include "utils.h"
#include "assessment.h"


// weighted C-index
// [[Rcpp::export]]
Rcpp::List rcpp_cIndex(arma::vec time,
                       arma::vec event,
                       arma::vec risk_score,
                       arma::vec weight)
{
    unsigned int nObs { time.n_elem };
    if (nObs <= 1) {
        throw std::logic_error(
            "The inputs must have length greater than one."
            );
    }
    if (weight.n_elem == 1) {
        weight = arma::ones(nObs);
    }
    if (event.n_elem != nObs || risk_score.n_elem != nObs ||
        weight.n_elem != nObs) {
        throw std::logic_error(
            "The inputs must have the same length."
            );
    }
    // create a concordance object
    Intsurv::Concordance c_obj {
        Intsurv::Concordance(time, event, risk_score, weight)
    };
    return Rcpp::List::create(
        Rcpp::Named("index") = c_obj.index,
        Rcpp::Named("concordant") = c_obj.concordant,
        Rcpp::Named("comparable") = c_obj.comparable,
        Rcpp::Named("tied_risk") = c_obj.tied_risk
        );
}
