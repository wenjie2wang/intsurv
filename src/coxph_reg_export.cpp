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

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <intsurv/coxph_reg.h>
#include <intsurv/utils.h>

// fit regular Cox model that allows non-integer "event" and tied events
// [[Rcpp::export]]
Rcpp::List rcpp_coxph(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& x,
    const arma::vec& offset,
    const arma::vec& start,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const unsigned int verbose = 0
    )
{
    // define object
    Intsurv::CoxphReg object { time, event, x };
    // set offset if it is not zero
    if (! offset.is_empty()) {
        object.set_offset(offset, false);
    }
    // model-fitting
    object.fit(start, max_iter, epsilon, verbose);
    // compute baseline hazard and survival function
    object.compute_haz_surv_time();
    object.compute_censor_haz_surv_time();
    object.est_haz_surv();
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.sample_size()
            ),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(object.unique_time_),
            Rcpp::Named("h0") = Intsurv::arma2rvec(object.h0_est_),
            Rcpp::Named("H0") = Intsurv::arma2rvec(object.H0_est_),
            Rcpp::Named("S0") = Intsurv::arma2rvec(object.S0_est_),
            Rcpp::Named("hc") = Intsurv::arma2rvec(object.hc_est_),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(object.Hc_est_),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(object.Sc_est_)
            )
        );
}


// fit regularized Cox model with coordinate-majorizatio-descent algorithm
// that allows non-integer "event" and tied events
// for particular lambda's
// [[Rcpp::export]]
Rcpp::List rcpp_coxnet1(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& x,
    const double l1_lambda,
    const double l2_lambda,
    const arma::vec& penalty_factor,
    const arma::vec& offset,
    const arma::vec& start,
    const bool varying_active = true,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const unsigned int verbose = 0
    )
{
    Intsurv::CoxphReg object { time, event, x };
    // set offset if it is not zero
    if (! Intsurv::isAlmostEqual(arma::sum(arma::abs(offset)), 0.0)) {
        object.set_offset(offset, false);
    }
    object.net_fit(l1_lambda, l2_lambda, penalty_factor,
                   start, varying_active, max_iter, epsilon, verbose);
    object.compute_haz_surv_time();
    object.compute_censor_haz_surv_time();
    object.est_haz_surv();
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        // Rcpp::Named("en_coef") = Intsurv::arma2rvec(object.en_coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.sample_size()
            ),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(object.unique_time_),
            Rcpp::Named("h0") = Intsurv::arma2rvec(object.h0_est_),
            Rcpp::Named("H0") = Intsurv::arma2rvec(object.H0_est_),
            Rcpp::Named("S0") = Intsurv::arma2rvec(object.S0_est_),
            Rcpp::Named("hc") = Intsurv::arma2rvec(object.hc_est_),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(object.Hc_est_),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(object.Sc_est_)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda") = object.l1_lambda_,
            Rcpp::Named("l2_lambda") = object.l2_lambda_,
            Rcpp::Named("penalty_factor") =
            Intsurv::arma2rvec(object.penalty_factor_)
            )
        );
}


// for a sequence of lambda's and a given alpha
// [[Rcpp::export]]
Rcpp::List rcpp_coxnet2(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& x,
    const arma::vec& lambda,
    const double alpha,
    const unsigned int nlambda,
    const double lambda_min_ratio,
    const arma::vec& penalty_factor,
    const arma::vec& offset,
    const bool varying_active = true,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const bool verbose = false
    )
{
    Intsurv::CoxphReg object { time, event, x };
    // set offset if it is not zero
    if (! Intsurv::isAlmostEqual(arma::sum(arma::abs(offset)), 0.0)) {
        object.set_offset(offset, false);
    }
    object.net_path(lambda, alpha, nlambda, lambda_min_ratio,
                    penalty_factor, varying_active,
                    max_iter, epsilon, verbose);
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef_mat_,
        // Rcpp::Named("en_coef") = object.en_coef_mat_,
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.sample_size()
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("lambda_max") = object.lambda_max_,
            Rcpp::Named("alpha") = object.alpha_,
            Rcpp::Named("lambda") = Intsurv::arma2rvec(object.lambda_vec_),
            Rcpp::Named("penalty_factor") =
            Intsurv::arma2rvec(object.penalty_factor_)
            )
        );
}
