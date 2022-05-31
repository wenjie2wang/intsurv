//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2022  Wenjie Wang <wang@wwenjie.org>
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
#include <intsurv/Control.h>
#include <intsurv/CoxphReg.h>
#include <intsurv/utils.h>


// fit regular Cox model that allows non-integer "event" and tied events
// [[Rcpp::export]]
Rcpp::List rcpp_coxph(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& x,
    const arma::vec& offset = 0,
    const arma::vec& start = 0,
    const bool standardize = true,
    const unsigned int maxit = 200,
    const double epsilon = 1e-4,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control { maxit, epsilon, standardize, verbose };
    intsurv::CoxphReg object { time, event, x, control };
    object.set_start(start)->set_offset(offset, false);
    object.fit();
    // compute baseline hazard and survival function
    object.compute_haz_surv_time();
    object.compute_censor_haz_surv_time();
    object.est_haz_surv();
    return Rcpp::List::create(
        Rcpp::Named("coef") = intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("offset") = object.control_.offset_
            ),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = intsurv::arma2rvec(object.unique_time_),
            Rcpp::Named("h0") = intsurv::arma2rvec(object.h0_est_),
            Rcpp::Named("H0") = intsurv::arma2rvec(object.H0_est_),
            Rcpp::Named("S0") = intsurv::arma2rvec(object.S0_est_),
            Rcpp::Named("hc") = intsurv::arma2rvec(object.hc_est_),
            Rcpp::Named("Hc") = intsurv::arma2rvec(object.Hc_est_),
            Rcpp::Named("Sc") = intsurv::arma2rvec(object.Sc_est_)
            ),
        Rcpp::Named("control") = Rcpp::List::create(
            Rcpp::Named("maxit") = object.control_.max_iter_,
            Rcpp::Named("epsilon") = object.control_.epsilon_,
            Rcpp::Named("standardize") = object.control_.standardize_,
            Rcpp::Named("start") = intsurv::arma2rvec(object.control_.start_),
            Rcpp::Named("verbose") = object.control_.verbose_
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
    const bool standardize = true,
    const unsigned int maxit = 200,
    const double epsilon = 1e-4,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control { maxit, epsilon, standardize, verbose };
    control.net(penalty_factor, varying_active)->
        net_fit(l1_lambda, l2_lambda);
    intsurv::CoxphReg object { time, event, x, control };
    object.set_start(start)->set_offset(offset, false);
    object.net_fit();
    object.compute_haz_surv_time();
    object.compute_censor_haz_surv_time();
    object.est_haz_surv();
    return Rcpp::List::create(
        Rcpp::Named("coef") = intsurv::arma2rvec(object.coef_),
        // Rcpp::Named("en_coef") = intsurv::arma2rvec(object.en_coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("offset") = object.control_.offset_
            ),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = intsurv::arma2rvec(object.unique_time_),
            Rcpp::Named("h0") = intsurv::arma2rvec(object.h0_est_),
            Rcpp::Named("H0") = intsurv::arma2rvec(object.H0_est_),
            Rcpp::Named("S0") = intsurv::arma2rvec(object.S0_est_),
            Rcpp::Named("hc") = intsurv::arma2rvec(object.hc_est_),
            Rcpp::Named("Hc") = intsurv::arma2rvec(object.Hc_est_),
            Rcpp::Named("Sc") = intsurv::arma2rvec(object.Sc_est_)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda") = object.control_.l1_lambda_,
            Rcpp::Named("l2_lambda") = object.control_.l2_lambda_,
            Rcpp::Named("penalty_factor") =
            intsurv::arma2rvec(object.control_.penalty_factor_)
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
    const bool standardize = true,
    const unsigned int maxit = 200,
    const double epsilon = 1e-4,
    const bool verbose = false
    )
{
    intsurv::Control control { maxit, epsilon, standardize, verbose };
    control.net(penalty_factor, varying_active)->
        net_path(nlambda, lambda_min_ratio, alpha, lambda);
    intsurv::CoxphReg object { time, event, x, control };
    object.set_offset(offset);
    object.net_path();
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef_mat_,
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("offset") = object.control_.offset_
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("lambda_max") = object.lambda_max_,
            Rcpp::Named("alpha") = object.control_.alpha_,
            Rcpp::Named("lambda") = intsurv::arma2rvec(object.control_.lambda_),
            Rcpp::Named("penalty_factor") =
            intsurv::arma2rvec(object.control_.penalty_factor_)
            )
        );
}
