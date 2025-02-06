//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2025  Wenjie Wang <wang@wwenjie.org>
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
#include <intsurv/LogisticReg.h>
#include <intsurv/utils.h>

// fitting regular logistic model by monotonic quadratic approximation algorithm
// non-integer y vector is allowed
// [[Rcpp::export]]
Rcpp::List rcpp_logistic(
    const arma::mat& x,
    const arma::vec& y,
    const bool intercept = true,
    const arma::vec& offset = 0,
    const arma::vec& start = 0,
    const bool standardize = true,
    const unsigned int maxit = 200,
    const double epsilon = 1e-4,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control { maxit, epsilon, standardize, verbose };
    control.logistic(intercept, pmin);
    intsurv::LogisticReg object { x, y, control };
    object.set_start(start)->set_offset(offset);
    object.fit();
    double neg_ll { object.objective() };
    arma::vec prob { object.predict() };
    return Rcpp::List::create(
        Rcpp::Named("coef") = intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("prob") = intsurv::arma2rvec(prob),
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("intercept") = object.control_.intercept_,
            Rcpp::Named("offset") = intsurv::arma2rvec(object.control_.offset_),
            Rcpp::Named("negLogL") = neg_ll
            ),
        Rcpp::Named("control") = Rcpp::List::create(
            Rcpp::Named("maxit") = object.control_.max_iter_,
            Rcpp::Named("epsilon") = object.control_.epsilon_,
            Rcpp::Named("standardize") = object.control_.standardize_,
            Rcpp::Named("start") = intsurv::arma2rvec(object.control_.start_),
            Rcpp::Named("verbose") = object.control_.verbose_,
            Rcpp::Named("pmin") = object.control_.pmin_
            )
        );
}


// regularized logistic model by coordinate-majorization-descent algorithm
// for perticular lambda's
// [[Rcpp::export]]
Rcpp::List rcpp_lognet1(
    const arma::mat& x,
    const arma::vec& y,
    const double l1_lambda,
    const double l2_lambda,
    const arma::vec& penalty_factor,
    const arma::vec& start,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec& offset = 0,
    const bool varying_active = true,
    const unsigned int maxit = 200,
    const double epsilon = 1e-4,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control { maxit, epsilon, standardize, verbose };
    control.logistic(intercept, pmin)->
        net(penalty_factor, varying_active)->
        net_fit(l1_lambda, l2_lambda);
    intsurv::LogisticReg object { x, y, control };
    object.set_start(start)->set_offset(offset);
    object.net_fit();
    double neg_ll { object.objective() };
    unsigned int coef_df { intsurv::compute_coef_df(object.coef_) };
    return Rcpp::List::create(
        Rcpp::Named("coef") = intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = neg_ll,
            Rcpp::Named("coef_df") = coef_df
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("l1_lambda") = object.control_.l1_lambda_,
            Rcpp::Named("l2_lambda") = object.control_.l2_lambda_,
            Rcpp::Named("penalty_factor") =
            intsurv::arma2rvec(object.control_.penalty_factor_)
            )
        );
}


// for a sequence of lambda's and a given alpha
// [[Rcpp::export]]
Rcpp::List rcpp_lognet2(
    const arma::mat& x,
    const arma::vec& y,
    const arma::vec& lambda,
    const unsigned int nlambda,
    const double lambda_min_ratio,
    const double alpha,
    const arma::vec& penalty_factor,
    const bool intercept = true,
    const arma::vec& offset = 0,
    const bool standardize = true,
    const bool varying_active = true,
    const unsigned int maxit = 200,
    const double epsilon = 1e-4,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control { maxit, epsilon, standardize, verbose };
    control.logistic(intercept, pmin)->
        net(penalty_factor, varying_active)->
        net_path(nlambda, lambda_min_ratio, alpha, lambda);
    intsurv::LogisticReg object { x, y, control };
    object.set_offset(offset);
    object.net_path();
    return Rcpp::List::create(
        Rcpp::Named("coef") = intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("neg_ll") = object.neg_ll_,
            Rcpp::Named("df") = object.coef_df_,
            Rcpp::Named("bic") = object.bic_
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("lambda_max") = object.lambda_max_,
            Rcpp::Named("lambda") = intsurv::arma2rvec(object.control_.lambda_),
            Rcpp::Named("alpha") = object.control_.alpha_,
            Rcpp::Named("penalty_factor") =
            intsurv::arma2rvec(object.control_.penalty_factor_)
            ),
        Rcpp::Named("path") = Rcpp::List::create(
            Rcpp::Named("coef_mat") = object.coef_mat_,
            Rcpp::Named("neg_ll") = intsurv::arma2rvec(object.neg_ll_path_),
            Rcpp::Named("df") = intsurv::arma2rvec(object.coef_df_path_),
            Rcpp::Named("bic") = intsurv::arma2rvec(object.bic_path_)
            )
        );
}
