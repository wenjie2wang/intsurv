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
#include <intsurv/logistic_reg.h>
#include <intsurv/utils.h>

// fitting regular logistic model by monotonic quadratic approximation algorithm
// non-integer y vector is allowed
// [[Rcpp::export]]
Rcpp::List rcpp_logistic(
    const arma::mat& x,
    const arma::vec& y,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec& offset = 0,
    const arma::vec& start = 0,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    Intsurv::LogisticReg object {
        x, y, intercept, standardize
    };
    object.set_offset(offset);
    object.set_pmin(pmin);
    object.fit(start, max_iter, epsilon, verbose);
    double neg_ll { object.objective() };
    arma::vec prob { object.predict() };
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("prob") = Intsurv::arma2rvec(prob),
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = neg_ll
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
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    Intsurv::LogisticReg object { x, y, intercept, standardize };
    object.set_offset(offset);
    object.set_pmin(pmin);
    object.net_fit(l1_lambda, l2_lambda, penalty_factor, start,
                   varying_active, max_iter, epsilon, verbose);
    double neg_ll { object.objective() };
    unsigned int coef_df { Intsurv::compute_coef_df(object.coef_) };
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = neg_ll,
            Rcpp::Named("coef_df") = coef_df
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("l1_lambda") = object.l1_lambda_,
            Rcpp::Named("l2_lambda") = object.l2_lambda_,
            Rcpp::Named("penalty_factor") =
            Intsurv::arma2rvec(object.penalty_factor_)
            )
        );
}


// for a sequence of lambda's and a given alpha
// [[Rcpp::export]]
Rcpp::List rcpp_lognet2(
    const arma::mat& x,
    const arma::vec& y,
    const arma::vec& lambda,
    const double alpha,
    const unsigned int nlambda,
    const double lambda_min_ratio,
    const arma::vec& penalty_factor,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec& offset = 0,
    const bool varying_active = true,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    Intsurv::LogisticReg object { x, y, intercept, standardize };
    object.set_offset(offset);
    object.set_pmin(pmin);
    object.net_path(lambda, alpha, nlambda, lambda_min_ratio,
                    penalty_factor, varying_active, max_iter, epsilon,
                    verbose);
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef_mat_,
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("lambda_max") = object.lambda_max_,
            Rcpp::Named("lambda") = Intsurv::arma2rvec(object.lambda_vec_),
            Rcpp::Named("alpha") = object.alpha_,
            Rcpp::Named("penalty_factor") =
            Intsurv::arma2rvec(object.penalty_factor_)
            )
        );
}
