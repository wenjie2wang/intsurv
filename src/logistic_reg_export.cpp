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
#include <intsurv.h>

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
    const unsigned int max_iter = 300,
    const double rel_tol = 1e-6,
    const double pmin = 1e-5,
    const bool early_stop = false,
    const bool verbose = false
    )
{
    Intsurv::LogisticReg object {
        x, y, intercept, standardize
    };
    object.set_offset(offset);
    object.fit(start, max_iter, rel_tol, pmin, early_stop, verbose);
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("fitted") = Intsurv::arma2rvec(object.prob_vec_),
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = object.neg_ll_,
            Rcpp::Named("coef_df") = object.coef_df_
            )
        );
}


// fit Firth logistic regression
// [[Rcpp::export]]
Rcpp::List rcpp_firth_logistic(
    const arma::mat& x,
    const arma::vec& y,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec& offset = 0,
    const arma::vec& start = 0,
    const unsigned int max_iter = 300,
    const double rel_tol = 1e-6
    )
{
    Intsurv::LogisticReg object {
        x, y, intercept, standardize
    };
    object.set_offset(offset);
    object.firth_fit(start, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("fitted") = Intsurv::arma2rvec(object.prob_vec_),
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = object.neg_ll_,
            Rcpp::Named("coef_df") = object.coef_df_
            )
        );
}


// regularized logistic model by coordinate-majorization-descent algorithm
// for perticular lambda's
// [[Rcpp::export]]
Rcpp::List rcpp_reg_logistic1(
    const arma::mat& x,
    const arma::vec& y,
    const double l1_lambda = 0,
    const double l2_lambda = 0,
    const arma::vec& l1_penalty_factor = 0,
    const arma::vec& start = 0,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec& offset = 0,
    const unsigned int max_iter = 300,
    const double rel_tol = 1e-5,
    const double pmin = 1e-5,
    const bool early_stop = false,
    const bool verbose = false
    )
{
    Intsurv::LogisticReg object { x, y, intercept, standardize };
    object.set_offset(offset);
    object.regularized_fit(l1_lambda, l2_lambda, l1_penalty_factor,
                           start, max_iter, rel_tol, pmin,
                           early_stop, verbose);
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef_),
        // Rcpp::Named("en_coef") = Intsurv::arma2rvec(object.en_coef_),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = object.neg_ll_,
            Rcpp::Named("coef_df") = object.coef_df_
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("l1_lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("l1_lambda") = object.l1_lambda_,
            Rcpp::Named("l2_lambda") = object.l2_lambda_,
            Rcpp::Named("l1_penalty_factor") =
            Intsurv::arma2rvec(object.l1_penalty_factor_)
            )
        );
}


// for a sequence of lambda's and a given alpha
// [[Rcpp::export]]
Rcpp::List rcpp_reg_logistic2(
    const arma::mat& x,
    const arma::vec& y,
    const arma::vec& lambda = 0,
    const double alpha = 1,
    const unsigned int nlambda = 1,
    const double lambda_min_ratio = 1e-4,
    const arma::vec& penalty_factor = 0,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec& offset = 0,
    const unsigned int max_iter = 300,
    const double rel_tol = 1e-5,
    const double pmin = 1e-5,
    const bool early_stop = false,
    const bool verbose = false
    )
{
    Intsurv::LogisticReg object { x, y, intercept, standardize };
    object.set_offset(offset);
    object.regularized_fit(lambda, alpha, nlambda, lambda_min_ratio,
                           penalty_factor, max_iter, rel_tol,
                           pmin, early_stop, verbose);
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef_mat_,
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = object.n_obs_,
            Rcpp::Named("negLogL") = Intsurv::arma2rvec(object.neg_ll_vec_),
            Rcpp::Named("coef_df") = Intsurv::arma2rvec(object.coef_df_vec_)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("lambda_max") = object.l1_lambda_max_,
            Rcpp::Named("lambda") = Intsurv::arma2rvec(object.lambda_vec_),
            Rcpp::Named("alpha") = object.alpha_,
            Rcpp::Named("l1_penalty_factor") =
            Intsurv::arma2rvec(object.l1_penalty_factor_)
            )
        );
}
