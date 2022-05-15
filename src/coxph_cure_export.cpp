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
#include <intsurv/coxph_cure.h>
#include <intsurv/cv_coxph_cure.h>
#include <intsurv/control.h>
#include <intsurv/subset.h>
#include <intsurv/utils.h>

// fit regular Cox cure rate model by EM algorithm
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& surv_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const unsigned int bootstrap = 0,
    const arma::vec& surv_start = 0,
    const arma::vec& cure_start = 0,
    const arma::vec& surv_offset = 0,
    const arma::vec& cure_offset = 0,
    const bool surv_standardize = true,
    const bool cure_standardize = true,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const unsigned int surv_max_iter = 100,
    const double surv_epsilon = 1e-4,
    const unsigned int cure_max_iter = 100,
    const double cure_epsilon = 1e-4,
    const unsigned int tail_completion = 1,
    const double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control0 { max_iter, epsilon };
    control0.cure(tail_completion, tail_tau)->
        set_verbose(verbose);
    intsurv::Control surv_control {
        surv_max_iter, surv_epsilon, surv_standardize,
        intsurv::less_verbose(verbose, 3)
    };
    surv_control.set_start(surv_start)->
        set_offset(surv_offset);
    intsurv::Control cure_control {
        cure_max_iter, cure_epsilon, cure_standardize,
        intsurv::less_verbose(verbose, 3)
    };
    cure_control.logistic(cure_intercept, pmin)->
        set_start(cure_start)->
        set_offset(cure_offset);
    // define object
    intsurv::CoxphCure obj {
        time, event, surv_x, cure_x,
        control0, surv_control, cure_control
    };
    // model-fitting
    obj.fit();
    // initialize bootstrap estimates
    arma::mat boot_surv_coef_mat, boot_cure_coef_mat;
    if (bootstrap > 0) {
        boot_surv_coef_mat = arma::zeros(obj.surv_coef_.n_elem, bootstrap);
        boot_cure_coef_mat = arma::zeros(obj.cure_coef_.n_elem, bootstrap);
        const arma::uvec& case1_ind { obj.case1_ind_ };
        const arma::uvec& case2_ind { obj.case2_ind_ };
        for (size_t i {0}; i < bootstrap; ++i) {
            // generate a bootstrap sample
            arma::uvec boot_ind {
                intsurv::vec_union(
                    intsurv::bootstrap_sample(case1_ind),
                    intsurv::bootstrap_sample(case2_ind)
                    )
            };
            intsurv::CoxphCure boot_obj {
                intsurv::subset(obj, boot_ind)
            };
            boot_obj.control_.set_verbose(0);
            boot_obj.surv_obj_.control_.set_verbose(0);
            boot_obj.cure_obj_.control_.set_verbose(0);
            // fit the bootstarp sample
            boot_obj.fit();
            boot_surv_coef_mat.col(i) = boot_obj.surv_coef_;
            boot_cure_coef_mat.col(i) = boot_obj.cure_coef_;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = intsurv::arma2rvec(obj.surv_coef_),
        Rcpp::Named("cure_coef") = intsurv::arma2rvec(obj.cure_coef_),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = intsurv::arma2rvec(obj.unique_time_),
            Rcpp::Named("h0") = intsurv::arma2rvec(obj.h0_est_),
            Rcpp::Named("H0") = intsurv::arma2rvec(obj.H0_est_),
            Rcpp::Named("S0") = intsurv::arma2rvec(obj.S0_est_)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("surv_xBeta") = intsurv::arma2rvec(obj.surv_xbeta_),
            Rcpp::Named("cure_xBeta") = intsurv::arma2rvec(obj.cure_xbeta_),
            Rcpp::Named("susceptible_prob") =
            intsurv::arma2rvec(obj.susceptible_prob_),
            Rcpp::Named("estep_cured") =
            intsurv::arma2rvec(obj.estep_cured_),
            Rcpp::Named("estep_susceptible") =
            intsurv::arma2rvec(obj.estep_susceptible_)
            ),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("surv_offset") =
            intsurv::arma2rvec(obj.surv_obj_.control_.offset_),
            Rcpp::Named("cure_offset") =
            intsurv::arma2rvec(obj.cure_obj_.control_.offset_),
            Rcpp::Named("nObs") = obj.n_obs_,
            Rcpp::Named("nEvent") = obj.n_event_,
            Rcpp::Named("coef_df") = obj.coef_df_,
            Rcpp::Named("negLogL") = obj.neg_ll_,
            Rcpp::Named("c_index") = obj.c_index_,
            Rcpp::Named("aic") = obj.aic_,
            Rcpp::Named("bic1") = obj.bic1_,
            Rcpp::Named("bic2") = obj.bic2_
            ),
        Rcpp::Named("bootstrap") = Rcpp::List::create(
            Rcpp::Named("B") = bootstrap,
            Rcpp::Named("surv_coef_mat") = boot_surv_coef_mat.t(),
            Rcpp::Named("cure_coef_mat") = boot_cure_coef_mat.t()
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.n_iter_
            )
        );
}


// fit regularized Cox cure rate model by EM algorithm,
// where the M-step utilized CMD algoritm
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure_reg(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& surv_x,
    const arma::mat& cure_x,
    const bool cure_intercept,
    const double surv_l1_lambda,
    const double surv_l2_lambda,
    const arma::vec& surv_penalty_factor,
    const double cure_l1_lambda,
    const double cure_l2_lambda,
    const arma::vec& cure_penalty_factor,
    const unsigned long cv_nfolds,
    const arma::vec& surv_start,
    const arma::vec& cure_start,
    const arma::vec& surv_offset,
    const arma::vec& cure_offset,
    const bool surv_standardize = true,
    const bool cure_standardize = true,
    const bool surv_varying_active = true,
    const bool cure_varying_active = true,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const unsigned int surv_max_iter = 100,
    const double surv_epsilon = 1e-4,
    const unsigned int cure_max_iter = 100,
    const double cure_epsilon = 1e-4,
    const unsigned int tail_completion = 1,
    const double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control0 { max_iter, epsilon };
    control0.cure(tail_completion, tail_tau)->
        set_verbose(verbose);
    intsurv::Control surv_control {
        surv_max_iter, surv_epsilon, surv_standardize,
        intsurv::less_verbose(verbose, 3)
    };
    surv_control.set_start(surv_start)->
        set_offset(surv_offset)->
        net(surv_penalty_factor, surv_varying_active)->
        net_fit(surv_l1_lambda, surv_l2_lambda);
    intsurv::Control cure_control {
        cure_max_iter, cure_epsilon, cure_standardize,
        intsurv::less_verbose(verbose, 3)
    };
    cure_control.logistic(cure_intercept, pmin)->
        net(cure_penalty_factor, cure_varying_active)->
        net_fit(cure_l1_lambda, cure_l2_lambda)->
        set_start(cure_start)->
        set_offset(cure_offset);
    // define object
    intsurv::CoxphCure obj {
        time, event, surv_x, cure_x,
        control0, surv_control, cure_control
    };
    obj.net_fit();
    // cross-validation
    arma::vec cv_vec;
    if (cv_nfolds > 1) {
        cv_vec = intsurv::cv_coxph_cure_reg(obj, cv_nfolds);
    }
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = intsurv::arma2rvec(obj.surv_coef_),
        Rcpp::Named("cure_coef") = intsurv::arma2rvec(obj.cure_coef_),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = intsurv::arma2rvec(obj.unique_time_),
            Rcpp::Named("h0_est") = intsurv::arma2rvec(obj.h0_est_),
            Rcpp::Named("H0_est") = intsurv::arma2rvec(obj.H0_est_),
            Rcpp::Named("S0_est") = intsurv::arma2rvec(obj.S0_est_)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("surv_xBeta") = intsurv::arma2rvec(obj.surv_xbeta_),
            Rcpp::Named("cure_xBeta") = intsurv::arma2rvec(obj.cure_xbeta_),
            Rcpp::Named("susceptible_prob") =
            intsurv::arma2rvec(obj.susceptible_prob_),
            Rcpp::Named("estep_cured") =
            intsurv::arma2rvec(obj.estep_cured_),
            Rcpp::Named("estep_susceptible") =
            intsurv::arma2rvec(obj.estep_susceptible_)
            ),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.n_obs_,
            Rcpp::Named("nEvent") = obj.n_event_,
            Rcpp::Named("coef_df") = obj.coef_df_,
            Rcpp::Named("negLogL") = obj.neg_ll_,
            Rcpp::Named("c_index") = obj.c_index_,
            Rcpp::Named("aic") = obj.aic_,
            Rcpp::Named("bic1") = obj.bic1_,
            Rcpp::Named("bic2") = obj.bic2_,
            Rcpp::Named("cv_logL") = intsurv::arma2rvec(cv_vec)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("surv_l1_lambda_max") = obj.surv_obj_.l1_lambda_max_,
            Rcpp::Named("surv_l1_lambda") = obj.surv_obj_.control_.l1_lambda_,
            Rcpp::Named("surv_l2_lambda") = obj.surv_obj_.control_.l2_lambda_,
            Rcpp::Named("surv_penalty_factor") =
            intsurv::arma2rvec(obj.surv_obj_.control_.penalty_factor_),
            Rcpp::Named("cure_l1_lambda_max") = obj.cure_obj_.l1_lambda_max_,
            Rcpp::Named("cure_l1_lambda") = obj.cure_obj_.control_.l1_lambda_,
            Rcpp::Named("cure_l2_lambda") = obj.cure_obj_.control_.l2_lambda_,
            Rcpp::Named("cure_penalty_factor") =
            intsurv::arma2rvec(obj.cure_obj_.control_.penalty_factor_)
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.n_iter_
            )
        );
}


// variable selection for Cox cure rate model by EM algorithm,
// where the M-step utilized CMD algoritm
// for a sequence of lambda's
// lambda * (penalty_factor * alpha * lasso + (1 - alpha) / 2 * ridge)
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure_vs(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& surv_x,
    const arma::mat& cure_x,
    const bool cure_intercept,
    const arma::vec& surv_lambda,
    const double surv_alpha,
    const unsigned int surv_nlambda,
    const double surv_lambda_min_ratio,
    const arma::vec& surv_penalty_factor,
    const arma::vec& cure_lambda,
    const double cure_alpha,
    const unsigned int cure_nlambda,
    const double cure_lambda_min_ratio,
    const arma::vec& cure_penalty_factor,
    const unsigned long cv_nfolds,
    const arma::vec& surv_start,
    const arma::vec& cure_start,
    const arma::vec& surv_offset,
    const arma::vec& cure_offset,
    const bool surv_standardize = true,
    const bool cure_standardize = true,
    const bool surv_varying_active = true,
    const bool cure_varying_active = true,
    const unsigned int max_iter = 200,
    const double epsilon = 1e-4,
    const unsigned int surv_max_iter = 100,
    const double surv_epsilon = 1e-4,
    const unsigned int cure_max_iter = 100,
    const double cure_epsilon = 1e-4,
    const unsigned int tail_completion = 1,
    const double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int verbose = 0
    )
{
    intsurv::Control control0 { max_iter, epsilon };
    control0.cure(tail_completion, tail_tau)->
        set_verbose(verbose);
    intsurv::Control surv_control {
        surv_max_iter, surv_epsilon, surv_standardize,
        intsurv::less_verbose(verbose, 3)
    };
    surv_control.set_start(surv_start)->
        set_offset(surv_offset)->
        net(surv_penalty_factor, surv_varying_active)->
        net_path(surv_nlambda, surv_lambda_min_ratio, surv_alpha, surv_lambda);
    intsurv::Control cure_control {
        cure_max_iter, cure_epsilon, cure_standardize,
        intsurv::less_verbose(verbose, 3)
    };
    cure_control.logistic(cure_intercept, pmin)->
        net(cure_penalty_factor, cure_varying_active)->
        net_path(cure_nlambda, cure_lambda_min_ratio, cure_alpha, cure_lambda)->
        set_start(cure_start)->
        set_offset(cure_offset);
    // define object
    intsurv::CoxphCure obj {
        time, event, surv_x, cure_x,
        control0, surv_control, cure_control
    };
    obj.surv_obj_.set_penalty_factor();
    obj.cure_obj_.set_penalty_factor();
    obj.surv_obj_.set_l1_lambda_max();
    obj.cure_obj_.set_l1_lambda_max();
    // construct lambda sequence
    arma::vec surv_lambda_seq, cure_lambda_seq;
    if (surv_lambda.is_empty()) {
        const double surv_lambda_max {
            obj.surv_obj_.l1_lambda_max_ / std::max(surv_alpha, 1e-2)
        };
        double log_lambda_max { std::log(surv_lambda_max) };
        surv_lambda_seq = arma::exp(
            arma::linspace(log_lambda_max,
                           log_lambda_max + std::log(surv_lambda_min_ratio),
                           surv_nlambda)
            );
    } else {
        // take unique lambda and sort descendingly
        surv_lambda_seq = arma::reverse(arma::unique(surv_lambda));
    }
    if (cure_lambda.is_empty()) {
        const double cure_lambda_max {
            obj.cure_obj_.l1_lambda_max_ / std::max(cure_alpha, 1e-2)
        };
        double log_lambda_max { std::log(cure_lambda_max) };
        cure_lambda_seq = arma::exp(
            arma::linspace(log_lambda_max,
                           log_lambda_max + std::log(cure_lambda_min_ratio),
                           cure_nlambda)
            );
    } else {
        // take unique lambda and sort descendingly
        cure_lambda_seq = arma::reverse(arma::unique(cure_lambda));
    }
    // get the length of lambdas
    const unsigned int n_surv_lambda { surv_lambda_seq.n_elem };
    const unsigned int n_cure_lambda { cure_lambda_seq.n_elem };
    const unsigned int n_lambda { n_surv_lambda * n_cure_lambda };
    // initialize the coef matrices
    const unsigned int surv_p { obj.surv_p_ };
    const unsigned int cure_p { obj.cure_p_ };
    arma::mat surv_coef_mat { arma::zeros(surv_p, n_lambda) };
    arma::mat cure_coef_mat { arma::zeros(cure_p, n_lambda) };
    arma::vec bic1 { arma::zeros(n_lambda) }, bic2 { bic1 }, aic { bic1 },
        coef_df { bic1 }, negLogL { bic1 };
    arma::mat lambda_mat { arma::zeros(n_lambda, 4) };
    arma::vec cv_loglik { arma::zeros(n_lambda) };
    // warm starts
    arma::vec surv_warm_start0 { surv_start };
    arma::vec cure_warm_start0 { cure_start };
    arma::vec surv_warm_start;
    arma::vec cure_warm_start;
    // for each lambda
    unsigned int iter {0};
    for (size_t i {0}; i < n_surv_lambda; ++i) {
        // get the specific lambda's
        double surv_l1_lambda { surv_lambda_seq(i) * surv_alpha };
        double surv_l2_lambda { surv_lambda_seq(i) * (1 - surv_alpha) / 2 };
        obj.surv_obj_.control_.net_fit(surv_l1_lambda, surv_l2_lambda);
        surv_warm_start = surv_warm_start0;
        cure_warm_start = cure_warm_start0;
        for (size_t j {0}; j < n_cure_lambda; ++j) {
            double cure_l1_lambda { cure_lambda_seq(j) * cure_alpha };
            double cure_l2_lambda { cure_lambda_seq(j) * (1 - cure_alpha) / 2 };
            obj.cure_obj_.control_.net_fit(surv_l1_lambda, surv_l2_lambda);
            obj.surv_obj_.control_.set_start(surv_warm_start);
            obj.cure_obj_.control_.set_start(cure_warm_start);
            // model-fitting
            obj.net_fit();
            // cross-validation
            arma::vec cv_vec { arma::datum::nan };
            if (cv_nfolds > 1) {
                cv_vec = intsurv::cv_coxph_cure_reg(obj, cv_nfolds);
            }
            // update starting value
            surv_warm_start = obj.surv_coef_;
            cure_warm_start = obj.cure_coef_;
            if (j == 0) {
                // save starting value for next i
                surv_warm_start0 = obj.surv_coef_;
                cure_warm_start0 = obj.cure_coef_;
            }
            // store results
            surv_coef_mat.col(iter) = obj.surv_coef_;
            cure_coef_mat.col(iter) = obj.cure_coef_;
            // surv_en_coef_mat.col(iter) = obj.surv_en_coef_;
            // cure_en_coef_mat.col(iter) = obj.cure_en_coef_;
            aic(iter) = obj.aic_;
            bic1(iter) = obj.bic1_;
            bic2(iter) = obj.bic2_;
            coef_df(iter) = obj.coef_df_;
            negLogL(iter) = obj.neg_ll_;
            lambda_mat(iter, 0) = surv_l1_lambda;
            lambda_mat(iter, 1) = surv_l2_lambda;
            lambda_mat(iter, 2) = cure_l1_lambda;
            lambda_mat(iter, 3) = cure_l2_lambda;
            cv_loglik(iter) = arma::mean(cv_vec);
            // update iterators
            iter++;
        }
    }
    // return results in a list
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = surv_coef_mat.t(),
        Rcpp::Named("cure_coef") = cure_coef_mat.t(),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.n_obs_,
            Rcpp::Named("nEvent") = obj.n_event_,
            Rcpp::Named("coef_df") = intsurv::arma2rvec(coef_df),
            Rcpp::Named("negLogL") = intsurv::arma2rvec(negLogL),
            Rcpp::Named("aic") = intsurv::arma2rvec(aic),
            Rcpp::Named("bic1") = intsurv::arma2rvec(bic1),
            Rcpp::Named("bic2") = intsurv::arma2rvec(bic2),
            Rcpp::Named("cv_logL") = intsurv::arma2rvec(cv_loglik)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("lambda_mat") = lambda_mat,
            Rcpp::Named("surv_alpha") = surv_alpha,
            Rcpp::Named("cure_alpha") = cure_alpha,
            Rcpp::Named("surv_l1_lambda_max") = obj.surv_obj_.l1_lambda_max_,
            Rcpp::Named("cure_l1_lambda_max") = obj.cure_obj_.l1_lambda_max_,
            Rcpp::Named("surv_penalty_factor") =
            intsurv::arma2rvec(obj.surv_obj_.control_.penalty_factor_),
            Rcpp::Named("cure_penalty_factor") =
            intsurv::arma2rvec(obj.cure_obj_.control_.penalty_factor_)
            )
        );
}
