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

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <intsurv.h>

// Cox cure model with uncertain events without regularization
// [[Rcpp::export]]
Rcpp::List coxph_cure_uncer(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const unsigned int& bootstrap = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 300,
    const double& em_rel_tol = 1e-5,
    const unsigned int& cox_mstep_max_iter = 100,
    const double& cox_mstep_rel_tol = 1e-5,
    const unsigned int& cure_mstep_max_iter = 100,
    const double& cure_mstep_rel_tol = 1e-5,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool& spline_start = false,
    const unsigned int& iSpline_num_knots = 3,
    const unsigned int& iSpline_degree = 2,
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    // define object
    Intsurv::CoxphCureUncer obj {
        time, event, cox_x, cure_x, cure_intercept,
        cox_standardize, cure_standardize
    };
    // model-fitting
    obj.fit(cox_start, cure_start,
            em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol,
            spline_start, iSpline_num_knots, iSpline_degree,
            tail_completion, tail_tau,
            pmin, early_stop, verbose
        );
    // initialize bootstrap estimates
    arma::mat boot_cox_coef_mat, boot_cure_coef_mat;
    if (bootstrap > 0) {
        boot_cox_coef_mat = arma::zeros(obj.cox_coef.n_elem, bootstrap);
        boot_cure_coef_mat = arma::zeros(obj.cure_coef.n_elem, bootstrap);
        arma::vec event0na { event };
        const double const4na { 0.5 };
        event0na.replace(arma::datum::nan, const4na);
        arma::uvec case1_ind = arma::find(event0na > const4na);
        arma::uvec case2_ind = arma::find(event0na < const4na);
        arma::uvec case3_ind = arma::find(event0na == const4na);
        for (size_t i {0}; i < bootstrap; ++i) {
            // generate a bootstrap sample
            arma::uvec boot_ind {
                Intsurv::vec_union(
                    Intsurv::bootstrap_sample(case1_ind),
                    Intsurv::bootstrap_sample(case2_ind)
                    )
            };
            boot_ind = Intsurv::vec_union(
                boot_ind, Intsurv::bootstrap_sample(case3_ind));
            Intsurv::CoxphCureUncer boot_obj {
                time.elem(boot_ind),
                event.elem(boot_ind),
                cox_x.rows(boot_ind),
                cure_x.rows(boot_ind),
                cure_intercept,
                cox_standardize,
                cure_standardize
            };
            boot_obj.fit(cox_start, cure_start,
                         em_max_iter, em_rel_tol,
                         cox_mstep_max_iter, cox_mstep_rel_tol,
                         cure_mstep_max_iter, cure_mstep_rel_tol,
                         spline_start, iSpline_num_knots, iSpline_degree,
                         tail_completion, tail_tau,
                         pmin, early_stop, 0);
            boot_cox_coef_mat.col(i) = boot_obj.cox_coef;
            boot_cure_coef_mat.col(i) = boot_obj.cure_coef;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est),
            Rcpp::Named("hc") = Intsurv::arma2rvec(obj.hc_est),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(obj.Hc_est),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(obj.Sc_est)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("surv_xBeta") = Intsurv::arma2rvec(obj.cox_xBeta),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xBeta),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob),
            Rcpp::Named("estep_cured") = Intsurv::arma2rvec(obj.estep_cured),
            Rcpp::Named("estep_event") = Intsurv::arma2rvec(obj.estep_event),
            Rcpp::Named("estep_censor") = Intsurv::arma2rvec(obj.estep_censor)
            ),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.nObs,
            Rcpp::Named("nEvent") = obj.nEvent,
            Rcpp::Named("coef_df") = obj.coef_df,
            Rcpp::Named("negLogL") = obj.negLogL,
            Rcpp::Named("c_index") = obj.c_index,
            Rcpp::Named("aic") = obj.aic,
            Rcpp::Named("bic1") = obj.bic1,
            Rcpp::Named("bic2") = obj.bic2
            ),
        Rcpp::Named("bootstrap") = Rcpp::List::create(
            Rcpp::Named("B") = bootstrap,
            Rcpp::Named("surv_coef_mat") = boot_cox_coef_mat.t(),
            Rcpp::Named("cure_coef_mat") = boot_cure_coef_mat.t()
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.num_iter
            )
        );
}


// fit regularized Cox cure rate model with uncertain events
// by EM algorithm, where the M-step utilized CMD algoritm
// [[Rcpp::export]]
Rcpp::List coxph_cure_uncer_reg(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const double& cox_l1_lambda = 0,
    const double& cox_l2_lambda = 0,
    const arma::vec& cox_l1_penalty_factor = 0,
    const double& cure_l1_lambda = 0,
    const double& cure_l2_lambda = 0,
    const arma::vec& cure_l1_penalty_factor = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 500,
    const double& em_rel_tol = 1e-4,
    const unsigned int& cox_mstep_max_iter = 200,
    const double& cox_mstep_rel_tol = 1e-4,
    const unsigned int& cure_mstep_max_iter = 200,
    const double& cure_mstep_rel_tol = 1e-4,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool& spline_start = false,
    const unsigned int& iSpline_num_knots = 3,
    const unsigned int& iSpline_degree = 2,
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    // define object
    Intsurv::CoxphCureUncer obj {
        time, event, cox_x, cure_x, cure_intercept,
        cox_standardize, cure_standardize
    };
    // model-fitting
    obj.regularized_fit(
        cox_l1_lambda, cox_l2_lambda,
        cure_l1_lambda, cure_l2_lambda,
        cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start, em_max_iter, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol,
        spline_start, iSpline_num_knots, iSpline_degree,
        tail_completion, tail_tau,
        pmin, early_stop, verbose
        );
    // return results in a list
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef),
        Rcpp::Named("cox_en_coef") = Intsurv::arma2rvec(obj.cox_en_coef),
        Rcpp::Named("cure_en_coef") = Intsurv::arma2rvec(obj.cure_en_coef),

        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est),
            Rcpp::Named("hc") = Intsurv::arma2rvec(obj.hc_est),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(obj.Hc_est),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(obj.Sc_est)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("cox_xBeta") = Intsurv::arma2rvec(obj.cox_xBeta),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xBeta),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob),
            Rcpp::Named("estep_cured") = Intsurv::arma2rvec(obj.estep_cured),
            Rcpp::Named("estep_event") = Intsurv::arma2rvec(obj.estep_event),
            Rcpp::Named("estep_censor") = Intsurv::arma2rvec(obj.estep_censor)
            ),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.nObs,
            Rcpp::Named("nEvent") = obj.nEvent,
            Rcpp::Named("coef_df") = obj.coef_df,
            Rcpp::Named("negLogL") = obj.negLogL,
            Rcpp::Named("c_index") = obj.c_index,
            Rcpp::Named("aic") = obj.aic,
            Rcpp::Named("bic1") = obj.bic1,
            Rcpp::Named("bic2") = obj.bic2
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("cox_l1_lambda_max") = obj.cox_l1_lambda_max,
            Rcpp::Named("cox_l1_lambda") = obj.cox_l1_lambda,
            Rcpp::Named("cox_l2_lambda") = obj.cox_l2_lambda,
            Rcpp::Named("cox_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cox_l1_penalty_factor),
            Rcpp::Named("cure_l1_lambda_max") = obj.cure_l1_lambda_max,
            Rcpp::Named("cure_l1_lambda") = obj.cure_l1_lambda,
            Rcpp::Named("cure_l2_lambda") = obj.cure_l2_lambda,
            Rcpp::Named("cure_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cure_l1_penalty_factor)
            ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.num_iter
            )
        );
}


// variable selection for the regularized Cox cure rate model with uncertain
// events by EM algorithm, where the M-step utilized CMD algoritm
// for a sequence of lambda's
// lambda * (penalty_factor * alpha * lasso + (1 - alpha) / 2 * ridge)
// [[Rcpp::export]]
Rcpp::List coxph_cure_uncer_vs(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_lambda = 0,
    const double& cox_alpha = 1,
    const unsigned int& cox_nlambda = 1,
    const double& cox_lambda_min_ratio = 1e-4,
    const arma::vec& cox_l1_penalty_factor = 0,
    const arma::vec& cure_lambda = 0,
    const double& cure_alpha = 1,
    const unsigned int& cure_nlambda = 1,
    const double& cure_lambda_min_ratio = 1e-4,
    const arma::vec& cure_l1_penalty_factor = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 500,
    const double& em_rel_tol = 1e-4,
    const unsigned int& cox_mstep_max_iter = 200,
    const double& cox_mstep_rel_tol = 1e-4,
    const unsigned int& cure_mstep_max_iter = 200,
    const double& cure_mstep_rel_tol = 1e-4,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool& spline_start = false,
    const unsigned int& iSpline_num_knots = 3,
    const unsigned int& iSpline_degree = 2,
    const unsigned int& tail_completion = 1,
    double tail_tau = -1,
    const double& pmin = 1e-5,
    const unsigned int& early_stop = 0,
    const unsigned int& verbose = 0
    )
{
    // define object
    Intsurv::CoxphCureUncer obj {
        time, event, cox_x, cure_x, cure_intercept,
        cox_standardize, cure_standardize
    };
    // get the maximum lambdas by setting em_max_iter = 0
    obj.regularized_fit(
        0, 0, 0, 0, cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start, 0, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol,
        spline_start, iSpline_num_knots, iSpline_degree,
        tail_completion, tail_tau,
        pmin, early_stop, 0
        );
    // already considered penalty factor
    const double cox_lambda_max {
        obj.cox_l1_lambda_max / std::max(cox_alpha, 1e-10)
    };
    const double cure_lambda_max {
        obj.cure_l1_lambda_max / std::max(cure_alpha, 1e-10)
    };
    // construct lambda sequence
    arma::vec cox_lambda_seq, cure_lambda_seq;
    // if cox_nlambda is not default value, generate the lambda sequence
    if (cox_nlambda > 1) {
        double log_lambda_max { std::log(cox_lambda_max) };
        cox_lambda_seq = arma::exp(
            arma::linspace(log_lambda_max,
                           log_lambda_max + std::log(cox_lambda_min_ratio),
                           cox_nlambda)
            );
    } else {
        // take unique lambda and sort descendingly
        cox_lambda_seq = arma::reverse(arma::unique(cox_lambda));
    }
    // if cure_nlambda is not default value, generate the lambda sequence
    if (cure_nlambda > 1) {
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
    const unsigned int n_cox_lambda { cox_lambda_seq.n_elem };
    const unsigned int n_cure_lambda { cure_lambda_seq.n_elem };
    const unsigned int n_lambda { n_cox_lambda * n_cure_lambda };

    // initialize the coef matrices
    const unsigned int cox_p { obj.get_cox_p() };
    const unsigned int cure_p { obj.get_cure_p() };
    arma::mat cox_coef_mat { arma::zeros(cox_p, n_lambda) };
    arma::mat cure_coef_mat { arma::zeros(cure_p, n_lambda) };
    arma::mat cox_en_coef_mat { arma::zeros(cox_p, n_lambda) };
    arma::mat cure_en_coef_mat { arma::zeros(cure_p, n_lambda) };
    arma::vec bic1 { arma::zeros(n_lambda) }, bic2 { bic1 }, aic { bic1 };
    arma::vec coef_df { bic1 }, negLogL { bic1 };
    arma::mat lambda_mat { arma::zeros(n_lambda, 4) };

    // warm starts
    arma::vec cox_warm_start0 { cox_start };
    arma::vec cure_warm_start0 { cure_start };
    arma::vec cox_warm_start;
    arma::vec cure_warm_start;

    // for each lambda
    unsigned int iter {0};
    for (size_t i {0}; i < n_cox_lambda; ++i) {
        // get the specific lambda's
        double cox_l1_lambda { cox_lambda_seq(i) * cox_alpha };
        double cox_l2_lambda { cox_lambda_seq(i) * (1 - cox_alpha) / 2 };
        cox_warm_start = cox_warm_start0;
        cure_warm_start = cure_warm_start0;
        for (size_t j {0}; j < n_cure_lambda; ++j) {
            double cure_l1_lambda { cure_lambda_seq(j) * cure_alpha };
            double cure_l2_lambda { cure_lambda_seq(j) * (1 - cure_alpha) / 2 };
            // model-fitting
            obj.regularized_fit(
                cox_l1_lambda, cox_l2_lambda,
                cure_l1_lambda, cure_l2_lambda,
                cox_l1_penalty_factor, cure_l1_penalty_factor,
                cox_warm_start, cure_warm_start,
                em_max_iter, em_rel_tol,
                cox_mstep_max_iter, cox_mstep_rel_tol,
                cure_mstep_max_iter, cure_mstep_rel_tol,
                spline_start, iSpline_num_knots, iSpline_degree,
                tail_completion, tail_tau,
                pmin, early_stop, verbose
                );
            // update starting value
            cox_warm_start = obj.cox_coef;
            cure_warm_start = obj.cure_coef;
            if (j == 0) {
                // save starting value for next i
                cox_warm_start0 = obj.cox_coef;
                cure_warm_start0 = obj.cure_coef;
            }
            // store results
            cox_coef_mat.col(iter) = obj.cox_coef;
            cure_coef_mat.col(iter) = obj.cure_coef;
            cox_en_coef_mat.col(iter) = obj.cox_en_coef;
            cure_en_coef_mat.col(iter) = obj.cure_en_coef;
            aic(iter) = obj.aic;
            bic1(iter) = obj.bic1;
            bic2(iter) = obj.bic2;
            coef_df(iter) = obj.coef_df;
            negLogL(iter) = obj.negLogL;
            lambda_mat(iter, 0) = cox_l1_lambda;
            lambda_mat(iter, 1) = cox_l2_lambda;
            lambda_mat(iter, 2) = cure_l1_lambda;
            lambda_mat(iter, 3) = cure_l2_lambda;
            // update iterators
            iter++;
        }
    }
    // return results in a list
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = cox_coef_mat.t(),
        Rcpp::Named("cure_coef") = cure_coef_mat.t(),
        Rcpp::Named("surv_en_coef") = cox_en_coef_mat.t(),
        Rcpp::Named("cure_en_coef") = cure_en_coef_mat.t(),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.nObs,
            Rcpp::Named("nEvent") = obj.nEvent,
            Rcpp::Named("coef_df") = Intsurv::arma2rvec(coef_df),
            Rcpp::Named("negLogL") = Intsurv::arma2rvec(negLogL),
            Rcpp::Named("aic") = Intsurv::arma2rvec(aic),
            Rcpp::Named("bic1") = Intsurv::arma2rvec(bic1),
            Rcpp::Named("bic2") = Intsurv::arma2rvec(bic2)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("lambda_mat") = lambda_mat,
            Rcpp::Named("surv_alpha") = cox_alpha,
            Rcpp::Named("cure_alpha") = cure_alpha,
            Rcpp::Named("surv_l1_lambda_max") = obj.cox_l1_lambda_max,
            Rcpp::Named("cure_l1_lambda_max") = obj.cure_l1_lambda_max,
            Rcpp::Named("surv_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cox_l1_penalty_factor),
            Rcpp::Named("cure_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cure_l1_penalty_factor)
            )
        );
}
