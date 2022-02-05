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

// fit regular Cox cure rate model by EM algorithm
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const unsigned int bootstrap = 0,
    const bool firth = false,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const arma::vec& cox_offset = 0,
    const arma::vec& cure_offset = 0,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const unsigned int em_max_iter = 1000,
    const double em_rel_tol = 1e-4,
    const unsigned int cox_mstep_max_iter = 200,
    const double cox_mstep_rel_tol = 1e-4,
    const unsigned int cure_mstep_max_iter = 200,
    const double cure_mstep_rel_tol = 1e-4,
    const unsigned int tail_completion = 1,
    const double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int early_stop = 0,
    const unsigned int verbose = 0
    )
{
    // define object
    Intsurv::CoxphCure obj {
        time, event, cox_x, cure_x, cure_intercept,
        cox_standardize, cure_standardize,
        cox_offset, cure_offset
    };
    // model-fitting
    obj.fit(cox_start, cure_start,
            em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol,
            firth, tail_completion, tail_tau,
            pmin, early_stop, verbose);
    // initialize bootstrap estimates
    arma::mat boot_cox_coef_mat, boot_cure_coef_mat;
    if (bootstrap > 0) {
        boot_cox_coef_mat = arma::zeros(obj.cox_coef_.n_elem, bootstrap);
        boot_cure_coef_mat = arma::zeros(obj.cure_coef_.n_elem, bootstrap);
        arma::uvec case1_ind { arma::find(event > 0) };
        arma::uvec case2_ind { arma::find(event < 1) };
        for (size_t i {0}; i < bootstrap; ++i) {
            // generate a bootstrap sample
            arma::uvec boot_ind {
                Intsurv::vec_union(
                    Intsurv::bootstrap_sample(case1_ind),
                    Intsurv::bootstrap_sample(case2_ind)
                    )
            };
            Intsurv::CoxphCure boot_obj {
                time.elem(boot_ind),
                event.elem(boot_ind),
                cox_x.rows(boot_ind),
                cure_x.rows(boot_ind),
                cure_intercept,
                cox_standardize,
                cure_standardize,
                cox_offset.elem(boot_ind),
                cure_offset.elem(boot_ind),
            };
            // fit the bootstarp sample
            boot_obj.fit(cox_start, cure_start,
                         em_max_iter, em_rel_tol,
                         cox_mstep_max_iter, cox_mstep_rel_tol,
                         cure_mstep_max_iter, cure_mstep_rel_tol,
                         firth, tail_completion, tail_tau,
                         pmin, early_stop, 0);
            boot_cox_coef_mat.col(i) = boot_obj.cox_coef_;
            boot_cure_coef_mat.col(i) = boot_obj.cure_coef_;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = Intsurv::arma2rvec(obj.cox_coef_),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef_),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time_),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est_),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est_),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est_)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("surv_xBeta") = Intsurv::arma2rvec(obj.cox_xbeta_),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xbeta_),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob_),
            Rcpp::Named("estep_cured") =
            Intsurv::arma2rvec(obj.estep_cured_),
            Rcpp::Named("estep_susceptible") =
            Intsurv::arma2rvec(obj.estep_susceptible_)
            ),
        Rcpp::Named("model") = Rcpp::List::create(
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
            Rcpp::Named("surv_coef_mat") = boot_cox_coef_mat.t(),
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
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const double cox_l1_lambda = 0,
    const double cox_l2_lambda = 0,
    const arma::vec& cox_l1_penalty_factor = 0,
    const double cure_l1_lambda = 0,
    const double cure_l2_lambda = 0,
    const arma::vec& cure_l1_penalty_factor = 0,
    const unsigned long cv_nfolds = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const arma::vec& cox_offset = 0,
    const arma::vec& cure_offset = 0,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const unsigned int em_max_iter = 200,
    const double em_rel_tol = 1e-5,
    const unsigned int cox_mstep_max_iter = 100,
    const double cox_mstep_rel_tol = 1e-4,
    const unsigned int cure_mstep_max_iter = 100,
    const double cure_mstep_rel_tol = 1e-4,
    const unsigned int tail_completion = 1,
    const double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int early_stop = 0,
    const unsigned int verbose = 0
    )
{
    Intsurv::CoxphCure obj {
        time, event, cox_x, cure_x, cure_intercept,
        cox_standardize, cure_standardize,
        cox_offset, cure_offset
    };
    obj.regularized_fit(
        cox_l1_lambda, cox_l2_lambda,
        cure_l1_lambda, cure_l2_lambda,
        cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start,
        em_max_iter, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol,
        tail_completion, tail_tau,
        pmin, early_stop, verbose
        );
    // cross-validation
    arma::vec cv_vec;
    if (cv_nfolds > 1) {
        cv_vec = Intsurv::cv_coxph_cure_reg(
            time, event, cox_x, cure_x, cure_intercept,
            cv_nfolds,
            cox_l1_lambda, cox_l2_lambda,
            cox_l1_penalty_factor,
            cure_l1_lambda, cure_l2_lambda,
            cure_l1_penalty_factor,
            cox_start, cure_start,
            cox_offset, cure_offset,
            cox_standardize, cure_standardize,
            em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol,
            tail_completion, tail_tau,
            pmin, early_stop, verbose
            );
    }
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = Intsurv::arma2rvec(obj.cox_coef_),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef_),
        // Rcpp::Named("cox_en_coef") = Intsurv::arma2rvec(obj.cox_en_coef_),
        // Rcpp::Named("cure_en_coef") = Intsurv::arma2rvec(obj.cure_en_coef_),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time_),
            Rcpp::Named("h0_est") = Intsurv::arma2rvec(obj.h0_est_),
            Rcpp::Named("H0_est") = Intsurv::arma2rvec(obj.H0_est_),
            Rcpp::Named("S0_est") = Intsurv::arma2rvec(obj.S0_est_)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("cox_xBeta") = Intsurv::arma2rvec(obj.cox_xbeta_),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xbeta_),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob_),
            Rcpp::Named("estep_cured") =
            Intsurv::arma2rvec(obj.estep_cured_),
            Rcpp::Named("estep_susceptible") =
            Intsurv::arma2rvec(obj.estep_susceptible_)
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
            Rcpp::Named("cv_logL") = Intsurv::arma2rvec(cv_vec)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("cox_l1_lambda_max") = obj.cox_l1_lambda_max_,
            Rcpp::Named("cox_l1_lambda") = obj.cox_l1_lambda_,
            Rcpp::Named("cox_l2_lambda") = obj.cox_l2_lambda_,
            Rcpp::Named("cox_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cox_l1_penalty_factor_),
            Rcpp::Named("cure_l1_lambda_max") = obj.cure_l1_lambda_max_,
            Rcpp::Named("cure_l1_lambda") = obj.cure_l1_lambda_,
            Rcpp::Named("cure_l2_lambda") = obj.cure_l2_lambda_,
            Rcpp::Named("cure_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cure_l1_penalty_factor_)
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
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_lambda = 0,
    const double cox_alpha = 1,
    const unsigned int cox_nlambda = 1,
    const double cox_lambda_min_ratio = 1e-4,
    const arma::vec& cox_l1_penalty_factor = 0,
    const arma::vec& cure_lambda = 0,
    const double cure_alpha = 1,
    const unsigned int cure_nlambda = 1,
    const double cure_lambda_min_ratio = 1e-4,
    const arma::vec& cure_l1_penalty_factor = 0,
    const unsigned long cv_nfolds = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const arma::vec& cox_offset = 0,
    const arma::vec& cure_offset = 0,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const unsigned int em_max_iter = 500,
    const double em_rel_tol = 1e-4,
    const unsigned int cox_mstep_max_iter = 200,
    const double cox_mstep_rel_tol = 1e-4,
    const unsigned int cure_mstep_max_iter = 200,
    const double cure_mstep_rel_tol = 1e-4,
    const unsigned int tail_completion = 1,
    const double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int early_stop = 0,
    const unsigned int verbose = 0
    )
{
    // define object
    Intsurv::CoxphCure obj {
        time, event, cox_x, cure_x, cure_intercept,
        cox_standardize, cure_standardize,
        cox_offset, cure_offset
    };
    // get the maximum lambdas by setting em_max_iter = 0
    obj.regularized_fit(0, 0, 0, 0,
                        cox_l1_penalty_factor, cure_l1_penalty_factor,
                        cox_start, cure_start,
                        0, em_rel_tol,
                        cox_mstep_max_iter, cox_mstep_rel_tol,
                        cure_mstep_max_iter, cure_mstep_rel_tol,
                        tail_completion, tail_tau,
                        pmin, early_stop, verbose);
    // already considered penalty factor
    const double cox_lambda_max {
        obj.cox_l1_lambda_max_ / std::max(cox_alpha, 1e-10)
    };
    const double cure_lambda_max {
        obj.cure_l1_lambda_max_ / std::max(cure_alpha, 1e-10)
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
    const unsigned int cox_p { obj.cox_p_ };
    const unsigned int cure_p { obj.cure_p_ };
    arma::mat cox_coef_mat { arma::zeros(cox_p, n_lambda) };
    arma::mat cure_coef_mat { arma::zeros(cure_p, n_lambda) };
    // arma::mat cox_en_coef_mat { arma::zeros(cox_p, n_lambda) };
    // arma::mat cure_en_coef_mat { arma::zeros(cure_p, n_lambda) };
    arma::vec bic1 { arma::zeros(n_lambda) }, bic2 { bic1 }, aic { bic1 };
    arma::vec coef_df { bic1 }, negLogL { bic1 };
    arma::mat lambda_mat { arma::zeros(n_lambda, 4) };
    arma::vec cv_loglik { arma::zeros(n_lambda) };

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
                tail_completion, tail_tau,
                pmin, early_stop, verbose
                );
            // cross-validation
            arma::vec cv_vec;
            if (cv_nfolds > 1) {
                cv_vec = Intsurv::cv_coxph_cure_reg(
                    time, event, cox_x, cure_x,
                    cure_intercept, cv_nfolds,
                    cox_l1_lambda, cox_l2_lambda,
                    cox_l1_penalty_factor,
                    cure_l1_lambda, cure_l2_lambda,
                    cure_l1_penalty_factor,
                    cox_warm_start, cure_warm_start,
                    cox_offset, cure_offset,
                    cox_standardize, cure_standardize,
                    em_max_iter, em_rel_tol,
                    cox_mstep_max_iter, cox_mstep_rel_tol,
                    cure_mstep_max_iter, cure_mstep_rel_tol,
                    tail_completion, tail_tau,
                    pmin, early_stop, 0
                    );
            }
            // update starting value
            cox_warm_start = obj.cox_coef_;
            cure_warm_start = obj.cure_coef_;
            if (j == 0) {
                // save starting value for next i
                cox_warm_start0 = obj.cox_coef_;
                cure_warm_start0 = obj.cure_coef_;
            }
            // store results
            cox_coef_mat.col(iter) = obj.cox_coef_;
            cure_coef_mat.col(iter) = obj.cure_coef_;
            // cox_en_coef_mat.col(iter) = obj.cox_en_coef_;
            // cure_en_coef_mat.col(iter) = obj.cure_en_coef_;
            aic(iter) = obj.aic_;
            bic1(iter) = obj.bic1_;
            bic2(iter) = obj.bic2_;
            coef_df(iter) = obj.coef_df_;
            negLogL(iter) = obj.neg_ll_;
            lambda_mat(iter, 0) = cox_l1_lambda;
            lambda_mat(iter, 1) = cox_l2_lambda;
            lambda_mat(iter, 2) = cure_l1_lambda;
            lambda_mat(iter, 3) = cure_l2_lambda;
            cv_loglik(iter) = arma::sum(cv_vec);
            // update iterators
            iter++;
        }
    }
    // return results in a list
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = cox_coef_mat.t(),
        Rcpp::Named("cure_coef") = cure_coef_mat.t(),
        // Rcpp::Named("surv_en_coef") = cox_en_coef_mat.t(),
        // Rcpp::Named("cure_en_coef") = cure_en_coef_mat.t(),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.n_obs_,
            Rcpp::Named("nEvent") = obj.n_event_,
            Rcpp::Named("coef_df") = Intsurv::arma2rvec(coef_df),
            Rcpp::Named("negLogL") = Intsurv::arma2rvec(negLogL),
            Rcpp::Named("aic") = Intsurv::arma2rvec(aic),
            Rcpp::Named("bic1") = Intsurv::arma2rvec(bic1),
            Rcpp::Named("bic2") = Intsurv::arma2rvec(bic2),
            Rcpp::Named("cv_logL") = Intsurv::arma2rvec(cv_loglik)
            ),
        Rcpp::Named("penalty") = Rcpp::List::create(
            Rcpp::Named("lambda_mat") = lambda_mat,
            Rcpp::Named("surv_alpha") = cox_alpha,
            Rcpp::Named("cure_alpha") = cure_alpha,
            Rcpp::Named("surv_l1_lambda_max") = obj.cox_l1_lambda_max_,
            Rcpp::Named("cure_l1_lambda_max") = obj.cure_l1_lambda_max_,
            Rcpp::Named("surv_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cox_l1_penalty_factor_),
            Rcpp::Named("cure_l1_penalty_factor") =
            Intsurv::arma2rvec(obj.cure_l1_penalty_factor_)
            )
        );
}
