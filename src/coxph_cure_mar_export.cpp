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
#include <intsurv.h>

// Cox cure model with uncertain events without regularization
// [[Rcpp::export]]
Rcpp::List rcpp_coxph_cure_mar(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const arma::mat& mar_x,
    const bool cure_intercept = true,
    const bool mar_intercept = true,
    const unsigned int bootstrap = 0,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const arma::vec& mar_start = 0,
    const arma::vec& cox_offset = 0,
    const arma::vec& cure_offset = 0,
    const arma::vec& mar_offset = 0,
    const bool cox_standardize = true,
    const bool cure_standardize = true,
    const bool mar_standardize = true,
    const unsigned int em_max_iter = 300,
    const double em_rel_tol = 1e-5,
    const unsigned int cox_mstep_max_iter = 100,
    const double cox_mstep_rel_tol = 1e-5,
    const unsigned int cure_mstep_max_iter = 100,
    const double cure_mstep_rel_tol = 1e-5,
    const unsigned int mar_mstep_max_iter = 100,
    const double mar_mstep_rel_tol = 1e-5,
    const unsigned int tail_completion = 1,
    double tail_tau = -1,
    const double pmin = 1e-5,
    const unsigned int early_stop = 0,
    const unsigned int verbose = 0
    )
{
    // define object
    Intsurv::CoxphCureMar obj {
        time, event, cox_x, cure_x, mar_x,
        cure_intercept, mar_intercept,
        cox_standardize, cure_standardize, mar_standardize,
        cox_offset, cure_offset, mar_offset
    };
    // model-fitting
    obj.fit(cox_start, cure_start, mar_start,
            em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol,
            mar_mstep_max_iter, mar_mstep_rel_tol,
            tail_completion, tail_tau,
            pmin, early_stop, verbose
        );
    // // initialize bootstrap estimates
    // arma::mat boot_cox_coef_mat, boot_cure_coef_mat;
    // if (bootstrap > 0) {
    //     boot_cox_coef_mat = arma::zeros(obj.cox_coef_.n_elem, bootstrap);
    //     boot_cure_coef_mat = arma::zeros(obj.cure_coef_.n_elem, bootstrap);
    //     arma::vec event0na { event };
    //     const double const4na { 0.5 };
    //     event0na.replace(arma::datum::nan, const4na);
    //     arma::uvec case1_ind = arma::find(event0na > const4na);
    //     arma::uvec case2_ind = arma::find(event0na < const4na);
    //     arma::uvec case3_ind = arma::find(event0na == const4na);
    //     for (size_t i {0}; i < bootstrap; ++i) {
    //         // generate a bootstrap sample
    //         arma::uvec boot_ind {
    //             Intsurv::vec_union(
    //                 Intsurv::bootstrap_sample(case1_ind),
    //                 Intsurv::bootstrap_sample(case2_ind)
    //                 )
    //         };
    //         boot_ind = Intsurv::vec_union(
    //             boot_ind, Intsurv::bootstrap_sample(case3_ind));
    //         Intsurv::CoxphCureMcar boot_obj {
    //             time.elem(boot_ind),
    //             event.elem(boot_ind),
    //             cox_x.rows(boot_ind),
    //             cure_x.rows(boot_ind),
    //             cure_intercept,
    //             cox_standardize,
    //             cure_standardize,
    //             cox_offset.elem(boot_ind),
    //             cure_offset.elem(boot_ind)
    //         };
    //         boot_obj.fit(cox_start, cure_start,
    //                      em_max_iter, em_rel_tol,
    //                      cox_mstep_max_iter, cox_mstep_rel_tol,
    //                      cure_mstep_max_iter, cure_mstep_rel_tol,
    //                      tail_completion, tail_tau,
    //                      pmin, early_stop, 0);
    //         boot_cox_coef_mat.col(i) = boot_obj.cox_coef_;
    //         boot_cure_coef_mat.col(i) = boot_obj.cure_coef_;
    //     }
    // }
    return Rcpp::List::create(
        Rcpp::Named("surv_coef") = Intsurv::arma2rvec(obj.cox_coef_),
        Rcpp::Named("cure_coef") = Intsurv::arma2rvec(obj.cure_coef_),
        Rcpp::Named("mar_coef") = Intsurv::arma2rvec(obj.mar_coef_),
        Rcpp::Named("baseline") = Rcpp::List::create(
            Rcpp::Named("time") = Intsurv::arma2rvec(obj.unique_time_),
            Rcpp::Named("h0") = Intsurv::arma2rvec(obj.h0_est_),
            Rcpp::Named("H0") = Intsurv::arma2rvec(obj.H0_est_),
            Rcpp::Named("S0") = Intsurv::arma2rvec(obj.S0_est_),
            Rcpp::Named("hc") = Intsurv::arma2rvec(obj.hc_est_),
            Rcpp::Named("Hc") = Intsurv::arma2rvec(obj.Hc_est_),
            Rcpp::Named("Sc") = Intsurv::arma2rvec(obj.Sc_est_)
            ),
        Rcpp::Named("fitted") = Rcpp::List::create(
            Rcpp::Named("surv_xBeta") = Intsurv::arma2rvec(obj.cox_xbeta_),
            Rcpp::Named("cure_xBeta") = Intsurv::arma2rvec(obj.cure_xbeta_),
            Rcpp::Named("susceptible_prob") =
            Intsurv::arma2rvec(obj.susceptible_prob_),
            Rcpp::Named("estep_cured") = Intsurv::arma2rvec(obj.estep_cured_),
            Rcpp::Named("estep_event") = Intsurv::arma2rvec(obj.estep_event_),
            Rcpp::Named("estep_censor") = Intsurv::arma2rvec(obj.estep_censor_)
            ),
        Rcpp::Named("model") = Rcpp::List::create(
            Rcpp::Named("nObs") = obj.cox_obj_.n_obs_,
            Rcpp::Named("nCase1") = obj.n_case1_,
            Rcpp::Named("coef_df") = obj.coef_df_,
            Rcpp::Named("negLogL") = obj.neg_ll_,
            Rcpp::Named("c_index") = obj.c_index_,
            Rcpp::Named("aic") = obj.aic_,
            Rcpp::Named("bic1") = obj.bic1_,
            Rcpp::Named("bic2") = obj.bic2_
            ),
        // Rcpp::Named("bootstrap") = Rcpp::List::create(
        //     Rcpp::Named("B") = bootstrap,
        //     Rcpp::Named("surv_coef_mat") = boot_cox_coef_mat.t(),
        //     Rcpp::Named("cure_coef_mat") = boot_cure_coef_mat.t()
        //     ),
        Rcpp::Named("convergence") = Rcpp::List::create(
            Rcpp::Named("num_iter") = obj.n_iter_
            )
        );
}


