#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "coxph_cure.h"
#include "utils.h"


// fit regular Cox cure rate model by EM algorithm
// [[Rcpp::export]]
Rcpp::List coxph_cure(
    const arma::vec& time,
    const arma::vec& event,
    const arma::mat& cox_x,
    const arma::mat& cure_x,
    const bool cure_intercept = true,
    const arma::vec& cox_start = 0,
    const arma::vec& cure_start = 0,
    const unsigned int& em_max_iter = 1000,
    const double& em_rel_tol = 1e-3,
    const unsigned int& cox_mstep_max_iter = 30,
    const double& cox_mstep_rel_tol = 1e-3,
    const unsigned int& cure_mstep_max_iter = 30,
    const double& cure_mstep_rel_tol = 1e-6,
    const bool cox_standardize = true,
    const bool cure_standardize = true
    )
{
    Intsurv::CoxphCure obj {
        Intsurv::CoxphCure(time, event, cox_x, cure_x, cure_intercept,
                           cox_standardize, cure_standardize)
    };
    obj.fit(cox_start, cure_start, em_max_iter, em_rel_tol,
            cox_mstep_max_iter, cox_mstep_rel_tol,
            cure_mstep_max_iter, cure_mstep_rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = obj.cox_coef,
        Rcpp::Named("cure_coef") = obj.cure_coef,
        Rcpp::Named("negLogL") = obj.negLogL,
        Rcpp::Named("nObs") = obj.nObs
        );
}


// fit regularized Cox cure rate model by EM algorithm,
// where the M-step utilized CMD algoritm
// [[Rcpp::export]]
Rcpp::List coxph_cure_reg(
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
    const unsigned int& em_max_iter = 200,
    const double& em_rel_tol = 1e-4,
    const unsigned int& cox_mstep_max_iter = 30,
    const double& cox_mstep_rel_tol = 1e-3,
    const unsigned int& cure_mstep_max_iter = 30,
    const double& cure_mstep_rel_tol = 1e-3
    )
{
    Intsurv::CoxphCure obj {
        Intsurv::CoxphCure(time, event, cox_x, cure_x, cure_intercept)
    };
    obj.regularized_fit(
        cox_l1_lambda, cox_l2_lambda,
        cure_l1_lambda, cure_l2_lambda,
        cox_l1_penalty_factor, cure_l1_penalty_factor,
        cox_start, cure_start, em_max_iter, em_rel_tol,
        cox_mstep_max_iter, cox_mstep_rel_tol,
        cure_mstep_max_iter, cure_mstep_rel_tol
        );
    return Rcpp::List::create(
        Rcpp::Named("cox_coef") = obj.cox_coef,
        Rcpp::Named("cure_coef") = obj.cure_coef,
        Rcpp::Named("en_cox_coef") = obj.en_cox_coef,
        Rcpp::Named("en_cure_coef") = obj.en_cure_coef,
        Rcpp::Named("negLogL") = obj.negLogL,
        Rcpp::Named("nObs") = obj.nObs,
        Rcpp::Named("cox_l1_lambda_max") = obj.cox_l1_lambda_max,
        Rcpp::Named("cox_l1_lambda") = obj.cox_l1_lambda,
        Rcpp::Named("cox_l2_lambda") = obj.cox_l2_lambda,
        Rcpp::Named("cox_l1_penalty_factor") = obj.cox_l1_penalty_factor,
        Rcpp::Named("cure_l1_lambda_max") = obj.cure_l1_lambda_max,
        Rcpp::Named("cure_l1_lambda") = obj.cure_l1_lambda,
        Rcpp::Named("cure_l2_lambda") = obj.cure_l2_lambda,
        Rcpp::Named("cure_l1_penalty_factor") = obj.cure_l1_penalty_factor
        );
}
