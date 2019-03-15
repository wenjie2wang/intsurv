//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
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

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// coxph_cure
Rcpp::List coxph_cure(const arma::vec& time, const arma::vec& event, const arma::mat& cox_x, const arma::mat& cure_x, const bool cure_intercept, const arma::vec& cox_start, const arma::vec& cure_start, const unsigned int& em_max_iter, const double& em_rel_tol, const unsigned int& cox_mstep_max_iter, const double& cox_mstep_rel_tol, const unsigned int& cure_mstep_max_iter, const double& cure_mstep_rel_tol);
RcppExport SEXP _intsurv_coxph_cure(SEXP timeSEXP, SEXP eventSEXP, SEXP cox_xSEXP, SEXP cure_xSEXP, SEXP cure_interceptSEXP, SEXP cox_startSEXP, SEXP cure_startSEXP, SEXP em_max_iterSEXP, SEXP em_rel_tolSEXP, SEXP cox_mstep_max_iterSEXP, SEXP cox_mstep_rel_tolSEXP, SEXP cure_mstep_max_iterSEXP, SEXP cure_mstep_rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cox_x(cox_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cure_x(cure_xSEXP);
    Rcpp::traits::input_parameter< const bool >::type cure_intercept(cure_interceptSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cox_start(cox_startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cure_start(cure_startSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type em_max_iter(em_max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type em_rel_tol(em_rel_tolSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type cox_mstep_max_iter(cox_mstep_max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type cox_mstep_rel_tol(cox_mstep_rel_tolSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type cure_mstep_max_iter(cure_mstep_max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type cure_mstep_rel_tol(cure_mstep_rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(coxph_cure(time, event, cox_x, cure_x, cure_intercept, cox_start, cure_start, em_max_iter, em_rel_tol, cox_mstep_max_iter, cox_mstep_rel_tol, cure_mstep_max_iter, cure_mstep_rel_tol));
    return rcpp_result_gen;
END_RCPP
}
// coxph_cure_reg
Rcpp::List coxph_cure_reg(const arma::vec& time, const arma::vec& event, const arma::mat& cox_x, const arma::mat& cure_x, const bool cure_intercept, const arma::vec& cox_lambda, const unsigned int& cox_nlambda, double cox_lambda_min_ratio, const arma::vec& cox_penalty_factor, const arma::vec& cure_lambda, const unsigned int& cure_nlambda, double cure_lambda_min_ratio, const arma::vec& cure_penalty_factor, const arma::vec& cox_start, const arma::vec& cure_start, const unsigned int& em_max_iter, const double& em_rel_tol, const unsigned int& cox_mstep_max_iter, const double& cox_mstep_rel_tol, const unsigned int& cure_mstep_max_iter, const double& cure_mstep_rel_tol);
RcppExport SEXP _intsurv_coxph_cure_reg(SEXP timeSEXP, SEXP eventSEXP, SEXP cox_xSEXP, SEXP cure_xSEXP, SEXP cure_interceptSEXP, SEXP cox_lambdaSEXP, SEXP cox_nlambdaSEXP, SEXP cox_lambda_min_ratioSEXP, SEXP cox_penalty_factorSEXP, SEXP cure_lambdaSEXP, SEXP cure_nlambdaSEXP, SEXP cure_lambda_min_ratioSEXP, SEXP cure_penalty_factorSEXP, SEXP cox_startSEXP, SEXP cure_startSEXP, SEXP em_max_iterSEXP, SEXP em_rel_tolSEXP, SEXP cox_mstep_max_iterSEXP, SEXP cox_mstep_rel_tolSEXP, SEXP cure_mstep_max_iterSEXP, SEXP cure_mstep_rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cox_x(cox_xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cure_x(cure_xSEXP);
    Rcpp::traits::input_parameter< const bool >::type cure_intercept(cure_interceptSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cox_lambda(cox_lambdaSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type cox_nlambda(cox_nlambdaSEXP);
    Rcpp::traits::input_parameter< double >::type cox_lambda_min_ratio(cox_lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cox_penalty_factor(cox_penalty_factorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cure_lambda(cure_lambdaSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type cure_nlambda(cure_nlambdaSEXP);
    Rcpp::traits::input_parameter< double >::type cure_lambda_min_ratio(cure_lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cure_penalty_factor(cure_penalty_factorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cox_start(cox_startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cure_start(cure_startSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type em_max_iter(em_max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type em_rel_tol(em_rel_tolSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type cox_mstep_max_iter(cox_mstep_max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type cox_mstep_rel_tol(cox_mstep_rel_tolSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type cure_mstep_max_iter(cure_mstep_max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type cure_mstep_rel_tol(cure_mstep_rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(coxph_cure_reg(time, event, cox_x, cure_x, cure_intercept, cox_lambda, cox_nlambda, cox_lambda_min_ratio, cox_penalty_factor, cure_lambda, cure_nlambda, cure_lambda_min_ratio, cure_penalty_factor, cox_start, cure_start, em_max_iter, em_rel_tol, cox_mstep_max_iter, cox_mstep_rel_tol, cure_mstep_max_iter, cure_mstep_rel_tol));
    return rcpp_result_gen;
END_RCPP
}
// aggregateSum
Rcpp::NumericVector aggregateSum(const arma::vec& x, const arma::vec& indices, const bool simplify, const bool cumulative, const bool reversely);
RcppExport SEXP _intsurv_aggregateSum(SEXP xSEXP, SEXP indicesSEXP, SEXP simplifySEXP, SEXP cumulativeSEXP, SEXP reverselySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const bool >::type simplify(simplifySEXP);
    Rcpp::traits::input_parameter< const bool >::type cumulative(cumulativeSEXP);
    Rcpp::traits::input_parameter< const bool >::type reversely(reverselySEXP);
    rcpp_result_gen = Rcpp::wrap(aggregateSum(x, indices, simplify, cumulative, reversely));
    return rcpp_result_gen;
END_RCPP
}
// revcumsum
Rcpp::NumericVector revcumsum(const arma::vec& x);
RcppExport SEXP _intsurv_revcumsum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(revcumsum(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_intsurv_coxph_cure", (DL_FUNC) &_intsurv_coxph_cure, 13},
    {"_intsurv_coxph_cure_reg", (DL_FUNC) &_intsurv_coxph_cure_reg, 21},
    {"_intsurv_aggregateSum", (DL_FUNC) &_intsurv_aggregateSum, 5},
    {"_intsurv_revcumsum", (DL_FUNC) &_intsurv_revcumsum, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_intsurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
