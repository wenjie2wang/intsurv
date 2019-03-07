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

// rcpp_coxph
arma::vec rcpp_coxph(const arma::vec& time, const arma::vec& event, const arma::mat& x, const arma::vec& start, const unsigned int max_iter, const double rel_tol);
RcppExport SEXP _intsurv_rcpp_coxph(SEXP timeSEXP, SEXP eventSEXP, SEXP xSEXP, SEXP startSEXP, SEXP max_iterSEXP, SEXP rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type rel_tol(rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_coxph(time, event, x, start, max_iter, rel_tol));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_reg_coxph
arma::vec rcpp_reg_coxph(const arma::vec& time, const arma::vec& event, const arma::mat& x, const double lambda, arma::vec penalty_factor, const arma::vec& start, const unsigned int max_iter, const double rel_tol);
RcppExport SEXP _intsurv_rcpp_reg_coxph(SEXP timeSEXP, SEXP eventSEXP, SEXP xSEXP, SEXP lambdaSEXP, SEXP penalty_factorSEXP, SEXP startSEXP, SEXP max_iterSEXP, SEXP rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type penalty_factor(penalty_factorSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type rel_tol(rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_reg_coxph(time, event, x, lambda, penalty_factor, start, max_iter, rel_tol));
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
// rcpp_logistic
Rcpp::NumericVector rcpp_logistic(const arma::mat& x, const arma::vec& y, const arma::vec start, const unsigned int max_iter, const double rel_tol);
RcppExport SEXP _intsurv_rcpp_logistic(SEXP xSEXP, SEXP ySEXP, SEXP startSEXP, SEXP max_iterSEXP, SEXP rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type rel_tol(rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_logistic(x, y, start, max_iter, rel_tol));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_reg_logistic
Rcpp::NumericVector rcpp_reg_logistic(const arma::mat& x, const arma::vec& y, const double lambda, arma::vec penalty_factor, const arma::vec start, const unsigned int max_iter, const double rel_tol);
RcppExport SEXP _intsurv_rcpp_reg_logistic(SEXP xSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP penalty_factorSEXP, SEXP startSEXP, SEXP max_iterSEXP, SEXP rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type penalty_factor(penalty_factorSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type rel_tol(rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_reg_logistic(x, y, lambda, penalty_factor, start, max_iter, rel_tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_intsurv_rcpp_coxph", (DL_FUNC) &_intsurv_rcpp_coxph, 6},
    {"_intsurv_rcpp_reg_coxph", (DL_FUNC) &_intsurv_rcpp_reg_coxph, 8},
    {"_intsurv_aggregateSum", (DL_FUNC) &_intsurv_aggregateSum, 5},
    {"_intsurv_revcumsum", (DL_FUNC) &_intsurv_revcumsum, 1},
    {"_intsurv_rcpp_logistic", (DL_FUNC) &_intsurv_rcpp_logistic, 5},
    {"_intsurv_rcpp_reg_logistic", (DL_FUNC) &_intsurv_rcpp_reg_logistic, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_intsurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
