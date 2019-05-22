#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "logistic_reg.h"
#include "utils.h"


// fitting regular logistic model by monotonic quadratic approximation algorithm
// non-integer y vector is allowed
// [[Rcpp::export]]
Rcpp::List rcpp_logistic(
    const arma::mat& x,
    const arma::vec& y,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec start = 0,
    const unsigned int max_iter = 300,
    const double rel_tol = 1e-6
    )
{
    Intsurv::LogisticReg object {
        Intsurv::LogisticReg(x, y, intercept, standardize)
    };
    object.fit(start, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef),
        Rcpp::Named("negLogL") = object.negLogL,
        Rcpp::Named("nObs") = object.nObs
        );
}


// fit Firth logistic regression
// [[Rcpp::export]]
Rcpp::List rcpp_firth_logistic(
    const arma::mat& x,
    const arma::vec& y,
    const bool intercept = true,
    const bool standardize = true,
    const arma::vec start = 0,
    const unsigned int max_iter = 300,
    const double rel_tol = 1e-6
    )
{
    Intsurv::LogisticReg object {
        Intsurv::LogisticReg(x, y, intercept, standardize)
    };
    object.firth_fit(start, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef),
        Rcpp::Named("negLogL") = object.negLogL,
        Rcpp::Named("nObs") = object.nObs
        );
}


// regularized logistic model by coordinate-majorization-descent algorithm
// for perticular lambda's
// [[Rcpp::export]]
Rcpp::List rcpp_reg_logistic1(const arma::mat& x,
                              const arma::vec& y,
                              const double& l1_lambda = 0,
                              const double& l2_lambda = 0,
                              const arma::vec& l1_penalty_factor = 0,
                              const arma::vec& start = 0,
                              const bool intercept = true,
                              const unsigned int max_iter = 1000,
                              const double rel_tol = 1e-6)
{
    Intsurv::LogisticReg object { Intsurv::LogisticReg(x, y, intercept) };
    object.regularized_fit(l1_lambda, l2_lambda, l1_penalty_factor,
                           start, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = Intsurv::arma2rvec(object.coef),
        Rcpp::Named("en_coef") = Intsurv::arma2rvec(object.en_coef),
        Rcpp::Named("negLogL") = object.negLogL,
        Rcpp::Named("nObs") = object.nObs,
        Rcpp::Named("l1_lambda_max") = object.l1_lambda_max,
        Rcpp::Named("l1_lambda") = object.l1_lambda,
        Rcpp::Named("l2_lambda") = object.l2_lambda,
        Rcpp::Named("l1_penalty_factor") = object.l1_penalty_factor
        );
}


// for a sequence of lambda's and a given alpha
// [[Rcpp::export]]
Rcpp::List rcpp_reg_logistic2(const arma::mat& x,
                              const arma::vec& y,
                              arma::vec lambda = 0,
                              const double alpha = 1,
                              const unsigned int& nlambda = 1,
                              double lambda_min_ratio = 1e-4,
                              const arma::vec& penalty_factor = 0,
                              const bool intercept = true,
                              const unsigned int max_iter = 1000,
                              const double rel_tol = 1e-6)
{
    Intsurv::LogisticReg object { Intsurv::LogisticReg(x, y, intercept) };
    object.regularized_fit(lambda, alpha, nlambda, lambda_min_ratio,
                           penalty_factor, max_iter, rel_tol);
    return Rcpp::List::create(
        Rcpp::Named("coef") = object.coef_mat,
        Rcpp::Named("en_coef") = object.en_coef_mat,
        Rcpp::Named("negLogL") = Intsurv::arma2rvec(object.negLogL_vec),
        Rcpp::Named("nObs") = object.nObs,
        Rcpp::Named("lambda_max") = object.l1_lambda_max,
        Rcpp::Named("lambda") = Intsurv::arma2rvec(object.lambda_vec),
        Rcpp::Named("alpha") = object.alpha,
        Rcpp::Named("l1_penalty_factor") = object.l1_penalty_factor
        );
}
