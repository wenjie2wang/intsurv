#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "../include/intsurv.h"

// check the subset method for CoxphCure objects
// [[Rcpp::export]]
Rcpp::List rt_cv(const unsigned int nobs,
                 const unsigned int nfolds,
                 const arma::uvec& static_train_index,
                 const arma::uvec& strata)
{
    intsurv::CrossValidation cv_obj {
        nobs, nfolds, static_train_index, strata
    };
    Rcpp::List train_list, valid_list;
    for (size_t i {0}; i < nfolds; ++i) {
        train_list.push_back(intsurv::arma2rvec(cv_obj.train_index_.at(i)));
        valid_list.push_back(intsurv::arma2rvec(cv_obj.test_index_.at(i)));
    }
    return Rcpp::List::create(
        Rcpp::Named("train_index") = train_list,
        Rcpp::Named("valid_index") = valid_list
        );
}

