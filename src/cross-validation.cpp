#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <intsurv/CrossValidation.h>
#include <intsurv/utils.h>

// [[Rcpp::export]]
Rcpp::List rcpp_gen_cv_index(const unsigned int nobs,
                             const unsigned int nfolds,
                             const arma::uvec& strata,
                             const arma::uvec& static_train_index)
{
    intsurv::CrossValidation cv_obj {
        nobs, nfolds, static_train_index, strata
    };
    Rcpp::List train_list, valid_list;
    arma::uvec tmp;
    for (size_t i {0}; i < nfolds; ++i) {
        tmp = cv_obj.train_index_.at(i) + 1;
        train_list.push_back(intsurv::arma2rvec(tmp));
        tmp = cv_obj.test_index_.at(i) + 1;
        valid_list.push_back(intsurv::arma2rvec(tmp));
    }
    return Rcpp::List::create(
        Rcpp::Named("train_index") = train_list,
        Rcpp::Named("valid_index") = valid_list
        );
}
